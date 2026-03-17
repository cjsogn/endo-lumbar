#!/usr/bin/env python3
"""
Supplementary analysis: G-computation on follow-up-eligible restricted samples.

Replicates the exact same PyMC Bayesian G-computation from 22_pymc_bayesian.py
but on date-restricted samples:
  - 3m outcomes: patients with surgery_date <= 2025-09-30
  - 12m outcomes: patients with surgery_date <= 2024-12-31

Same models, same priors, same covariates, same MCMC settings.
Results stored separately in supplementary_followup_eligible/.
"""

import numpy as np
import pandas as pd
import pymc as pm
import pytensor.tensor as pt
import arviz as az
from scipy.special import expit
import os
import time
import gc
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# 1. CONFIGURATION (identical to 22_pymc_bayesian.py)
# =============================================================================

SEED = 20260204
N_CORES = 14
N_CHAINS = 4
N_DRAWS = 1000
N_TUNE = 1000
TARGET_ACCEPT = 0.95
MAX_TREEDEPTH = 12

NI_MARGINS = {
    "odi": 7.0,
    "nrs": 1.0,
    "eq5d": 0.05,
    "binary": 0.10,
}
NI_THRESHOLD = 0.95

# Paths
BASE_DIR = "/Users/cjsogn/ENDO_LUMBAR_FINAL/supplementary_followup_eligible"
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
MODEL_DIR = os.path.join(BASE_DIR, "models")

for d in [RESULTS_DIR, MODEL_DIR]:
    os.makedirs(d, exist_ok=True)

# Covariate specification (identical to main analysis)
CONTINUOUS_COVS = ["age", "bmi", "odi_baseline", "eq5d_baseline",
                   "nrs_back_baseline", "nrs_leg_baseline", "n_prior_surgeries"]
FACTOR_COVS = ["sex", "smoking", "education", "asa_cat"]
BINARY_COVS = [
    "employed_baseline", "sick_leave", "disability", "analgesic_baseline",
    "symptom_duration_back", "symptom_duration_leg", "motor_deficit",
    "depression_anxiety", "chronic_pain", "spondylolisthesis", "scoliosis",
    "prior_surgery_any", "multilevel", "prolapse_intraforaminal",
    "prolapse_extralateral", "stenosis_central",
]

# Global variable for column names (used by model functions)
col_names_global = []


# =============================================================================
# 2. DATA LOADING & PREPARATION (identical logic)
# =============================================================================

def load_and_prepare_data():
    """Load restricted CSVs and build design matrices."""
    print("=" * 80)
    print("LOADING FOLLOW-UP ELIGIBLE RESTRICTED DATA")
    print("=" * 80)

    df_3m = pd.read_csv(os.path.join(DATA_DIR, "df_3m_eligible.csv"))
    df_12m = pd.read_csv(os.path.join(DATA_DIR, "df_12m_eligible.csv"))

    print(f"3m-eligible: N={len(df_3m)} "
          f"(ELD={df_3m['treatment_num'].sum():.0f}, "
          f"MSD={(1-df_3m['treatment_num']).sum():.0f})")
    print(f"12m-eligible: N={len(df_12m)} "
          f"(ELD={df_12m['treatment_num'].sum():.0f}, "
          f"MSD={(1-df_12m['treatment_num']).sum():.0f})")

    return df_3m, df_12m


def build_design_matrix(df):
    """Build standardized design matrix X from dataframe."""
    cols = []
    col_names = []
    std_params = {}

    for var in CONTINUOUS_COVS:
        vals = df[var].values.astype(float)
        m = np.nanmean(vals)
        s = max(np.nanstd(vals), 0.01)
        std_params[var] = {"mean": float(m), "std": float(s)}
        z = (vals - m) / s
        z = np.where(np.isnan(z), 0.0, z)
        cols.append(z)
        col_names.append(var + "_z")

    for var in FACTOR_COVS:
        vals = df[var].values
        unique_vals = sorted(pd.Series(vals).dropna().unique())
        if len(unique_vals) <= 2:
            col = (vals == unique_vals[1]).astype(float)
            cols.append(col)
            col_names.append(f"{var}_{unique_vals[1]}")
        else:
            for level in unique_vals[1:]:
                col = (vals == level).astype(float)
                cols.append(col)
                col_names.append(f"{var}_{level}")

    for var in BINARY_COVS:
        vals = df[var].values.astype(float)
        vals = np.where(np.isnan(vals), 0.0, vals)
        cols.append(vals)
        col_names.append(var)

    X = np.column_stack(cols)
    print(f"Design matrix: {X.shape[0]} x {X.shape[1]} "
          f"({len(col_names)} covariates)")
    return X, col_names, std_params


# =============================================================================
# 3. ZOIB MODEL (identical to 22_pymc_bayesian.py)
# =============================================================================

def fit_zoib_model(df, X, outcome_var, scale_upper, zi_baseline_var_z,
                   outcome_label):
    """Fit custom ZOIB model using pm.Potential for manual log-likelihood."""
    print(f"\n--- ZOIB Model: {outcome_label} ---")

    N = len(df)
    treatment = df["treatment_num"].values.astype(float)
    y_raw = df[outcome_var].values.astype(float)

    y_01 = y_raw / scale_upper
    y_01 = np.where(np.isnan(y_01), np.nan, np.minimum(y_01, 1.0 - 1e-6))

    obs_mask = ~np.isnan(y_01)
    obs_idx = np.where(obs_mask)[0]
    y_obs = y_01[obs_mask]
    n_obs = len(y_obs)

    zero_mask_obs = y_obs == 0.0
    pos_mask_obs = y_obs > 0.0
    zero_idx_obs = np.where(zero_mask_obs)[0]
    pos_idx_obs = np.where(pos_mask_obs)[0]
    y_pos = y_obs[pos_mask_obs]

    n_zeros = zero_mask_obs.sum()
    n_pos = pos_mask_obs.sum()

    print(f"  N={N}, Observed={n_obs} ({100*n_obs/N:.1f}%), "
          f"Zeros={n_zeros} ({100*n_zeros/n_obs:.1f}%), "
          f"Positive={n_pos}")

    zi_base_col = None
    for i, name in enumerate(col_names_global):
        if name == zi_baseline_var_z:
            zi_base_col = i
            break
    if zi_base_col is None:
        raise ValueError(f"ZI baseline variable {zi_baseline_var_z} not found "
                         f"in design matrix columns")

    n_covs = X.shape[1]

    with pm.Model() as model:
        intercept_mu = pm.Normal("intercept_mu", mu=0, sigma=3)
        beta_treatment = pm.Normal("beta_treatment", mu=0, sigma=1)
        beta_covs = pm.Normal("beta_covs", mu=0, sigma=0.5, shape=n_covs)

        eta_mu = (intercept_mu
                  + beta_treatment * pt.as_tensor_variable(treatment)
                  + pt.dot(pt.as_tensor_variable(X), beta_covs))
        mu_all = pm.math.sigmoid(eta_mu)

        intercept_zi = pm.Normal("intercept_zi", mu=0, sigma=1.5)
        beta_zi_treatment = pm.Normal("beta_zi_treatment", mu=0, sigma=1)
        beta_zi_baseline = pm.Normal("beta_zi_baseline", mu=0, sigma=1)

        zi_baseline_vals = pt.as_tensor_variable(X[:, zi_base_col])
        eta_zi = (intercept_zi
                  + beta_zi_treatment * pt.as_tensor_variable(treatment)
                  + beta_zi_baseline * zi_baseline_vals)
        p_zero_all = pm.math.sigmoid(eta_zi)

        phi = pm.Gamma("phi", alpha=2, beta=0.1)

        mu_obs = mu_all[obs_idx]
        p_zero_obs = p_zero_all[obs_idx]

        alpha_beta = mu_obs * phi
        beta_beta = (1.0 - mu_obs) * phi

        logp_zeros = pt.sum(pt.log(p_zero_obs[zero_idx_obs] + 1e-12))
        logp_pos_zi = pt.sum(pt.log(1.0 - p_zero_obs[pos_idx_obs] + 1e-12))

        alpha_pos = alpha_beta[pos_idx_obs]
        beta_pos = beta_beta[pos_idx_obs]
        y_pos_tensor = pt.as_tensor_variable(y_pos)

        logp_beta = pt.sum(
            pt.gammaln(alpha_pos + beta_pos)
            - pt.gammaln(alpha_pos)
            - pt.gammaln(beta_pos)
            + (alpha_pos - 1.0) * pt.log(y_pos_tensor + 1e-12)
            + (beta_pos - 1.0) * pt.log(1.0 - y_pos_tensor + 1e-12)
        )

        pm.Potential("zoib_logp", logp_zeros + logp_pos_zi + logp_beta)

        idata = pm.sample(
            draws=N_DRAWS, tune=N_TUNE, chains=N_CHAINS,
            target_accept=TARGET_ACCEPT, cores=N_CORES,
            random_seed=SEED, return_inferencedata=True, progressbar=True,
            nuts_sampler_kwargs={"max_treedepth": MAX_TREEDEPTH},
        )

    diag = check_diagnostics(idata, outcome_label)

    # G-computation
    posterior = idata.posterior
    S = N_CHAINS * N_DRAWS

    intercept_mu_s = posterior["intercept_mu"].values.flatten()
    beta_t_s = posterior["beta_treatment"].values.flatten()
    beta_covs_s = posterior["beta_covs"].values.reshape(S, -1)

    intercept_zi_s = posterior["intercept_zi"].values.flatten()
    beta_zi_t_s = posterior["beta_zi_treatment"].values.flatten()
    beta_zi_b_s = posterior["beta_zi_baseline"].values.flatten()

    cov_contrib = beta_covs_s @ X.T
    zi_base_vals_np = X[:, zi_base_col]

    eta_mu_eld = intercept_mu_s[:, None] + beta_t_s[:, None] * 1.0 + cov_contrib
    mu_eld = expit(eta_mu_eld)
    eta_zi_eld = (intercept_zi_s[:, None]
                  + beta_zi_t_s[:, None] * 1.0
                  + beta_zi_b_s[:, None] * zi_base_vals_np[None, :])
    p_zero_eld = expit(eta_zi_eld)
    ey_eld = (1.0 - p_zero_eld) * mu_eld * scale_upper

    eta_mu_msd = intercept_mu_s[:, None] + beta_t_s[:, None] * 0.0 + cov_contrib
    mu_msd = expit(eta_mu_msd)
    eta_zi_msd = (intercept_zi_s[:, None]
                  + beta_zi_t_s[:, None] * 0.0
                  + beta_zi_b_s[:, None] * zi_base_vals_np[None, :])
    p_zero_msd = expit(eta_zi_msd)
    ey_msd = (1.0 - p_zero_msd) * mu_msd * scale_upper

    ate_samples = ey_msd.mean(axis=1) - ey_eld.mean(axis=1)

    return idata, ate_samples, diag, N, n_obs


# =============================================================================
# 4. GAUSSIAN MODEL (identical to 22_pymc_bayesian.py)
# =============================================================================

def fit_gaussian_model(df, X, outcome_var, outcome_label,
                       prior_intercept_mu=30, prior_intercept_sigma=20,
                       prior_treatment_sigma=10, prior_sigma_scale=15,
                       prior_odi_baseline_mu=1):
    """Fit Gaussian model with masked array for missing outcome handling."""
    print(f"\n--- Gaussian Model: {outcome_label} ---")

    N = len(df)
    treatment = df["treatment_num"].values.astype(float)
    y = df[outcome_var].values.astype(float)

    obs_mask = ~np.isnan(y)
    n_obs = obs_mask.sum()
    print(f"  N={N}, Observed={n_obs} ({100*n_obs/N:.1f}%)")

    odi_base_col = None
    for i, name in enumerate(col_names_global):
        if name == "odi_baseline_z":
            odi_base_col = i
            break

    n_covs = X.shape[1]
    y_masked = np.ma.masked_invalid(y)

    with pm.Model() as model:
        intercept = pm.Normal("intercept", mu=prior_intercept_mu,
                              sigma=prior_intercept_sigma)
        beta_treatment = pm.Normal("beta_treatment", mu=0,
                                   sigma=prior_treatment_sigma)

        if odi_base_col is not None:
            beta_odi_baseline = pm.Normal("beta_odi_baseline",
                                          mu=prior_odi_baseline_mu, sigma=0.5)
            other_idx = [i for i in range(n_covs) if i != odi_base_col]
            beta_others = pm.Normal("beta_others", mu=0, sigma=2,
                                    shape=len(other_idx))
            X_others = X[:, other_idx]
            odi_z = X[:, odi_base_col]

            mu = (intercept
                  + beta_treatment * pt.as_tensor_variable(treatment)
                  + beta_odi_baseline * pt.as_tensor_variable(odi_z)
                  + pt.dot(pt.as_tensor_variable(X_others), beta_others))
        else:
            beta_others = pm.Normal("beta_others", mu=0, sigma=2, shape=n_covs)
            mu = (intercept
                  + beta_treatment * pt.as_tensor_variable(treatment)
                  + pt.dot(pt.as_tensor_variable(X), beta_others))

        sigma = pm.HalfStudentT("sigma", nu=3,
                                lam=1.0 / (prior_sigma_scale**2))

        pm.Normal("y_obs", mu=mu, sigma=sigma, observed=y_masked)

        idata = pm.sample(
            draws=N_DRAWS, tune=N_TUNE, chains=N_CHAINS,
            target_accept=TARGET_ACCEPT, cores=N_CORES,
            random_seed=SEED, return_inferencedata=True, progressbar=True,
            nuts_sampler_kwargs={"max_treedepth": MAX_TREEDEPTH},
        )

    diag = check_diagnostics(idata, outcome_label)

    # G-computation
    posterior = idata.posterior
    S = N_CHAINS * N_DRAWS

    intercept_s = posterior["intercept"].values.flatten()
    beta_t_s = posterior["beta_treatment"].values.flatten()

    if odi_base_col is not None:
        beta_odi_s = posterior["beta_odi_baseline"].values.flatten()
        beta_oth_s = posterior["beta_others"].values.reshape(S, -1)
        odi_z_np = X[:, odi_base_col]
        X_others_np = X[:, other_idx]

        cov_contrib = (beta_odi_s[:, None] * odi_z_np[None, :]
                       + beta_oth_s @ X_others_np.T)
    else:
        beta_oth_s = posterior["beta_others"].values.reshape(S, -1)
        cov_contrib = beta_oth_s @ X.T

    mu_eld = intercept_s[:, None] + beta_t_s[:, None] * 1.0 + cov_contrib
    mu_msd = intercept_s[:, None] + beta_t_s[:, None] * 0.0 + cov_contrib

    # Higher-is-better (EQ-5D)
    ate_samples = mu_eld.mean(axis=1) - mu_msd.mean(axis=1)

    return idata, ate_samples, diag, N, n_obs


# =============================================================================
# 5. BERNOULLI MODEL (identical to 22_pymc_bayesian.py)
# =============================================================================

def fit_bernoulli_model(df, X, outcome_var, outcome_label, lower_is_better):
    """Fit Bernoulli model with logit link and observed-case likelihood."""
    print(f"\n--- Bernoulli Model: {outcome_label} ---")

    N = len(df)
    treatment = df["treatment_num"].values.astype(float)
    y = df[outcome_var].values.astype(float)

    obs_mask = ~np.isnan(y)
    obs_idx = np.where(obs_mask)[0]
    y_obs = y[obs_mask]
    n_obs = len(y_obs)
    print(f"  N={N}, Observed={n_obs} ({100*n_obs/N:.1f}%)")

    n_covs = X.shape[1]

    with pm.Model() as model:
        intercept = pm.Normal("intercept", mu=0, sigma=5)
        beta_treatment = pm.Normal("beta_treatment", mu=0, sigma=1)
        beta_covs = pm.Normal("beta_covs", mu=0, sigma=0.5,
                              shape=n_covs)

        logit_p = (intercept
                   + beta_treatment * pt.as_tensor_variable(treatment)
                   + pt.dot(pt.as_tensor_variable(X), beta_covs))

        pm.Bernoulli("y_obs", logit_p=logit_p[obs_idx], observed=y_obs)

        idata = pm.sample(
            draws=N_DRAWS, tune=N_TUNE, chains=N_CHAINS,
            target_accept=TARGET_ACCEPT, cores=N_CORES,
            random_seed=SEED, return_inferencedata=True, progressbar=True,
            nuts_sampler_kwargs={"max_treedepth": MAX_TREEDEPTH},
        )

    diag = check_diagnostics(idata, outcome_label)

    # G-computation
    posterior = idata.posterior
    S = N_CHAINS * N_DRAWS

    intercept_s = posterior["intercept"].values.flatten()
    beta_t_s = posterior["beta_treatment"].values.flatten()
    beta_c_s = posterior["beta_covs"].values.reshape(S, -1)

    cov_contrib = beta_c_s @ X.T

    p_eld = expit(intercept_s[:, None] + beta_t_s[:, None] * 1.0 + cov_contrib)
    p_msd = expit(intercept_s[:, None] + beta_t_s[:, None] * 0.0 + cov_contrib)

    if lower_is_better:
        ate_samples = p_msd.mean(axis=1) - p_eld.mean(axis=1)
    else:
        ate_samples = p_eld.mean(axis=1) - p_msd.mean(axis=1)

    return idata, ate_samples, diag, N, n_obs


# =============================================================================
# 6. CONVERGENCE DIAGNOSTICS (identical)
# =============================================================================

def check_diagnostics(idata, label):
    """Check MCMC diagnostics: R-hat, ESS, divergences."""
    rhat = az.rhat(idata)
    max_rhat = max(float(rhat[var].values.max()) for var in rhat.data_vars)

    ess_bulk = az.ess(idata, method="bulk")
    ess_tail = az.ess(idata, method="tail")
    min_bulk = min(float(ess_bulk[var].values.min())
                   for var in ess_bulk.data_vars)
    min_tail = min(float(ess_tail[var].values.min())
                   for var in ess_tail.data_vars)

    n_divergent = 0
    if hasattr(idata, "sample_stats") and "diverging" in idata.sample_stats:
        n_divergent = int(idata.sample_stats["diverging"].sum().values)

    valid = (max_rhat < 1.01 and min_bulk > 400 and min_tail > 400
             and n_divergent == 0)

    status = "PASS" if valid else "ISSUE"
    print(f"  Convergence: {status} (Rhat_max={max_rhat:.4f}, "
          f"ESS_bulk_min={min_bulk:.0f}, ESS_tail_min={min_tail:.0f}, "
          f"Divergences={n_divergent})")

    return {
        "outcome": label,
        "rhat_max": float(max_rhat),
        "ess_bulk_min": float(min_bulk),
        "ess_tail_min": float(min_tail),
        "n_divergent": n_divergent,
        "valid": valid,
    }


# =============================================================================
# 7. OUTCOME DEFINITIONS (identical specs, different dataset labels)
# =============================================================================

def get_outcome_specs():
    """Define all 17 outcomes with their specifications.

    Key difference from main: 3m outcomes use '3m_eligible' dataset,
    12m outcomes use '12m_eligible' dataset.
    """
    outcomes = [
        # ZOIB outcomes (lower-is-better)
        {"var": "odi_3m", "label": "ODI 3 months", "family": "zoib",
         "scale": 100, "zi_baseline_z": "odi_baseline_z",
         "dataset": "3m", "margin": NI_MARGINS["odi"],
         "lower_is_better": True},
        {"var": "odi_12m", "label": "ODI 12 months", "family": "zoib",
         "scale": 100, "zi_baseline_z": "odi_baseline_z",
         "dataset": "12m", "margin": NI_MARGINS["odi"],
         "lower_is_better": True},
        {"var": "nrs_back_3m", "label": "NRS back 3 months", "family": "zoib",
         "scale": 10, "zi_baseline_z": "nrs_back_baseline_z",
         "dataset": "3m", "margin": NI_MARGINS["nrs"],
         "lower_is_better": True},
        {"var": "nrs_back_12m", "label": "NRS back 12 months", "family": "zoib",
         "scale": 10, "zi_baseline_z": "nrs_back_baseline_z",
         "dataset": "12m", "margin": NI_MARGINS["nrs"],
         "lower_is_better": True},
        {"var": "nrs_leg_3m", "label": "NRS leg 3 months", "family": "zoib",
         "scale": 10, "zi_baseline_z": "nrs_leg_baseline_z",
         "dataset": "3m", "margin": NI_MARGINS["nrs"],
         "lower_is_better": True},
        {"var": "nrs_leg_12m", "label": "NRS leg 12 months", "family": "zoib",
         "scale": 10, "zi_baseline_z": "nrs_leg_baseline_z",
         "dataset": "12m", "margin": NI_MARGINS["nrs"],
         "lower_is_better": True},
        # Gaussian outcomes (higher-is-better)
        {"var": "eq5d_3m", "label": "EQ-5D 3 months", "family": "gaussian",
         "dataset": "3m", "margin": NI_MARGINS["eq5d"],
         "lower_is_better": False,
         "prior_intercept_mu": 0.5, "prior_intercept_sigma": 0.5,
         "prior_treatment_sigma": 0.5, "prior_sigma_scale": 0.3,
         "prior_odi_baseline_mu": 0},
        {"var": "eq5d_12m", "label": "EQ-5D 12 months", "family": "gaussian",
         "dataset": "12m", "margin": NI_MARGINS["eq5d"],
         "lower_is_better": False,
         "prior_intercept_mu": 0.5, "prior_intercept_sigma": 0.5,
         "prior_treatment_sigma": 0.5, "prior_sigma_scale": 0.3,
         "prior_odi_baseline_mu": 0},
        # Bernoulli outcomes
        {"var": "responder_3m", "label": "Responder 3 months",
         "family": "bernoulli", "dataset": "3m",
         "margin": NI_MARGINS["binary"], "lower_is_better": False},
        {"var": "rtw_3m", "label": "RTW 3 months",
         "family": "bernoulli", "dataset": "3m",
         "margin": NI_MARGINS["binary"], "lower_is_better": False},
        {"var": "rtw_12m", "label": "RTW 12 months",
         "family": "bernoulli", "dataset": "12m",
         "margin": NI_MARGINS["binary"], "lower_is_better": False},
        {"var": "analgesic_3m", "label": "Analgesic 3 months",
         "family": "bernoulli", "dataset": "3m",
         "margin": NI_MARGINS["binary"], "lower_is_better": True},
        {"var": "analgesic_12m", "label": "Analgesic 12 months",
         "family": "bernoulli", "dataset": "12m",
         "margin": NI_MARGINS["binary"], "lower_is_better": True},
        {"var": "satisfied_3m", "label": "Satisfaction 3 months",
         "family": "bernoulli", "dataset": "3m",
         "margin": NI_MARGINS["binary"], "lower_is_better": False},
        {"var": "satisfied_12m", "label": "Satisfaction 12 months",
         "family": "bernoulli", "dataset": "12m",
         "margin": NI_MARGINS["binary"], "lower_is_better": False},
        {"var": "gpe_success_3m", "label": "GPE 3 months",
         "family": "bernoulli", "dataset": "3m",
         "margin": NI_MARGINS["binary"], "lower_is_better": False},
        {"var": "gpe_success_12m", "label": "GPE 12 months",
         "family": "bernoulli", "dataset": "12m",
         "margin": NI_MARGINS["binary"], "lower_is_better": False},
    ]
    return outcomes


# =============================================================================
# 8. MAIN
# =============================================================================

if __name__ == "__main__":
    total_start = time.time()

    # Load data
    df_3m, df_12m = load_and_prepare_data()

    # Build design matrices (separate for each sample)
    print("\nBuilding design matrix for 3m-eligible sample...")
    X_3m, col_names_3m, std_params_3m = build_design_matrix(df_3m)
    print("\nBuilding design matrix for 12m-eligible sample...")
    X_12m, col_names_12m, std_params_12m = build_design_matrix(df_12m)

    # Fit all 17 outcomes
    outcomes = get_outcome_specs()
    all_results = []
    all_diags = []
    all_ate_samples = {}

    for i, spec in enumerate(outcomes):
        print(f"\n{'='*80}")
        print(f"OUTCOME {i+1}/17: {spec['label']} ({spec['family']})")
        print(f"{'='*80}")
        t0 = time.time()

        # Select dataset
        if spec["dataset"] == "12m":
            df = df_12m
            X = X_12m
            col_names_global = col_names_12m
        else:
            df = df_3m
            X = X_3m
            col_names_global = col_names_3m

        try:
            if spec["family"] == "zoib":
                idata, ate, diag, N, n_obs = fit_zoib_model(
                    df, X, spec["var"], spec["scale"],
                    spec["zi_baseline_z"], spec["label"])
            elif spec["family"] == "gaussian":
                gauss_kwargs = {k: spec[k] for k in
                    ["prior_intercept_mu", "prior_intercept_sigma",
                     "prior_treatment_sigma", "prior_sigma_scale",
                     "prior_odi_baseline_mu"] if k in spec}
                idata, ate, diag, N, n_obs = fit_gaussian_model(
                    df, X, spec["var"], spec["label"], **gauss_kwargs)
            elif spec["family"] == "bernoulli":
                idata, ate, diag, N, n_obs = fit_bernoulli_model(
                    df, X, spec["var"], spec["label"],
                    spec["lower_is_better"])

            # Summarize ATE
            ate_mean = float(ate.mean())
            ate_ci = [float(np.percentile(ate, 2.5)),
                      float(np.percentile(ate, 97.5))]
            p_ni = float((ate > -spec["margin"]).mean())
            p_sup = float((ate > 0).mean())

            ni_conclusion = "NI demonstrated" if p_ni > NI_THRESHOLD else (
                "Inconclusive" if p_ni >= 0.80 else "Concern")

            print(f"\n  ATE = {ate_mean:.4f} "
                  f"(95% CrI: [{ate_ci[0]:.4f}, {ate_ci[1]:.4f}])")
            print(f"  P(NI) = {p_ni:.4f}, P(Superiority) = {p_sup:.4f}")
            print(f"  Conclusion: {ni_conclusion}")

            elapsed = time.time() - t0
            print(f"  Time: {elapsed/60:.1f} minutes")

            result = {
                "Outcome": spec["label"],
                "Family": spec["family"].upper(),
                "Dataset": spec["dataset"],
                "N_total": N,
                "N_observed": n_obs,
                "N_gcomp": N,
                "ATE": round(ate_mean, 4),
                "CrI_lo": round(ate_ci[0], 4),
                "CrI_hi": round(ate_ci[1], 4),
                "P_NI": round(p_ni, 4),
                "P_Superiority": round(p_sup, 4),
                "NI_Margin": spec["margin"],
                "NI_Conclusion": ni_conclusion,
                "Convergence": "OK" if diag["valid"] else "Issue",
            }
            all_results.append(result)
            all_diags.append(diag)
            all_ate_samples[spec["var"]] = ate

            # Save ATE samples
            np.save(os.path.join(MODEL_DIR,
                                 f"ate_{spec['var']}.npy"), ate)

            # Free memory (don't save full idata to save space)
            del idata
            gc.collect()

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            all_results.append({
                "Outcome": spec["label"],
                "Family": spec["family"].upper(),
                "Dataset": spec["dataset"],
                "N_total": len(df),
                "N_observed": int((~np.isnan(
                    df[spec["var"]].values.astype(float))).sum()),
                "N_gcomp": len(df),
                "ATE": np.nan, "CrI_lo": np.nan, "CrI_hi": np.nan,
                "P_NI": np.nan, "P_Superiority": np.nan,
                "NI_Margin": spec["margin"],
                "NI_Conclusion": "FAILED",
                "Convergence": "FAILED",
            })
            all_diags.append({
                "outcome": spec["label"],
                "rhat_max": np.nan, "ess_bulk_min": np.nan,
                "ess_tail_min": np.nan, "n_divergent": np.nan,
                "valid": False,
            })

    # ==========================================================================
    # SAVE RESULTS
    # ==========================================================================
    print(f"\n{'='*80}")
    print("SAVING RESULTS")
    print(f"{'='*80}")

    # Results table
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(os.path.join(RESULTS_DIR,
                                    "restricted_gcomp_results.csv"),
                      index=False)
    print(f"\nResults saved to: {RESULTS_DIR}/restricted_gcomp_results.csv")

    # Convergence table
    diag_df = pd.DataFrame(all_diags)
    diag_df.to_csv(os.path.join(RESULTS_DIR,
                                 "restricted_convergence.csv"),
                   index=False)

    # Load main analysis results for comparison
    main_results_path = os.path.join(
        "/Users/cjsogn/ENDO_LUMBAR_FINAL/tables",
        "table_pymc_results.csv")
    if os.path.exists(main_results_path):
        main_df = pd.read_csv(main_results_path)
        # Build comparison table
        comparison = []
        for _, row in results_df.iterrows():
            main_row = main_df[main_df["Outcome"] == row["Outcome"]]
            if len(main_row) > 0:
                mr = main_row.iloc[0]
                comparison.append({
                    "Outcome": row["Outcome"],
                    "Main_N": mr.get("N_total", mr.get("N_gcomp", "")),
                    "Restricted_N": row["N_total"],
                    "Main_ATE": mr["ATE"],
                    "Restricted_ATE": row["ATE"],
                    "Main_CrI": f"[{mr['CrI_lo']:.2f}, {mr['CrI_hi']:.2f}]",
                    "Restricted_CrI": f"[{row['CrI_lo']:.2f}, {row['CrI_hi']:.2f}]",
                    "Main_P_NI": mr["P_NI"],
                    "Restricted_P_NI": row["P_NI"],
                    "Main_P_Sup": mr["P_Superiority"],
                    "Restricted_P_Sup": row["P_Superiority"],
                })
        comp_df = pd.DataFrame(comparison)
        comp_df.to_csv(os.path.join(RESULTS_DIR,
                                     "comparison_main_vs_restricted.csv"),
                       index=False)
        print(f"Comparison saved to: {RESULTS_DIR}/comparison_main_vs_restricted.csv")

        # Print comparison summary
        print("\n" + "=" * 100)
        print("COMPARISON: MAIN ANALYSIS vs FOLLOW-UP ELIGIBLE RESTRICTED")
        print("=" * 100)
        print(f"{'Outcome':<25} {'Main N':>7} {'Restr N':>7} "
              f"{'Main ATE':>10} {'Restr ATE':>10} "
              f"{'Main P(NI)':>10} {'Restr P(NI)':>10}")
        print("-" * 100)
        for c in comparison:
            print(f"{c['Outcome']:<25} {c['Main_N']:>7} {c['Restricted_N']:>7} "
                  f"{c['Main_ATE']:>10.4f} {c['Restricted_ATE']:>10.4f} "
                  f"{c['Main_P_NI']:>10.4f} {c['Restricted_P_NI']:>10.4f}")

    total_elapsed = time.time() - total_start
    print(f"\nTotal runtime: {total_elapsed/60:.1f} minutes")
    print("Done.")
