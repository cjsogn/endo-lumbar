#!/usr/bin/env python3
"""
22_pymc_bayesian.py - Supplementary PyMC Bayesian Analysis
ENDO-LUMBAR: Endoscopic vs microsurgical lumbar discectomy

Fits all 17 outcomes using PyMC with native missing data handling:
  - ZOIB (custom via pm.Potential): ODI 3m/12m, NRS back/leg 3m/12m
  - Gaussian (masked array): EQ-5D 3m/12m
  - Bernoulli (observed-case): All binary outcomes

Key advantage over brms: PyMC computes G-computation ATE over ALL patients
(not just complete cases), providing valid MAR inference for ZOIB/Bernoulli
families where brms cannot use mi().
"""

import numpy as np
import pandas as pd
import pymc as pm
import pytensor.tensor as pt
import arviz as az
from scipy.special import expit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os
import time
import pickle
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

SEED = 20260204
N_CORES = 14
N_CHAINS = 4
N_DRAWS = 1000
N_TUNE = 1000
TARGET_ACCEPT = 0.95
MAX_TREEDEPTH = 12

# Non-inferiority margins
NI_MARGINS = {
    "odi": 7.0,
    "nrs": 1.0,
    "eq5d": 0.05,
    "binary": 0.10,
}
NI_THRESHOLD = 0.95

# Paths
BASE_DIR = "/Users/cjsogn/ENDO_LUMBAR"
DATA_DIR = os.path.join(BASE_DIR, "data")
TABLE_DIR = os.path.join(BASE_DIR, "tables")
FIG_DIR = os.path.join(BASE_DIR, "figures")
MODEL_DIR = os.path.join(BASE_DIR, "models")

# Legacy paths for copying outputs
LEGACY_TABLE_DIR = "/Users/cjsogn/endo_studies/lumbar/analysis/output/tables"
LEGACY_FIG_DIR = "/Users/cjsogn/endo_studies/lumbar/analysis/output/figures"

for d in [TABLE_DIR, FIG_DIR, MODEL_DIR]:
    os.makedirs(d, exist_ok=True)

# Covariate specification
CONTINUOUS_COVS = ["age", "bmi", "odi_baseline", "eq5d_baseline",
                   "nrs_back_baseline", "nrs_leg_baseline", "n_prior_surgeries"]
FACTOR_COVS = ["sex", "smoking", "education", "asa_cat"]
BINARY_COVS = [
    "employed_baseline", "sick_leave", "disability", "analgesic_baseline",
    "motor_deficit",
    "depression_anxiety", "chronic_pain", "spondylolisthesis", "scoliosis",
    "prior_surgery_any", "multilevel", "prolapse_intraforaminal",
    "prolapse_extralateral", "stenosis_central",
]
# Ordinal covariates (integer-coded 1-5, treated as linear)
ORDINAL_COVS = ["symptom_duration_back", "symptom_duration_leg"]

# Publication plot settings
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 10,
    "axes.linewidth": 0.8,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 300,
})

# =============================================================================
# 2. DATA LOADING & PREPARATION
# =============================================================================

def load_and_prepare_data():
    """Load CSVs and build design matrices."""
    print("=" * 80)
    print("LOADING DATA")
    print("=" * 80)

    df_full = pd.read_csv(os.path.join(DATA_DIR, "df_disc_full.csv"))
    df_12m = pd.read_csv(os.path.join(DATA_DIR, "df_disc_12m.csv"))

    print(f"Full sample: N={len(df_full)} "
          f"(ELD={df_full['treatment_num'].sum():.0f}, "
          f"MSD={(1-df_full['treatment_num']).sum():.0f})")
    print(f"12m-eligible: N={len(df_12m)} "
          f"(ELD={df_12m['treatment_num'].sum():.0f}, "
          f"MSD={(1-df_12m['treatment_num']).sum():.0f})")

    return df_full, df_12m


def build_design_matrix(df):
    """Build standardized design matrix X from dataframe.

    Returns X (ndarray), column names list, and standardization params dict.
    """
    cols = []
    col_names = []
    std_params = {}

    # Continuous covariates: z-score
    for var in CONTINUOUS_COVS:
        vals = df[var].values.astype(float)
        m = np.nanmean(vals)
        s = max(np.nanstd(vals), 0.01)
        std_params[var] = {"mean": float(m), "std": float(s)}
        z = (vals - m) / s
        z = np.where(np.isnan(z), 0.0, z)  # covariates pre-imputed, but safety
        cols.append(z)
        col_names.append(var + "_z")

    # Factor covariates: one-hot encode (drop first level as reference)
    for var in FACTOR_COVS:
        vals = df[var].values
        unique_vals = sorted(pd.Series(vals).dropna().unique())
        if len(unique_vals) <= 2:
            # Binary factor: encode second level as 1 (first level = reference)
            ref_level = unique_vals[0]
            col = (vals == unique_vals[1]).astype(float)
            cols.append(col)
            col_names.append(f"{var}_{unique_vals[1]}")
        else:
            # Multi-level: drop first level (reference)
            for level in unique_vals[1:]:
                col = (vals == level).astype(float)
                cols.append(col)
                col_names.append(f"{var}_{level}")

    # Binary covariates
    for var in BINARY_COVS:
        vals = df[var].values.astype(float)
        vals = np.where(np.isnan(vals), 0.0, vals)
        cols.append(vals)
        col_names.append(var)

    # Ordinal covariates (integer-coded, treated as linear)
    for var in ORDINAL_COVS:
        vals = df[var].values.astype(float)
        vals = np.where(np.isnan(vals), 0.0, vals)
        cols.append(vals)
        col_names.append(var)

    X = np.column_stack(cols)
    print(f"Design matrix: {X.shape[0]} x {X.shape[1]} "
          f"({len(col_names)} covariates)")
    return X, col_names, std_params


# =============================================================================
# 3. ZOIB MODEL (ODI, NRS outcomes)
# =============================================================================

def fit_zoib_model(df, X, outcome_var, scale_upper, zi_baseline_var_z,
                   outcome_label):
    """Fit custom ZOIB model using pm.Potential for manual log-likelihood.

    Parameters
    ----------
    df : DataFrame
    X : ndarray, design matrix (all covariates)
    outcome_var : str, column name of outcome
    scale_upper : float, upper bound of scale (100 for ODI, 10 for NRS)
    zi_baseline_var_z : str, z-scored baseline var name for ZI submodel
    outcome_label : str, display name
    """
    print(f"\n--- ZOIB Model: {outcome_label} ---")

    N = len(df)
    treatment = df["treatment_num"].values.astype(float)
    y_raw = df[outcome_var].values.astype(float)

    # Transform to [0, 1): divide by upper, clamp at 1 - 1e-6, preserve zeros
    y_01 = y_raw / scale_upper
    y_01 = np.where(np.isnan(y_01), np.nan, np.minimum(y_01, 1.0 - 1e-6))

    # Observed-case indices
    obs_mask = ~np.isnan(y_01)
    obs_idx = np.where(obs_mask)[0]
    y_obs = y_01[obs_mask]
    n_obs = len(y_obs)

    # Separate zeros and positives among observed
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

    # Find ZI baseline variable index in design matrix
    # The zi_baseline_var_z is already in the design matrix
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
        # --- Mu submodel (logit link) ---
        intercept_mu = pm.Normal("intercept_mu", mu=0, sigma=3)
        beta_treatment = pm.Normal("beta_treatment", mu=0, sigma=1)
        beta_covs = pm.Normal("beta_covs", mu=0, sigma=0.5, shape=n_covs)

        # Linear predictor for mu (all N patients)
        eta_mu = (intercept_mu
                  + beta_treatment * pt.as_tensor_variable(treatment)
                  + pt.dot(pt.as_tensor_variable(X), beta_covs))
        mu_all = pm.math.sigmoid(eta_mu)

        # --- ZI submodel (logit link) ---
        intercept_zi = pm.Normal("intercept_zi", mu=0, sigma=1.5)
        beta_zi_treatment = pm.Normal("beta_zi_treatment", mu=0, sigma=1)
        beta_zi_baseline = pm.Normal("beta_zi_baseline", mu=0, sigma=1)

        zi_baseline_vals = pt.as_tensor_variable(X[:, zi_base_col])
        eta_zi = (intercept_zi
                  + beta_zi_treatment * pt.as_tensor_variable(treatment)
                  + beta_zi_baseline * zi_baseline_vals)
        p_zero_all = pm.math.sigmoid(eta_zi)

        # --- Precision ---
        phi = pm.Gamma("phi", alpha=2, beta=0.1)

        # --- Manual log-likelihood via Potential ---
        # Only observed values contribute to the likelihood.
        mu_obs = mu_all[obs_idx]
        p_zero_obs = p_zero_all[obs_idx]

        # Beta distribution parameters
        alpha_beta = mu_obs * phi
        beta_beta = (1.0 - mu_obs) * phi

        # Log-likelihood for zeros: log(p_zero)
        logp_zeros = pt.sum(pt.log(p_zero_obs[zero_idx_obs] + 1e-12))

        # Log-likelihood for positive values: log(1 - p_zero) + Beta_logp
        logp_pos_zi = pt.sum(pt.log(1.0 - p_zero_obs[pos_idx_obs] + 1e-12))

        # Beta log-pdf for positive observations
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

        # Sample
        idata = pm.sample(
            draws=N_DRAWS, tune=N_TUNE, chains=N_CHAINS,
            target_accept=TARGET_ACCEPT, cores=N_CORES,
            random_seed=SEED, return_inferencedata=True, progressbar=True,
            nuts_sampler_kwargs={"max_treedepth": MAX_TREEDEPTH},
        )

    # --- Diagnostics ---
    diag = check_diagnostics(idata, outcome_label)

    # --- G-computation ---
    # For each posterior draw: compute expected outcome for all N patients
    # E[Y | X, A] = (1 - p_zero) * mu * scale_upper
    posterior = idata.posterior
    S = N_CHAINS * N_DRAWS

    intercept_mu_s = posterior["intercept_mu"].values.flatten()
    beta_t_s = posterior["beta_treatment"].values.flatten()
    beta_covs_s = posterior["beta_covs"].values.reshape(S, -1)

    intercept_zi_s = posterior["intercept_zi"].values.flatten()
    beta_zi_t_s = posterior["beta_zi_treatment"].values.flatten()
    beta_zi_b_s = posterior["beta_zi_baseline"].values.flatten()

    # Covariate contributions (same for both counterfactuals)
    cov_contrib = beta_covs_s @ X.T  # (S, N)
    zi_base_vals_np = X[:, zi_base_col]

    # Under ELD (treatment=1)
    eta_mu_eld = intercept_mu_s[:, None] + beta_t_s[:, None] * 1.0 + cov_contrib
    mu_eld = expit(eta_mu_eld)
    eta_zi_eld = (intercept_zi_s[:, None]
                  + beta_zi_t_s[:, None] * 1.0
                  + beta_zi_b_s[:, None] * zi_base_vals_np[None, :])
    p_zero_eld = expit(eta_zi_eld)
    ey_eld = (1.0 - p_zero_eld) * mu_eld * scale_upper  # (S, N)

    # Under MSD (treatment=0)
    eta_mu_msd = intercept_mu_s[:, None] + beta_t_s[:, None] * 0.0 + cov_contrib
    mu_msd = expit(eta_mu_msd)
    eta_zi_msd = (intercept_zi_s[:, None]
                  + beta_zi_t_s[:, None] * 0.0
                  + beta_zi_b_s[:, None] * zi_base_vals_np[None, :])
    p_zero_msd = expit(eta_zi_msd)
    ey_msd = (1.0 - p_zero_msd) * mu_msd * scale_upper  # (S, N)

    # ATE: lower-is-better, positive = ELD superior
    ate_samples = ey_msd.mean(axis=1) - ey_eld.mean(axis=1)

    return idata, ate_samples, diag, N, n_obs


# =============================================================================
# 4. GAUSSIAN MODEL (EQ-5D outcomes)
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

    # brms priors_continuous assigns N(1, 0.5) specifically to odi_baseline_z
    # (hardcoded coef name), not to the outcome-matched baseline. To replicate
    # brms exactly, we apply the informative prior to odi_baseline_z.
    odi_base_col = None
    for i, name in enumerate(col_names_global):
        if name == "odi_baseline_z":
            odi_base_col = i
            break

    n_covs = X.shape[1]

    # Use masked array for PyMC's native missing data handling
    y_masked = np.ma.masked_invalid(y)

    with pm.Model() as model:
        intercept = pm.Normal("intercept", mu=prior_intercept_mu,
                              sigma=prior_intercept_sigma)
        beta_treatment = pm.Normal("beta_treatment", mu=0,
                                   sigma=prior_treatment_sigma)

        # Separate prior for odi_baseline_z
        if odi_base_col is not None:
            beta_odi_baseline = pm.Normal("beta_odi_baseline",
                                          mu=prior_odi_baseline_mu, sigma=0.5)
            # Remaining covariates get N(0, 2)
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

        # HalfStudentT with lam parameterization (PyMC 5.26.1 bug workaround)
        sigma = pm.HalfStudentT("sigma", nu=3,
                                lam=1.0 / (prior_sigma_scale**2))

        pm.Normal("y_obs", mu=mu, sigma=sigma, observed=y_masked)

        idata = pm.sample(
            draws=N_DRAWS, tune=N_TUNE, chains=N_CHAINS,
            target_accept=TARGET_ACCEPT, cores=N_CORES,
            random_seed=SEED, return_inferencedata=True, progressbar=True,
            nuts_sampler_kwargs={"max_treedepth": MAX_TREEDEPTH},
        )

    # --- Diagnostics ---
    diag = check_diagnostics(idata, outcome_label)

    # --- G-computation ---
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

    # Under ELD (treatment=1)
    mu_eld = intercept_s[:, None] + beta_t_s[:, None] * 1.0 + cov_contrib
    # Under MSD (treatment=0)
    mu_msd = intercept_s[:, None] + beta_t_s[:, None] * 0.0 + cov_contrib

    # ATE: higher-is-better (EQ-5D), positive = ELD superior
    ate_samples = mu_eld.mean(axis=1) - mu_msd.mean(axis=1)

    return idata, ate_samples, diag, N, n_obs


# =============================================================================
# 5. BERNOULLI MODEL (binary outcomes)
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
        # Normal(0, 0.5) for regularized shrinkage (stabilizes convergence
        # for binary outcomes with many covariates relative to events)
        beta_covs = pm.Normal("beta_covs", mu=0, sigma=0.5,
                              shape=n_covs)

        logit_p = (intercept
                   + beta_treatment * pt.as_tensor_variable(treatment)
                   + pt.dot(pt.as_tensor_variable(X), beta_covs))

        # Observed-case likelihood
        pm.Bernoulli("y_obs", logit_p=logit_p[obs_idx], observed=y_obs)

        idata = pm.sample(
            draws=N_DRAWS, tune=N_TUNE, chains=N_CHAINS,
            target_accept=TARGET_ACCEPT, cores=N_CORES,
            random_seed=SEED, return_inferencedata=True, progressbar=True,
            nuts_sampler_kwargs={"max_treedepth": MAX_TREEDEPTH},
        )

    # --- Diagnostics ---
    diag = check_diagnostics(idata, outcome_label)

    # --- G-computation ---
    posterior = idata.posterior
    S = N_CHAINS * N_DRAWS

    intercept_s = posterior["intercept"].values.flatten()
    beta_t_s = posterior["beta_treatment"].values.flatten()
    beta_c_s = posterior["beta_covs"].values.reshape(S, -1)

    cov_contrib = beta_c_s @ X.T  # (S, N)

    # Under ELD (treatment=1)
    p_eld = expit(intercept_s[:, None] + beta_t_s[:, None] * 1.0 + cov_contrib)
    # Under MSD (treatment=0)
    p_msd = expit(intercept_s[:, None] + beta_t_s[:, None] * 0.0 + cov_contrib)

    # ATE: sign convention depends on direction
    if lower_is_better:
        ate_samples = p_msd.mean(axis=1) - p_eld.mean(axis=1)
    else:
        ate_samples = p_eld.mean(axis=1) - p_msd.mean(axis=1)

    return idata, ate_samples, diag, N, n_obs


# =============================================================================
# 6. CONVERGENCE DIAGNOSTICS
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
# 7. OUTCOME DEFINITIONS
# =============================================================================

def get_outcome_specs():
    """Define all 17 outcomes with their specifications."""
    outcomes = [
        # ZOIB outcomes (lower-is-better)
        {"var": "odi_3m", "label": "ODI 3 months", "family": "zoib",
         "scale": 100, "zi_baseline_z": "odi_baseline_z",
         "dataset": "full", "margin": NI_MARGINS["odi"],
         "lower_is_better": True},
        {"var": "odi_12m", "label": "ODI 12 months", "family": "zoib",
         "scale": 100, "zi_baseline_z": "odi_baseline_z",
         "dataset": "12m", "margin": NI_MARGINS["odi"],
         "lower_is_better": True},
        {"var": "nrs_back_3m", "label": "NRS back 3 months", "family": "zoib",
         "scale": 10, "zi_baseline_z": "nrs_back_baseline_z",
         "dataset": "full", "margin": NI_MARGINS["nrs"],
         "lower_is_better": True},
        {"var": "nrs_back_12m", "label": "NRS back 12 months", "family": "zoib",
         "scale": 10, "zi_baseline_z": "nrs_back_baseline_z",
         "dataset": "12m", "margin": NI_MARGINS["nrs"],
         "lower_is_better": True},
        {"var": "nrs_leg_3m", "label": "NRS leg 3 months", "family": "zoib",
         "scale": 10, "zi_baseline_z": "nrs_leg_baseline_z",
         "dataset": "full", "margin": NI_MARGINS["nrs"],
         "lower_is_better": True},
        {"var": "nrs_leg_12m", "label": "NRS leg 12 months", "family": "zoib",
         "scale": 10, "zi_baseline_z": "nrs_leg_baseline_z",
         "dataset": "12m", "margin": NI_MARGINS["nrs"],
         "lower_is_better": True},
        # Gaussian outcomes (higher-is-better)
        {"var": "eq5d_3m", "label": "EQ-5D 3 months", "family": "gaussian",
         "dataset": "full", "margin": NI_MARGINS["eq5d"],
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
         "family": "bernoulli", "dataset": "full",
         "margin": NI_MARGINS["binary"], "lower_is_better": False},
        {"var": "rtw_3m", "label": "RTW 3 months",
         "family": "bernoulli", "dataset": "full",
         "margin": NI_MARGINS["binary"], "lower_is_better": False},
        {"var": "rtw_12m", "label": "RTW 12 months",
         "family": "bernoulli", "dataset": "12m",
         "margin": NI_MARGINS["binary"], "lower_is_better": False},
        {"var": "analgesic_3m", "label": "Analgesic 3 months",
         "family": "bernoulli", "dataset": "full",
         "margin": NI_MARGINS["binary"], "lower_is_better": True},
        {"var": "analgesic_12m", "label": "Analgesic 12 months",
         "family": "bernoulli", "dataset": "12m",
         "margin": NI_MARGINS["binary"], "lower_is_better": True},
        {"var": "satisfied_3m", "label": "Satisfaction 3 months",
         "family": "bernoulli", "dataset": "full",
         "margin": NI_MARGINS["binary"], "lower_is_better": False},
        {"var": "satisfied_12m", "label": "Satisfaction 12 months",
         "family": "bernoulli", "dataset": "12m",
         "margin": NI_MARGINS["binary"], "lower_is_better": False},
        {"var": "gpe_success_3m", "label": "GPE 3 months",
         "family": "bernoulli", "dataset": "full",
         "margin": NI_MARGINS["binary"], "lower_is_better": False},
        {"var": "gpe_success_12m", "label": "GPE 12 months",
         "family": "bernoulli", "dataset": "12m",
         "margin": NI_MARGINS["binary"], "lower_is_better": False},
    ]
    return outcomes


# =============================================================================
# 8. MAIN FITTING LOOP
# =============================================================================

def fit_all_outcomes(df_full, df_12m, X_full, X_12m,
                     col_names_full, col_names_12m):
    """Fit all 17 outcomes sequentially."""
    global col_names_global
    outcomes = get_outcome_specs()
    all_results = []
    all_diags = []
    all_ate_samples = {}

    for i, spec in enumerate(outcomes):
        print(f"\n{'='*80}")
        print(f"OUTCOME {i+1}/17: {spec['label']} ({spec['family']})")
        print(f"{'='*80}")
        t0 = time.time()

        # Select dataset and update col_names_global for model functions
        if spec["dataset"] == "12m":
            df = df_12m
            X = X_12m
            col_names_global = col_names_12m
        else:
            df = df_full
            X = X_full
            col_names_global = col_names_full

        # Fit model
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

            # Store results
            result = {
                "Outcome": spec["label"],
                "Family": spec["family"].upper(),
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

            # Save InferenceData and ATE samples
            nc_path = os.path.join(MODEL_DIR,
                                   f"pymc_idata_{spec['var']}.nc")
            idata.to_netcdf(nc_path)
            np.save(os.path.join(MODEL_DIR,
                                 f"pymc_ate_{spec['var']}.npy"), ate)

            # Free memory
            del idata
            import gc
            gc.collect()

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            all_results.append({
                "Outcome": spec["label"],
                "Family": spec["family"].upper(),
                "N_total": len(df),
                "N_observed": int((~np.isnan(
                    df[spec["var"]].values.astype(float))).sum()),
                "N_gcomp": len(df),
                "ATE": np.nan,
                "CrI_lo": np.nan,
                "CrI_hi": np.nan,
                "P_NI": np.nan,
                "P_Superiority": np.nan,
                "NI_Margin": spec["margin"],
                "NI_Conclusion": "FAILED",
                "Convergence": "FAILED",
            })
            all_diags.append({
                "outcome": spec["label"],
                "rhat_max": np.nan,
                "ess_bulk_min": np.nan,
                "ess_tail_min": np.nan,
                "n_divergent": np.nan,
                "valid": False,
            })

    return all_results, all_diags, all_ate_samples


# =============================================================================
# 9. COMPARISON WITH BRMS RESULTS
# =============================================================================

def build_comparison_table(pymc_results):
    """Build side-by-side comparison with brms results."""
    print("\n" + "=" * 80)
    print("COMPARISON WITH BRMS RESULTS")
    print("=" * 80)

    # Load brms tables
    brms_primary_path = os.path.join(TABLE_DIR, "table2_primary_results.csv")
    brms_tier2_path = os.path.join(TABLE_DIR, "table3_tier2_effectiveness.csv")

    brms_rows = []

    # Primary outcome
    if os.path.exists(brms_primary_path):
        brms_primary = pd.read_csv(brms_primary_path)
        for _, row in brms_primary.iterrows():
            cri = row["95% CrI"].strip("[]").split(",")
            brms_rows.append({
                "Outcome": row["Outcome"],
                "brms_ATE": float(row["ATE"]),
                "brms_CrI_lo": float(cri[0].strip()),
                "brms_CrI_hi": float(cri[1].strip()),
                "brms_P_NI": float(row["P(NI)"]),
                "brms_N_gcomp": int(row["N_ELD"]) + int(row["N_MSD"]),
            })

    # Tier 2 outcomes
    if os.path.exists(brms_tier2_path):
        brms_tier2 = pd.read_csv(brms_tier2_path)
        for _, row in brms_tier2.iterrows():
            brms_rows.append({
                "Outcome": row["Outcome"],
                "brms_ATE": float(row["ATE"]),
                "brms_CrI_lo": float(row["CrI_lo"]),
                "brms_CrI_hi": float(row["CrI_hi"]),
                "brms_P_NI": float(row["P_NI"]),
                "brms_N_gcomp": int(row["N"]),
            })

    brms_df = pd.DataFrame(brms_rows)
    pymc_df = pd.DataFrame(pymc_results)

    # Match outcomes by label
    # Create mapping from brms outcome names to PyMC outcome names
    name_map = {
        "ODI 3 months": "ODI 3 months",
        "ODI 12 months": "ODI 12 months",
        "NRS back pain 3 months": "NRS back 3 months",
        "NRS back pain 12 months": "NRS back 12 months",
        "NRS leg pain 3 months": "NRS leg 3 months",
        "NRS leg pain 12 months": "NRS leg 12 months",
        "EQ-5D 3 months": "EQ-5D 3 months",
        "EQ-5D 12 months": "EQ-5D 12 months",
        "Responder 3 months (>=30% or >=10pt)": "Responder 3 months",
        "Return to work 3 months": "RTW 3 months",
        "Return to work 12 months": "RTW 12 months",
        "Analgesic use 3 months": "Analgesic 3 months",
        "Analgesic use 12 months": "Analgesic 12 months",
        "Satisfaction 3 months": "Satisfaction 3 months",
        "Satisfaction 12 months": "Satisfaction 12 months",
        "GPE success 3 months": "GPE 3 months",
        "GPE success 12 months": "GPE 12 months",
    }

    comparison_rows = []
    for _, brms_row in brms_df.iterrows():
        pymc_name = name_map.get(brms_row["Outcome"], brms_row["Outcome"])
        pymc_match = pymc_df[pymc_df["Outcome"] == pymc_name]

        if len(pymc_match) > 0:
            pm_row = pymc_match.iloc[0]
            # Concordance: both agree on NI conclusion
            brms_ni = brms_row["brms_P_NI"] > NI_THRESHOLD
            pymc_ni = pm_row["P_NI"] > NI_THRESHOLD
            concordance = "Concordant" if brms_ni == pymc_ni else "Discordant"

            comparison_rows.append({
                "Outcome": brms_row["Outcome"],
                "brms_ATE": round(brms_row["brms_ATE"], 4),
                "brms_CrI_lo": round(brms_row["brms_CrI_lo"], 4),
                "brms_CrI_hi": round(brms_row["brms_CrI_hi"], 4),
                "brms_P_NI": round(brms_row["brms_P_NI"], 4),
                "brms_N_gcomp": brms_row["brms_N_gcomp"],
                "PyMC_ATE": pm_row["ATE"],
                "PyMC_CrI_lo": pm_row["CrI_lo"],
                "PyMC_CrI_hi": pm_row["CrI_hi"],
                "PyMC_P_NI": pm_row["P_NI"],
                "PyMC_N_gcomp": pm_row["N_gcomp"],
                "Concordance": concordance,
            })

    comp_df = pd.DataFrame(comparison_rows)

    n_concordant = (comp_df["Concordance"] == "Concordant").sum()
    print(f"\nConcordance: {n_concordant}/{len(comp_df)} outcomes agree on NI")

    for _, row in comp_df.iterrows():
        print(f"  {row['Outcome']}: brms ATE={row['brms_ATE']:.3f} "
              f"(N={row['brms_N_gcomp']}), PyMC ATE={row['PyMC_ATE']:.3f} "
              f"(N={row['PyMC_N_gcomp']}), {row['Concordance']}")

    return comp_df


# =============================================================================
# 10. FIGURES
# =============================================================================

def plot_forest_comparison(comparison_df, pymc_results):
    """Publication forest plot comparing brms and PyMC estimates."""
    print("\n--- Creating forest comparison plot ---")

    comp = comparison_df.copy()

    # Separate continuous and binary outcomes
    continuous_names = ["ODI 3 months", "ODI 12 months",
                        "NRS back pain 3 months", "NRS back pain 12 months",
                        "NRS leg pain 3 months", "NRS leg pain 12 months",
                        "EQ-5D 3 months", "EQ-5D 12 months"]
    continuous = comp[comp["Outcome"].isin(continuous_names)].copy()
    binary = comp[~comp["Outcome"].isin(continuous_names)].copy()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 10),
                                    gridspec_kw={"width_ratios": [1, 1]})

    # --- Panel A: Continuous outcomes ---
    if len(continuous) > 0:
        continuous = continuous.iloc[::-1].reset_index(drop=True)
        n = len(continuous)
        y_positions = np.arange(n) * 2.5

        for i, (_, row) in enumerate(continuous.iterrows()):
            y = y_positions[i]
            # brms (blue)
            ax1.errorbar(row["brms_ATE"], y + 0.4,
                         xerr=[[row["brms_ATE"] - row["brms_CrI_lo"]],
                               [row["brms_CrI_hi"] - row["brms_ATE"]]],
                         fmt="o", color="#2166AC", markersize=6,
                         capsize=3, linewidth=1.2, label="brms" if i == 0 else "")
            # PyMC (orange)
            ax1.errorbar(row["PyMC_ATE"], y - 0.4,
                         xerr=[[row["PyMC_ATE"] - row["PyMC_CrI_lo"]],
                               [row["PyMC_CrI_hi"] - row["PyMC_ATE"]]],
                         fmt="s", color="#E66100", markersize=6,
                         capsize=3, linewidth=1.2, label="PyMC" if i == 0 else "")

        ax1.axvline(x=0, color="0.30", linestyle="--", linewidth=0.8)
        ax1.set_yticks(y_positions)
        ax1.set_yticklabels(continuous["Outcome"].values)
        ax1.set_xlabel("ATE (positive = ELD superior)")
        ax1.set_title("A. Continuous Outcomes", fontweight="bold", fontsize=12)
        ax1.legend(loc="lower right", framealpha=0.9)
        ax1.spines["right"].set_visible(False)
        ax1.spines["top"].set_visible(False)

    # --- Panel B: Binary outcomes ---
    if len(binary) > 0:
        binary = binary.iloc[::-1].reset_index(drop=True)
        n = len(binary)
        y_positions = np.arange(n) * 2.5

        for i, (_, row) in enumerate(binary.iterrows()):
            y = y_positions[i]
            # brms (blue)
            ax2.errorbar(row["brms_ATE"], y + 0.4,
                         xerr=[[row["brms_ATE"] - row["brms_CrI_lo"]],
                               [row["brms_CrI_hi"] - row["brms_ATE"]]],
                         fmt="o", color="#2166AC", markersize=6,
                         capsize=3, linewidth=1.2, label="brms" if i == 0 else "")
            # PyMC (orange)
            ax2.errorbar(row["PyMC_ATE"], y - 0.4,
                         xerr=[[row["PyMC_ATE"] - row["PyMC_CrI_lo"]],
                               [row["PyMC_CrI_hi"] - row["PyMC_ATE"]]],
                         fmt="s", color="#E66100", markersize=6,
                         capsize=3, linewidth=1.2, label="PyMC" if i == 0 else "")

        ax2.axvline(x=0, color="0.30", linestyle="--", linewidth=0.8)
        ax2.axvline(x=-NI_MARGINS["binary"], color="red", linestyle=":",
                    linewidth=0.8, alpha=0.7)
        ax2.set_yticks(y_positions)
        ax2.set_yticklabels(binary["Outcome"].values)
        ax2.set_xlabel("ATE (positive = ELD superior)")
        ax2.set_title("B. Binary Outcomes", fontweight="bold", fontsize=12)
        ax2.legend(loc="lower right", framealpha=0.9)
        ax2.spines["right"].set_visible(False)
        ax2.spines["top"].set_visible(False)

    plt.tight_layout()
    path = os.path.join(FIG_DIR, "fig_brms_vs_pymc_forest.png")
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {path}")

    # Copy to legacy path
    import shutil
    legacy_path = os.path.join(LEGACY_FIG_DIR, "fig_brms_vs_pymc_forest.png")
    try:
        shutil.copy2(path, legacy_path)
    except Exception:
        pass


def plot_posterior_primary(ate_samples_dict):
    """Plot posterior density for the primary outcome (ODI 3m)."""
    print("\n--- Creating primary posterior plot ---")

    if "odi_3m" not in ate_samples_dict:
        print("  ODI 3m ATE samples not available, skipping")
        return

    ate = ate_samples_dict["odi_3m"]
    ate_mean = ate.mean()
    ate_ci = [np.percentile(ate, 2.5), np.percentile(ate, 97.5)]
    p_ni = (ate > -NI_MARGINS["odi"]).mean()

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(ate, bins=80, density=True, alpha=0.4, color="#92C5DE",
            edgecolor="none")

    # KDE overlay
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(ate)
    x_grid = np.linspace(ate.min() - 2, ate.max() + 2, 500)
    ax.plot(x_grid, kde(x_grid), color="#2166AC", linewidth=1.5)

    # Reference lines
    ax.axvline(x=0, color="0.30", linestyle="-", linewidth=0.8)
    ax.axvline(x=-NI_MARGINS["odi"], color="red", linestyle="--",
               linewidth=1.2)

    # Annotations
    ax.annotate(f"NI margin = {-NI_MARGINS['odi']:.0f}",
                xy=(-NI_MARGINS["odi"], 0), xytext=(-NI_MARGINS["odi"] - 1, 0.02),
                fontsize=9, color="red", ha="right")
    ax.text(0.97, 0.95,
            f"ATE = {ate_mean:.2f}\n"
            f"95% CrI: [{ate_ci[0]:.2f}, {ate_ci[1]:.2f}]\n"
            f"P(NI) = {p_ni:.4f}",
            transform=ax.transAxes, fontsize=10, va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="0.70", alpha=0.9))

    ax.set_xlabel("ATE: $\\mu_{MSD} - \\mu_{ELD}$ (ODI points)")
    ax.set_ylabel("Posterior Density")
    ax.set_title("PyMC Posterior: Treatment Effect on ODI at 3 Months",
                 fontweight="bold")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    path = os.path.join(FIG_DIR, "fig_pymc_posterior_primary.png")
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_diagnostics(ate_samples_dict):
    """Plot trace-style diagnostics for key outcomes."""
    print("\n--- Creating diagnostics plot ---")

    key_outcomes = ["odi_3m", "nrs_back_3m", "eq5d_3m", "responder_3m"]
    key_labels = ["ODI 3m", "NRS back 3m", "EQ-5D 3m", "Responder 3m"]

    # Load InferenceData for key outcomes
    available = []
    for var, label in zip(key_outcomes, key_labels):
        nc_path = os.path.join(MODEL_DIR, f"pymc_idata_{var}.nc")
        if os.path.exists(nc_path):
            available.append((var, label, nc_path))

    if not available:
        print("  No InferenceData files available, skipping")
        return

    n_plots = len(available)
    fig, axes = plt.subplots(n_plots, 2, figsize=(12, 3 * n_plots))
    if n_plots == 1:
        axes = axes.reshape(1, -1)

    for i, (var, label, nc_path) in enumerate(available):
        idata = az.from_netcdf(nc_path)
        post = idata.posterior["beta_treatment"]

        # Trace plot
        for chain in range(post.shape[0]):
            axes[i, 0].plot(post.values[chain, :], alpha=0.5, linewidth=0.3)
        axes[i, 0].set_ylabel(f"beta_treatment")
        axes[i, 0].set_title(f"{label}: Trace", fontsize=10)
        axes[i, 0].spines["right"].set_visible(False)
        axes[i, 0].spines["top"].set_visible(False)

        # Posterior density
        vals = post.values.flatten()
        axes[i, 1].hist(vals, bins=60, density=True, alpha=0.6,
                        color="#92C5DE", edgecolor="none")
        axes[i, 1].axvline(x=0, color="0.30", linestyle="--", linewidth=0.8)
        axes[i, 1].set_title(f"{label}: Posterior", fontsize=10)
        axes[i, 1].spines["right"].set_visible(False)
        axes[i, 1].spines["top"].set_visible(False)

        del idata

    axes[-1, 0].set_xlabel("Iteration")
    axes[-1, 1].set_xlabel("beta_treatment")
    plt.tight_layout()

    path = os.path.join(FIG_DIR, "fig_pymc_diagnostics.png")
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {path}")


# =============================================================================
# 11. MAIN
# =============================================================================

def main():
    global col_names_global  # needed by ZOIB model to find ZI baseline var

    print("=" * 80)
    print("ENDO-LUMBAR: SUPPLEMENTARY PyMC BAYESIAN ANALYSIS")
    print("=" * 80)
    total_start = time.time()

    # Load data
    df_full, df_12m = load_and_prepare_data()

    # Build design matrices
    print("\n--- Building design matrix (full sample) ---")
    X_full, col_names_full, std_params_full = build_design_matrix(df_full)

    print("\n--- Building design matrix (12m sample) ---")
    X_12m, col_names_12m, std_params_12m = build_design_matrix(df_12m)

    # Store globally for model functions to access
    col_names_global = col_names_full

    # Report missingness
    print("\nOutcome missingness:")
    outcomes = get_outcome_specs()
    for spec in outcomes:
        df = df_full if spec["dataset"] == "full" else df_12m
        y = df[spec["var"]].values.astype(float)
        n_miss = np.isnan(y).sum()
        N = len(df)
        print(f"  {spec['label']:30s}: {n_miss:3d}/{N} missing ({100*n_miss/N:.1f}%)")

    # Fit all 17 outcomes
    all_results, all_diags, ate_samples_dict = fit_all_outcomes(
        df_full, df_12m, X_full, X_12m, col_names_full, col_names_12m)

    # --- Save results tables ---
    print("\n" + "=" * 80)
    print("SAVING RESULTS")
    print("=" * 80)

    # Table 1: PyMC results
    results_df = pd.DataFrame(all_results)
    results_path = os.path.join(TABLE_DIR, "table_pymc_results.csv")
    results_df.to_csv(results_path, index=False)
    print(f"Saved: {results_path}")

    # Table 2: Convergence diagnostics
    diag_df = pd.DataFrame(all_diags)
    diag_path = os.path.join(TABLE_DIR, "table_pymc_convergence.csv")
    diag_df.to_csv(diag_path, index=False)
    print(f"Saved: {diag_path}")

    # Table 3: brms vs PyMC comparison
    comparison_df = build_comparison_table(all_results)
    comp_path = os.path.join(TABLE_DIR, "table_brms_vs_pymc.csv")
    comparison_df.to_csv(comp_path, index=False)
    print(f"Saved: {comp_path}")

    # --- Figures ---
    print("\n" + "=" * 80)
    print("CREATING FIGURES")
    print("=" * 80)

    if len(comparison_df) > 0:
        plot_forest_comparison(comparison_df, all_results)
    plot_posterior_primary(ate_samples_dict)
    plot_diagnostics(ate_samples_dict)

    # --- Save results pickle ---
    results_pkl = {
        "pymc_results": all_results,
        "diagnostics": all_diags,
        "ate_samples": {k: v.tolist() for k, v in ate_samples_dict.items()},
        "comparison": comparison_df.to_dict() if len(comparison_df) > 0 else {},
        "std_params_full": std_params_full,
        "std_params_12m": std_params_12m,
        "col_names": col_names_full,
    }
    pkl_path = os.path.join(MODEL_DIR, "pymc_results.pkl")
    with open(pkl_path, "wb") as f:
        pickle.dump(results_pkl, f)
    print(f"Saved: {pkl_path}")

    # --- Copy key outputs to legacy paths ---
    import shutil
    for src_name in ["table_pymc_results.csv", "table_brms_vs_pymc.csv",
                     "table_pymc_convergence.csv"]:
        src = os.path.join(TABLE_DIR, src_name)
        dst = os.path.join(LEGACY_TABLE_DIR, src_name)
        try:
            shutil.copy2(src, dst)
        except Exception:
            pass

    # --- Final summary ---
    total_elapsed = time.time() - total_start
    print(f"\n{'='*80}")
    print(f"ANALYSIS COMPLETE ({total_elapsed/60:.1f} minutes)")
    print(f"{'='*80}")

    print(f"\n{'Outcome':<30s} {'ATE':>8s} {'95% CrI':>20s} "
          f"{'P(NI)':>8s} {'N_gcomp':>8s}")
    print("-" * 80)
    for r in all_results:
        cri_str = f"[{r['CrI_lo']:.3f}, {r['CrI_hi']:.3f}]"
        print(f"{r['Outcome']:<30s} {r['ATE']:>8.3f} {cri_str:>20s} "
              f"{r['P_NI']:>8.4f} {r['N_gcomp']:>8d}")

    n_ni = sum(1 for r in all_results
               if r["NI_Conclusion"] == "NI demonstrated")
    print(f"\nNI demonstrated: {n_ni}/{len(all_results)} outcomes")
    print(f"All outputs saved to: {BASE_DIR}")

    return all_results, comparison_df


if __name__ == "__main__":
    all_results, comparison_df = main()
