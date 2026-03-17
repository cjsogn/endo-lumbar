#!/usr/bin/env python3
"""
22b_refit_bernoulli.py - Refit Bernoulli models with convergence issues.

The Cauchy(0,1) prior on 29 covariates causes divergences in several binary
outcomes. Fix: increase tune to 2000, raise target_accept to 0.99.
Updates the saved results tables and figures.
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
import os
import time
import pickle
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# CONFIGURATION (same as 22_pymc_bayesian.py)
# =============================================================================

SEED = 20260204
N_CORES = 14
N_CHAINS = 4
N_DRAWS = 1000
N_TUNE = 2000       # Doubled from 1000
TARGET_ACCEPT = 0.99  # Raised from 0.95
MAX_TREEDEPTH = 12
NI_THRESHOLD = 0.95

BASE_DIR = "/Users/cjsogn/ENDO_LUMBAR"
DATA_DIR = os.path.join(BASE_DIR, "data")
TABLE_DIR = os.path.join(BASE_DIR, "tables")
FIG_DIR = os.path.join(BASE_DIR, "figures")
MODEL_DIR = os.path.join(BASE_DIR, "models")
LEGACY_TABLE_DIR = "/Users/cjsogn/endo_studies/lumbar/analysis/output/tables"
LEGACY_FIG_DIR = "/Users/cjsogn/endo_studies/lumbar/analysis/output/figures"

NI_MARGINS = {"odi": 7.0, "nrs": 1.0, "eq5d": 0.05, "binary": 0.10}

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

    # Ordinal covariates (integer-coded, treated as linear)
    for var in ORDINAL_COVS:
        vals = df[var].values.astype(float)
        vals = np.where(np.isnan(vals), 0.0, vals)
        cols.append(vals)
        col_names.append(var)

    X = np.column_stack(cols)
    return X, col_names, std_params


def fit_bernoulli_model(df, X, outcome_var, outcome_label, lower_is_better):
    """Fit Bernoulli model with improved sampling settings."""
    print(f"\n--- Bernoulli Model (refit): {outcome_label} ---")

    N = len(df)
    treatment = df["treatment_num"].values.astype(float)
    y = df[outcome_var].values.astype(float)

    obs_mask = ~np.isnan(y)
    obs_idx = np.where(obs_mask)[0]
    y_obs = y[obs_mask]
    n_obs = len(y_obs)
    print(f"  N={N}, Observed={n_obs} ({100*n_obs/N:.1f}%)")
    print(f"  Settings: tune={N_TUNE}, target_accept={TARGET_ACCEPT}")

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

    # Diagnostics
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
    valid = max_rhat < 1.01 and min_bulk > 400 and min_tail > 400 and n_divergent == 0

    status = "PASS" if valid else "ISSUE"
    print(f"  Convergence: {status} (Rhat_max={max_rhat:.4f}, "
          f"ESS_bulk_min={min_bulk:.0f}, ESS_tail_min={min_tail:.0f}, "
          f"Divergences={n_divergent})")

    diag = {
        "outcome": outcome_label,
        "rhat_max": float(max_rhat),
        "ess_bulk_min": float(min_bulk),
        "ess_tail_min": float(min_tail),
        "n_divergent": n_divergent,
        "valid": valid,
    }

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


def main():
    print("=" * 80)
    print("REFITTING BERNOULLI MODELS WITH IMPROVED SAMPLING")
    print(f"  tune={N_TUNE} (was 1000), target_accept={TARGET_ACCEPT} (was 0.95)")
    print("=" * 80)
    total_start = time.time()

    # Load data
    df_full = pd.read_csv(os.path.join(DATA_DIR, "df_disc_full.csv"))
    df_12m = pd.read_csv(os.path.join(DATA_DIR, "df_disc_12m.csv"))

    # Build design matrices
    X_full, col_names_full, _ = build_design_matrix(df_full)
    X_12m, col_names_12m, _ = build_design_matrix(df_12m)

    # Models to refit (those with convergence issues)
    to_refit = [
        {"var": "rtw_3m", "label": "RTW 3 months",
         "dataset": "full", "margin": NI_MARGINS["binary"],
         "lower_is_better": False},
        {"var": "rtw_12m", "label": "RTW 12 months",
         "dataset": "12m", "margin": NI_MARGINS["binary"],
         "lower_is_better": False},
        {"var": "analgesic_3m", "label": "Analgesic 3 months",
         "dataset": "full", "margin": NI_MARGINS["binary"],
         "lower_is_better": True},
        {"var": "analgesic_12m", "label": "Analgesic 12 months",
         "dataset": "12m", "margin": NI_MARGINS["binary"],
         "lower_is_better": True},
        {"var": "satisfied_3m", "label": "Satisfaction 3 months",
         "dataset": "full", "margin": NI_MARGINS["binary"],
         "lower_is_better": False},
        {"var": "satisfied_12m", "label": "Satisfaction 12 months",
         "dataset": "12m", "margin": NI_MARGINS["binary"],
         "lower_is_better": False},
        {"var": "gpe_success_12m", "label": "GPE 12 months",
         "dataset": "12m", "margin": NI_MARGINS["binary"],
         "lower_is_better": False},
    ]

    refit_results = {}

    for spec in to_refit:
        if spec["dataset"] == "12m":
            df, X = df_12m, X_12m
        else:
            df, X = df_full, X_full

        t0 = time.time()
        idata, ate, diag, N, n_obs = fit_bernoulli_model(
            df, X, spec["var"], spec["label"], spec["lower_is_better"])

        ate_mean = float(ate.mean())
        ate_ci = [float(np.percentile(ate, 2.5)),
                  float(np.percentile(ate, 97.5))]
        p_ni = float((ate > -spec["margin"]).mean())
        p_sup = float((ate > 0).mean())

        ni_conclusion = "NI demonstrated" if p_ni > NI_THRESHOLD else (
            "Inconclusive" if p_ni >= 0.80 else "Concern")

        print(f"  ATE = {ate_mean:.4f} (95% CrI: [{ate_ci[0]:.4f}, {ate_ci[1]:.4f}])")
        print(f"  P(NI) = {p_ni:.4f}, P(Superiority) = {p_sup:.4f}")
        print(f"  Time: {(time.time()-t0)/60:.1f} min")

        refit_results[spec["var"]] = {
            "result": {
                "Outcome": spec["label"],
                "Family": "BERNOULLI",
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
            },
            "diag": diag,
            "ate": ate,
        }

        # Save updated InferenceData and ATE
        idata.to_netcdf(os.path.join(MODEL_DIR, f"pymc_idata_{spec['var']}.nc"))
        np.save(os.path.join(MODEL_DIR, f"pymc_ate_{spec['var']}.npy"), ate)

        del idata
        import gc
        gc.collect()

    # Update results table
    print("\n--- Updating results tables ---")
    results_df = pd.read_csv(os.path.join(TABLE_DIR, "table_pymc_results.csv"))
    diag_df = pd.read_csv(os.path.join(TABLE_DIR, "table_pymc_convergence.csv"))

    for var, data in refit_results.items():
        r = data["result"]
        d = data["diag"]

        # Update results table
        mask = results_df["Outcome"] == r["Outcome"]
        if mask.any():
            for col, val in r.items():
                results_df.loc[mask, col] = val

        # Update convergence table
        mask_d = diag_df["outcome"] == d["outcome"]
        if mask_d.any():
            for col, val in d.items():
                diag_df.loc[mask_d, col] = val

    results_df.to_csv(os.path.join(TABLE_DIR, "table_pymc_results.csv"),
                      index=False)
    diag_df.to_csv(os.path.join(TABLE_DIR, "table_pymc_convergence.csv"),
                   index=False)

    # Rebuild comparison table
    print("\n--- Rebuilding comparison table ---")
    brms_primary = pd.read_csv(os.path.join(TABLE_DIR, "table2_primary_results.csv"))
    brms_tier2 = pd.read_csv(os.path.join(TABLE_DIR, "table3_tier2_effectiveness.csv"))

    brms_rows = []
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
        pymc_match = results_df[results_df["Outcome"] == pymc_name]
        if len(pymc_match) > 0:
            pm_row = pymc_match.iloc[0]
            brms_ni = brms_row["brms_P_NI"] > NI_THRESHOLD
            pymc_ni = pm_row["P_NI"] > NI_THRESHOLD
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
                "Concordance": "Concordant" if brms_ni == pymc_ni else "Discordant",
            })

    comp_df = pd.DataFrame(comparison_rows)
    comp_df.to_csv(os.path.join(TABLE_DIR, "table_brms_vs_pymc.csv"), index=False)

    # Rebuild forest plot
    print("\n--- Rebuilding forest plot ---")
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

    continuous_names = ["ODI 3 months", "ODI 12 months",
                        "NRS back pain 3 months", "NRS back pain 12 months",
                        "NRS leg pain 3 months", "NRS leg pain 12 months",
                        "EQ-5D 3 months", "EQ-5D 12 months"]
    continuous = comp_df[comp_df["Outcome"].isin(continuous_names)].iloc[::-1].reset_index(drop=True)
    binary = comp_df[~comp_df["Outcome"].isin(continuous_names)].iloc[::-1].reset_index(drop=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 10),
                                    gridspec_kw={"width_ratios": [1, 1]})

    if len(continuous) > 0:
        y_pos = np.arange(len(continuous)) * 2.5
        for i, (_, row) in enumerate(continuous.iterrows()):
            y = y_pos[i]
            ax1.errorbar(row["brms_ATE"], y + 0.4,
                         xerr=[[row["brms_ATE"] - row["brms_CrI_lo"]],
                               [row["brms_CrI_hi"] - row["brms_ATE"]]],
                         fmt="o", color="#2166AC", markersize=6,
                         capsize=3, linewidth=1.2, label="brms" if i == 0 else "")
            ax1.errorbar(row["PyMC_ATE"], y - 0.4,
                         xerr=[[row["PyMC_ATE"] - row["PyMC_CrI_lo"]],
                               [row["PyMC_CrI_hi"] - row["PyMC_ATE"]]],
                         fmt="s", color="#E66100", markersize=6,
                         capsize=3, linewidth=1.2, label="PyMC" if i == 0 else "")
        ax1.axvline(x=0, color="0.30", linestyle="--", linewidth=0.8)
        ax1.set_yticks(y_pos)
        ax1.set_yticklabels(continuous["Outcome"].values)
        ax1.set_xlabel("ATE (positive = ELD superior)")
        ax1.set_title("A. Continuous Outcomes", fontweight="bold", fontsize=12)
        ax1.legend(loc="lower right", framealpha=0.9)
        ax1.spines["right"].set_visible(False)
        ax1.spines["top"].set_visible(False)

    if len(binary) > 0:
        y_pos = np.arange(len(binary)) * 2.5
        for i, (_, row) in enumerate(binary.iterrows()):
            y = y_pos[i]
            ax2.errorbar(row["brms_ATE"], y + 0.4,
                         xerr=[[row["brms_ATE"] - row["brms_CrI_lo"]],
                               [row["brms_CrI_hi"] - row["brms_ATE"]]],
                         fmt="o", color="#2166AC", markersize=6,
                         capsize=3, linewidth=1.2, label="brms" if i == 0 else "")
            ax2.errorbar(row["PyMC_ATE"], y - 0.4,
                         xerr=[[row["PyMC_ATE"] - row["PyMC_CrI_lo"]],
                               [row["PyMC_CrI_hi"] - row["PyMC_ATE"]]],
                         fmt="s", color="#E66100", markersize=6,
                         capsize=3, linewidth=1.2, label="PyMC" if i == 0 else "")
        ax2.axvline(x=0, color="0.30", linestyle="--", linewidth=0.8)
        ax2.axvline(x=-NI_MARGINS["binary"], color="red", linestyle=":",
                    linewidth=0.8, alpha=0.7)
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(binary["Outcome"].values)
        ax2.set_xlabel("ATE (positive = ELD superior)")
        ax2.set_title("B. Binary Outcomes", fontweight="bold", fontsize=12)
        ax2.legend(loc="lower right", framealpha=0.9)
        ax2.spines["right"].set_visible(False)
        ax2.spines["top"].set_visible(False)

    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "fig_brms_vs_pymc_forest.png"),
                dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # Copy to legacy
    import shutil
    for src_name in ["table_pymc_results.csv", "table_brms_vs_pymc.csv",
                     "table_pymc_convergence.csv"]:
        try:
            shutil.copy2(os.path.join(TABLE_DIR, src_name),
                         os.path.join(LEGACY_TABLE_DIR, src_name))
        except Exception:
            pass
    try:
        shutil.copy2(os.path.join(FIG_DIR, "fig_brms_vs_pymc_forest.png"),
                     os.path.join(LEGACY_FIG_DIR, "fig_brms_vs_pymc_forest.png"))
    except Exception:
        pass

    elapsed = time.time() - total_start
    print(f"\n{'='*80}")
    print(f"REFIT COMPLETE ({elapsed/60:.1f} minutes)")
    print(f"{'='*80}")

    # Summary
    print(f"\n{'Outcome':<30s} {'Div_old':>8s} {'Div_new':>8s} {'Rhat':>8s} {'ESS_tail':>10s} {'Status':>8s}")
    print("-" * 80)

    old_diag = {
        "RTW 3 months": {"div": 0, "rhat": 1.012},
        "RTW 12 months": {"div": 8, "rhat": 1.005},
        "Analgesic 3 months": {"div": 0, "rhat": 1.013},
        "Analgesic 12 months": {"div": 91, "rhat": 1.011},
        "Satisfaction 3 months": {"div": 69, "rhat": 1.014},
        "Satisfaction 12 months": {"div": 122, "rhat": 1.007},
        "GPE 12 months": {"div": 57, "rhat": 1.007},
    }

    for var, data in refit_results.items():
        d = data["diag"]
        old = old_diag.get(d["outcome"], {"div": "?", "rhat": "?"})
        status = "FIXED" if d["valid"] else "STILL"
        print(f"{d['outcome']:<30s} {old['div']:>8} {d['n_divergent']:>8d} "
              f"{d['rhat_max']:>8.4f} {d['ess_tail_min']:>10.0f} {status:>8s}")


if __name__ == "__main__":
    main()
