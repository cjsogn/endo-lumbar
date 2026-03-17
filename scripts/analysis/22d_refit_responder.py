#!/usr/bin/env python3
"""
22d_refit_responder.py - Refit Responder 3m (ESS_tail=223 < 400 threshold).

The original run had 0 divergences but ESS_tail=223, which is below the 400
minimum required. The validation bug (missing ESS_tail check) masked this.
Fix: increase tune to 2000, raise target_accept to 0.99.
Also regenerates tables with corrected validation logic (includes ESS_tail).
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
import shutil
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# CONFIGURATION
# =============================================================================

SEED = 20260204
N_CORES = 14
N_CHAINS = 4
N_DRAWS = 1000
N_TUNE = 2000           # Increased from 1000
TARGET_ACCEPT = 0.99    # Increased from 0.95
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


def main():
    print("=" * 80)
    print("REFITTING RESPONDER 3M (ESS_tail=223 < 400 threshold)")
    print(f"  tune={N_TUNE}, target_accept={TARGET_ACCEPT}")
    print("=" * 80)
    total_start = time.time()

    # Load data
    df_full = pd.read_csv(os.path.join(DATA_DIR, "df_disc_full.csv"))
    X_full, col_names_full, _ = build_design_matrix(df_full)

    # Fit Responder 3m
    outcome_var = "responder_3m"
    outcome_label = "Responder 3 months"
    lower_is_better = False
    margin = NI_MARGINS["binary"]

    N = len(df_full)
    treatment = df_full["treatment_num"].values.astype(float)
    y = df_full[outcome_var].values.astype(float)

    obs_mask = ~np.isnan(y)
    obs_idx = np.where(obs_mask)[0]
    y_obs = y[obs_mask]
    n_obs = len(y_obs)
    n_covs = X_full.shape[1]

    print(f"\n--- Bernoulli Model (refit): {outcome_label} ---")
    print(f"  N={N}, Observed={n_obs} ({100*n_obs/N:.1f}%)")

    t0 = time.time()

    with pm.Model() as model:
        intercept = pm.Normal("intercept", mu=0, sigma=5)
        beta_treatment = pm.Normal("beta_treatment", mu=0, sigma=1)
        beta_covs = pm.Normal("beta_covs", mu=0, sigma=0.5,
                              shape=n_covs)

        logit_p = (intercept
                   + beta_treatment * pt.as_tensor_variable(treatment)
                   + pt.dot(pt.as_tensor_variable(X_full), beta_covs))

        pm.Bernoulli("y_obs", logit_p=logit_p[obs_idx], observed=y_obs)

        idata = pm.sample(
            draws=N_DRAWS, tune=N_TUNE, chains=N_CHAINS,
            target_accept=TARGET_ACCEPT, cores=N_CORES,
            random_seed=SEED, return_inferencedata=True, progressbar=True,
            nuts_sampler_kwargs={"max_treedepth": MAX_TREEDEPTH},
        )

    # Diagnostics (with corrected ESS_tail check)
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

    # Corrected validation: includes ESS_tail >= 400
    valid = (max_rhat < 1.01 and min_bulk > 400 and min_tail > 400
             and n_divergent == 0)

    status = "PASS" if valid else "ISSUE"
    print(f"  Convergence: {status}")
    print(f"    Rhat_max={max_rhat:.4f}")
    print(f"    ESS_bulk_min={min_bulk:.0f}")
    print(f"    ESS_tail_min={min_tail:.0f} (was 223, threshold=400)")
    print(f"    Divergences={n_divergent}")
    print(f"  Fit time: {(time.time()-t0)/60:.1f} min")

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

    cov_contrib = beta_c_s @ X_full.T
    p_eld = expit(intercept_s[:, None] + beta_t_s[:, None] * 1.0 + cov_contrib)
    p_msd = expit(intercept_s[:, None] + beta_t_s[:, None] * 0.0 + cov_contrib)

    ate_samples = p_eld.mean(axis=1) - p_msd.mean(axis=1)

    ate_mean = float(ate_samples.mean())
    ate_ci = [float(np.percentile(ate_samples, 2.5)),
              float(np.percentile(ate_samples, 97.5))]
    p_ni = float((ate_samples > -margin).mean())
    p_sup = float((ate_samples > 0).mean())
    ni_conclusion = ("NI demonstrated" if p_ni > NI_THRESHOLD
                     else ("Inconclusive" if p_ni >= 0.80 else "Concern"))

    print(f"  ATE = {ate_mean:.4f} (95% CrI: [{ate_ci[0]:.4f}, {ate_ci[1]:.4f}])")
    print(f"  P(NI) = {p_ni:.4f}, P(Superiority) = {p_sup:.4f}")

    # Save model and ATE
    idata.to_netcdf(os.path.join(MODEL_DIR, f"pymc_idata_{outcome_var}.nc"))
    np.save(os.path.join(MODEL_DIR, f"pymc_ate_{outcome_var}.npy"), ate_samples)

    # =========================================================================
    # UPDATE TABLES (with corrected ESS_tail validation for ALL models)
    # =========================================================================
    print("\n--- Updating tables with corrected validation logic ---")

    # Update results table
    results_df = pd.read_csv(os.path.join(TABLE_DIR, "table_pymc_results.csv"))
    mask = results_df["Outcome"] == outcome_label
    if mask.any():
        results_df.loc[mask, "ATE"] = round(ate_mean, 4)
        results_df.loc[mask, "CrI_lo"] = round(ate_ci[0], 4)
        results_df.loc[mask, "CrI_hi"] = round(ate_ci[1], 4)
        results_df.loc[mask, "P_NI"] = round(p_ni, 4)
        results_df.loc[mask, "P_Superiority"] = round(p_sup, 4)
        results_df.loc[mask, "Convergence"] = "OK" if valid else "Issue"

    # Update convergence table
    diag_df = pd.read_csv(os.path.join(TABLE_DIR, "table_pymc_convergence.csv"))
    mask_d = diag_df["outcome"] == outcome_label
    if mask_d.any():
        diag_df.loc[mask_d, "rhat_max"] = float(max_rhat)
        diag_df.loc[mask_d, "ess_bulk_min"] = float(min_bulk)
        diag_df.loc[mask_d, "ess_tail_min"] = float(min_tail)
        diag_df.loc[mask_d, "n_divergent"] = n_divergent
        diag_df.loc[mask_d, "valid"] = valid

    # Also fix validation for ALL models in convergence table using corrected logic
    for idx, row in diag_df.iterrows():
        corrected_valid = (row["rhat_max"] < 1.01
                           and row["ess_bulk_min"] > 400
                           and row["ess_tail_min"] > 400
                           and row["n_divergent"] == 0)
        diag_df.loc[idx, "valid"] = corrected_valid

        # Also update convergence status in results table
        outcome_name = row["outcome"]
        r_mask = results_df["Outcome"] == outcome_name
        if r_mask.any():
            results_df.loc[r_mask, "Convergence"] = (
                "OK" if corrected_valid else "Issue")

    # Save updated tables
    results_df.to_csv(os.path.join(TABLE_DIR, "table_pymc_results.csv"),
                      index=False)
    diag_df.to_csv(os.path.join(TABLE_DIR, "table_pymc_convergence.csv"),
                   index=False)

    print("\nUpdated convergence flags (corrected ESS_tail check):")
    print(f"  {'Outcome':<30s} {'Rhat':>8s} {'ESS_bulk':>10s} {'ESS_tail':>10s} {'Div':>5s} {'Valid':>6s}")
    print("  " + "-" * 72)
    for _, row in diag_df.iterrows():
        print(f"  {row['outcome']:<30s} {row['rhat_max']:>8.4f} "
              f"{row['ess_bulk_min']:>10.0f} {row['ess_tail_min']:>10.0f} "
              f"{int(row['n_divergent']):>5d} {str(row['valid']):>6s}")

    # Rebuild comparison table
    print("\n--- Rebuilding comparison table ---")
    brms_primary = pd.read_csv(
        os.path.join(TABLE_DIR, "table2_primary_results.csv"))
    brms_tier2 = pd.read_csv(
        os.path.join(TABLE_DIR, "table3_tier2_effectiveness.csv"))

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
                "Concordance": ("Concordant" if brms_ni == pymc_ni
                                else "Discordant"),
            })

    comp_df = pd.DataFrame(comparison_rows)
    comp_df.to_csv(os.path.join(TABLE_DIR, "table_brms_vs_pymc.csv"),
                   index=False)

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

    continuous_names = [
        "ODI 3 months", "ODI 12 months",
        "NRS back pain 3 months", "NRS back pain 12 months",
        "NRS leg pain 3 months", "NRS leg pain 12 months",
        "EQ-5D 3 months", "EQ-5D 12 months",
    ]
    continuous = comp_df[comp_df["Outcome"].isin(continuous_names)
                         ].iloc[::-1].reset_index(drop=True)
    binary = comp_df[~comp_df["Outcome"].isin(continuous_names)
                     ].iloc[::-1].reset_index(drop=True)

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
                         capsize=3, linewidth=1.2,
                         label="brms" if i == 0 else "")
            ax1.errorbar(row["PyMC_ATE"], y - 0.4,
                         xerr=[[row["PyMC_ATE"] - row["PyMC_CrI_lo"]],
                               [row["PyMC_CrI_hi"] - row["PyMC_ATE"]]],
                         fmt="s", color="#E66100", markersize=6,
                         capsize=3, linewidth=1.2,
                         label="PyMC" if i == 0 else "")
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
                         capsize=3, linewidth=1.2,
                         label="brms" if i == 0 else "")
            ax2.errorbar(row["PyMC_ATE"], y - 0.4,
                         xerr=[[row["PyMC_ATE"] - row["PyMC_CrI_lo"]],
                               [row["PyMC_CrI_hi"] - row["PyMC_ATE"]]],
                         fmt="s", color="#E66100", markersize=6,
                         capsize=3, linewidth=1.2,
                         label="PyMC" if i == 0 else "")
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

    # Copy to legacy paths
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

    # Update pickle
    pkl_path = os.path.join(MODEL_DIR, "pymc_results.pkl")
    if os.path.exists(pkl_path):
        with open(pkl_path, "rb") as f:
            all_results = pickle.load(f)
        all_results["ate_samples"]["responder_3m"] = ate_samples
        # Update pymc_results list entry
        for i, r in enumerate(all_results.get("pymc_results", [])):
            if r.get("Outcome") == outcome_label:
                all_results["pymc_results"][i].update({
                    "ATE": round(ate_mean, 4),
                    "CrI_lo": round(ate_ci[0], 4),
                    "CrI_hi": round(ate_ci[1], 4),
                    "P_NI": round(p_ni, 4),
                    "P_Superiority": round(p_sup, 4),
                    "Convergence": "OK" if valid else "Issue",
                })
                break
        # Update diagnostics list entry
        for i, d in enumerate(all_results.get("diagnostics", [])):
            if d.get("outcome") == outcome_label:
                all_results["diagnostics"][i] = diag
                break
        with open(pkl_path, "wb") as f:
            pickle.dump(all_results, f)

    elapsed = time.time() - total_start
    print(f"\n{'='*80}")
    print(f"REFIT COMPLETE ({elapsed/60:.1f} minutes)")
    print(f"{'='*80}")
    print(f"\nResponder 3m:")
    print(f"  Before: ESS_tail=223, valid=True (BUG: ESS_tail not checked)")
    print(f"  After:  ESS_tail={min_tail:.0f}, valid={valid}")


if __name__ == "__main__":
    main()
