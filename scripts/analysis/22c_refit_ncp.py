#!/usr/bin/env python3
"""
22c_refit_ncp.py - Refit remaining Bernoulli models with non-centered
parameterization (NCP) for the Cauchy prior.

The Cauchy (StudentT nu=1) prior creates funnel geometry that NUTS struggles
with in centered parameterization. NCP decomposes it as:
  beta_raw ~ Normal(0, 1)
  tau ~ HalfCauchy(1)
  beta = beta_raw * tau

This is mathematically equivalent but avoids the funnel.
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

SEED = 20260204
N_CORES = 14
N_CHAINS = 4
N_DRAWS = 1000
N_TUNE = 2000
TARGET_ACCEPT = 0.99
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
    cols, col_names, std_params = [], [], {}
    for var in CONTINUOUS_COVS:
        vals = df[var].values.astype(float)
        m, s = np.nanmean(vals), max(np.nanstd(vals), 0.01)
        std_params[var] = {"mean": float(m), "std": float(s)}
        z = np.where(np.isnan(vals), 0.0, (vals - m) / s)
        cols.append(z); col_names.append(var + "_z")
    for var in FACTOR_COVS:
        vals = df[var].values
        unique_vals = sorted(pd.Series(vals).dropna().unique())
        if len(unique_vals) <= 2:
            cols.append((vals == unique_vals[1]).astype(float))
            col_names.append(f"{var}_{unique_vals[1]}")
        else:
            for level in unique_vals[1:]:
                cols.append((vals == level).astype(float))
                col_names.append(f"{var}_{level}")
    for var in BINARY_COVS:
        vals = np.where(np.isnan(df[var].values.astype(float)), 0.0,
                        df[var].values.astype(float))
        cols.append(vals); col_names.append(var)

    # Ordinal covariates (integer-coded, treated as linear)
    for var in ORDINAL_COVS:
        vals = np.where(np.isnan(df[var].values.astype(float)), 0.0,
                        df[var].values.astype(float))
        cols.append(vals); col_names.append(var)
    return np.column_stack(cols), col_names, std_params


def fit_bernoulli_ncp(df, X, outcome_var, outcome_label, lower_is_better):
    """Bernoulli model with non-centered Cauchy parameterization."""
    print(f"\n--- Bernoulli NCP: {outcome_label} ---")

    N = len(df)
    treatment = df["treatment_num"].values.astype(float)
    y = df[outcome_var].values.astype(float)
    obs_mask = ~np.isnan(y)
    obs_idx = np.where(obs_mask)[0]
    y_obs = y[obs_mask]
    n_obs = len(y_obs)
    n_covs = X.shape[1]

    print(f"  N={N}, Observed={n_obs}, Covariates={n_covs}")

    with pm.Model() as model:
        intercept = pm.Normal("intercept", mu=0, sigma=5)
        beta_treatment = pm.Normal("beta_treatment", mu=0, sigma=1)

        # Non-centered Cauchy: beta = beta_raw * tau
        # StudentT(nu=1, mu=0, sigma=1) = Cauchy(0, 1)
        # Decomposed as: raw ~ Normal(0, 1), tau ~ HalfCauchy(1)
        # Then beta = raw * tau (per-coefficient local scale)
        beta_raw = pm.Normal("beta_raw", mu=0, sigma=1, shape=n_covs)
        tau = pm.HalfCauchy("tau", beta=1, shape=n_covs)
        beta_covs = pm.Deterministic("beta_covs", beta_raw * tau)

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

    # Diagnostics (check on sampled params, not deterministics)
    rhat = az.rhat(idata)
    max_rhat = max(float(rhat[var].values.max()) for var in rhat.data_vars
                   if var != "beta_covs")  # skip deterministic
    ess_bulk = az.ess(idata, method="bulk")
    ess_tail = az.ess(idata, method="tail")
    min_bulk = min(float(ess_bulk[var].values.min())
                   for var in ess_bulk.data_vars if var != "beta_covs")
    min_tail = min(float(ess_tail[var].values.min())
                   for var in ess_tail.data_vars if var != "beta_covs")

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
    print("REFITTING WITH NON-CENTERED CAUCHY PARAMETERIZATION")
    print("=" * 80)
    total_start = time.time()

    df_full = pd.read_csv(os.path.join(DATA_DIR, "df_disc_full.csv"))
    df_12m = pd.read_csv(os.path.join(DATA_DIR, "df_disc_12m.csv"))
    X_full, _, _ = build_design_matrix(df_full)
    X_12m, _, _ = build_design_matrix(df_12m)

    # Only refit the 5 models still with convergence issues
    to_refit = [
        {"var": "rtw_12m", "label": "RTW 12 months",
         "dataset": "12m", "lower_is_better": False},
        {"var": "analgesic_12m", "label": "Analgesic 12 months",
         "dataset": "12m", "lower_is_better": True},
        {"var": "satisfied_3m", "label": "Satisfaction 3 months",
         "dataset": "full", "lower_is_better": False},
        {"var": "satisfied_12m", "label": "Satisfaction 12 months",
         "dataset": "12m", "lower_is_better": False},
        {"var": "gpe_success_12m", "label": "GPE 12 months",
         "dataset": "12m", "lower_is_better": False},
    ]

    refit_results = {}

    for spec in to_refit:
        df = df_12m if spec["dataset"] == "12m" else df_full
        X = X_12m if spec["dataset"] == "12m" else X_full

        t0 = time.time()
        idata, ate, diag, N, n_obs = fit_bernoulli_ncp(
            df, X, spec["var"], spec["label"], spec["lower_is_better"])

        ate_mean = float(ate.mean())
        ate_ci = [float(np.percentile(ate, 2.5)),
                  float(np.percentile(ate, 97.5))]
        p_ni = float((ate > -NI_MARGINS["binary"]).mean())
        p_sup = float((ate > 0).mean())
        ni_conclusion = "NI demonstrated" if p_ni > NI_THRESHOLD else (
            "Inconclusive" if p_ni >= 0.80 else "Concern")

        print(f"  ATE = {ate_mean:.4f} (95% CrI: [{ate_ci[0]:.4f}, {ate_ci[1]:.4f}])")
        print(f"  P(NI) = {p_ni:.4f}, Conclusion: {ni_conclusion}")
        print(f"  Time: {(time.time()-t0)/60:.1f} min")

        refit_results[spec["var"]] = {
            "result": {
                "Outcome": spec["label"],
                "Family": "BERNOULLI",
                "N_total": N, "N_observed": n_obs, "N_gcomp": N,
                "ATE": round(ate_mean, 4),
                "CrI_lo": round(ate_ci[0], 4),
                "CrI_hi": round(ate_ci[1], 4),
                "P_NI": round(p_ni, 4),
                "P_Superiority": round(p_sup, 4),
                "NI_Margin": NI_MARGINS["binary"],
                "NI_Conclusion": ni_conclusion,
                "Convergence": "OK" if diag["valid"] else "Issue",
            },
            "diag": diag, "ate": ate,
        }

        idata.to_netcdf(os.path.join(MODEL_DIR, f"pymc_idata_{spec['var']}.nc"))
        np.save(os.path.join(MODEL_DIR, f"pymc_ate_{spec['var']}.npy"), ate)
        del idata
        import gc; gc.collect()

    # Update tables
    print("\n--- Updating tables ---")
    results_df = pd.read_csv(os.path.join(TABLE_DIR, "table_pymc_results.csv"))
    diag_df = pd.read_csv(os.path.join(TABLE_DIR, "table_pymc_convergence.csv"))

    for var, data in refit_results.items():
        r = data["result"]
        d = data["diag"]
        mask = results_df["Outcome"] == r["Outcome"]
        if mask.any():
            for col, val in r.items():
                results_df.loc[mask, col] = val
        mask_d = diag_df["outcome"] == d["outcome"]
        if mask_d.any():
            for col, val in d.items():
                diag_df.loc[mask_d, col] = val

    results_df.to_csv(os.path.join(TABLE_DIR, "table_pymc_results.csv"), index=False)
    diag_df.to_csv(os.path.join(TABLE_DIR, "table_pymc_convergence.csv"), index=False)

    # Rebuild comparison table
    brms_primary = pd.read_csv(os.path.join(TABLE_DIR, "table2_primary_results.csv"))
    brms_tier2 = pd.read_csv(os.path.join(TABLE_DIR, "table3_tier2_effectiveness.csv"))
    brms_rows = []
    for _, row in brms_primary.iterrows():
        cri = row["95% CrI"].strip("[]").split(",")
        brms_rows.append({
            "Outcome": row["Outcome"], "brms_ATE": float(row["ATE"]),
            "brms_CrI_lo": float(cri[0].strip()), "brms_CrI_hi": float(cri[1].strip()),
            "brms_P_NI": float(row["P(NI)"]),
            "brms_N_gcomp": int(row["N_ELD"]) + int(row["N_MSD"]),
        })
    for _, row in brms_tier2.iterrows():
        brms_rows.append({
            "Outcome": row["Outcome"], "brms_ATE": float(row["ATE"]),
            "brms_CrI_lo": float(row["CrI_lo"]), "brms_CrI_hi": float(row["CrI_hi"]),
            "brms_P_NI": float(row["P_NI"]), "brms_N_gcomp": int(row["N"]),
        })
    brms_df = pd.DataFrame(brms_rows)

    name_map = {
        "ODI 3 months": "ODI 3 months", "ODI 12 months": "ODI 12 months",
        "NRS back pain 3 months": "NRS back 3 months",
        "NRS back pain 12 months": "NRS back 12 months",
        "NRS leg pain 3 months": "NRS leg 3 months",
        "NRS leg pain 12 months": "NRS leg 12 months",
        "EQ-5D 3 months": "EQ-5D 3 months", "EQ-5D 12 months": "EQ-5D 12 months",
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
    comp_rows = []
    for _, br in brms_df.iterrows():
        pymc_name = name_map.get(br["Outcome"], br["Outcome"])
        pm_match = results_df[results_df["Outcome"] == pymc_name]
        if len(pm_match) > 0:
            pm = pm_match.iloc[0]
            comp_rows.append({
                "Outcome": br["Outcome"],
                "brms_ATE": round(br["brms_ATE"], 4),
                "brms_CrI_lo": round(br["brms_CrI_lo"], 4),
                "brms_CrI_hi": round(br["brms_CrI_hi"], 4),
                "brms_P_NI": round(br["brms_P_NI"], 4),
                "brms_N_gcomp": br["brms_N_gcomp"],
                "PyMC_ATE": pm["ATE"], "PyMC_CrI_lo": pm["CrI_lo"],
                "PyMC_CrI_hi": pm["CrI_hi"], "PyMC_P_NI": pm["P_NI"],
                "PyMC_N_gcomp": pm["N_gcomp"],
                "Concordance": "Concordant" if (br["brms_P_NI"] > NI_THRESHOLD) == (pm["P_NI"] > NI_THRESHOLD) else "Discordant",
            })
    comp_df = pd.DataFrame(comp_rows)
    comp_df.to_csv(os.path.join(TABLE_DIR, "table_brms_vs_pymc.csv"), index=False)

    # Rebuild forest plot
    print("--- Rebuilding forest plot ---")
    plt.rcParams.update({
        "font.family": "sans-serif", "font.size": 10,
        "axes.linewidth": 0.8, "axes.labelsize": 11,
        "xtick.labelsize": 9, "ytick.labelsize": 9,
        "legend.fontsize": 9, "figure.dpi": 300,
    })
    continuous_names = ["ODI 3 months", "ODI 12 months",
                        "NRS back pain 3 months", "NRS back pain 12 months",
                        "NRS leg pain 3 months", "NRS leg pain 12 months",
                        "EQ-5D 3 months", "EQ-5D 12 months"]
    continuous = comp_df[comp_df["Outcome"].isin(continuous_names)].iloc[::-1].reset_index(drop=True)
    binary = comp_df[~comp_df["Outcome"].isin(continuous_names)].iloc[::-1].reset_index(drop=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 10),
                                    gridspec_kw={"width_ratios": [1, 1]})
    for panel_data, ax, title in [(continuous, ax1, "A. Continuous Outcomes"),
                                   (binary, ax2, "B. Binary Outcomes")]:
        if len(panel_data) == 0:
            continue
        y_pos = np.arange(len(panel_data)) * 2.5
        for i, (_, row) in enumerate(panel_data.iterrows()):
            y = y_pos[i]
            ax.errorbar(row["brms_ATE"], y + 0.4,
                        xerr=[[row["brms_ATE"] - row["brms_CrI_lo"]],
                              [row["brms_CrI_hi"] - row["brms_ATE"]]],
                        fmt="o", color="#2166AC", markersize=6,
                        capsize=3, linewidth=1.2, label="brms" if i == 0 else "")
            ax.errorbar(row["PyMC_ATE"], y - 0.4,
                        xerr=[[row["PyMC_ATE"] - row["PyMC_CrI_lo"]],
                              [row["PyMC_CrI_hi"] - row["PyMC_ATE"]]],
                        fmt="s", color="#E66100", markersize=6,
                        capsize=3, linewidth=1.2, label="PyMC" if i == 0 else "")
        ax.axvline(x=0, color="0.30", linestyle="--", linewidth=0.8)
        if "Binary" in title:
            ax.axvline(x=-NI_MARGINS["binary"], color="red", linestyle=":",
                       linewidth=0.8, alpha=0.7)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(panel_data["Outcome"].values)
        ax.set_xlabel("ATE (positive = ELD superior)")
        ax.set_title(title, fontweight="bold", fontsize=12)
        ax.legend(loc="lower right", framealpha=0.9)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "fig_brms_vs_pymc_forest.png"),
                dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # Copy outputs
    for f in ["table_pymc_results.csv", "table_brms_vs_pymc.csv",
              "table_pymc_convergence.csv"]:
        try: shutil.copy2(os.path.join(TABLE_DIR, f), os.path.join(LEGACY_TABLE_DIR, f))
        except: pass
    try: shutil.copy2(os.path.join(FIG_DIR, "fig_brms_vs_pymc_forest.png"),
                      os.path.join(LEGACY_FIG_DIR, "fig_brms_vs_pymc_forest.png"))
    except: pass

    elapsed = time.time() - total_start
    print(f"\n{'='*80}")
    print(f"NCP REFIT COMPLETE ({elapsed/60:.1f} minutes)")
    print(f"{'='*80}")

    print(f"\n{'Outcome':<30s} {'Div_before':>10s} {'Div_NCP':>8s} {'Rhat':>8s} {'ESS_tail':>10s} {'Status':>8s}")
    print("-" * 80)
    old_divs = {"RTW 12 months": 10, "Analgesic 12 months": 146,
                "Satisfaction 3 months": 50, "Satisfaction 12 months": 114,
                "GPE 12 months": 62}
    for var, data in refit_results.items():
        d = data["diag"]
        old_d = old_divs.get(d["outcome"], "?")
        status = "FIXED" if d["valid"] else "PERSIST"
        print(f"{d['outcome']:<30s} {old_d:>10} {d['n_divergent']:>8d} "
              f"{d['rhat_max']:>8.4f} {d['ess_tail_min']:>10.0f} {status:>8s}")


if __name__ == "__main__":
    main()
