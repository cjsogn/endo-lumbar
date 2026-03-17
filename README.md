# ENDO-LUMBAR: Endoscopic vs Microsurgical Lumbar Discectomy

Analysis code for:

**Endoscopic vs Microsurgical Discectomy for Lumbar Disc Herniation During Technique Adoption: A Bayesian Non-Inferiority Study**

## Study overview

Registry-based cohort study emulating a target trial, using data from the Norwegian
Registry for Spine Surgery (NORspine). Compares endoscopic lumbar discectomy (ELD,
n=124) with microsurgical discectomy (MSD, n=297) at Oslo University Hospital during
the introduction of the endoscopic technique (October 2023 to December 2025).

The primary outcome was the Oswestry Disability Index (ODI) at 3 months, assessed
using Bayesian G-computation with a pre-specified non-inferiority margin of 7 ODI
points.

Pre-registered on OSF Registries: https://osf.io/82qra/

## Repository structure

```
scripts/analysis/
├── 00_config.R              # Shared configuration, paths, parameters, seed
├── 01_data_preparation.R    # NORspine import, cohort assembly, covariate coding
├── 02_descriptive_table1.R  # Baseline characteristics (Table 1)
├── 03_propensity_balance.R  # Propensity score estimation, covariate balance
├── 04_primary_analysis.R    # Primary outcome: brms ZIB, G-computation ATE
├── 05_mcmc_diagnostics.R    # Convergence checks (Rhat, ESS, divergences)
├── 06_secondary_effectiveness.R  # Tier 2 outcomes (NRS, EQ-5D, responder, RTW)
├── 07_perioperative_superiority.R # Tier 3 outcomes (day surgery, LOS, complications)
├── 08_descriptive_tier4.R   # Operating time, negative control outcome
├── 09_prior_sensitivity.R   # Skeptical, reference, diffuse prior comparison
├── 10_missing_data_sensitivity.R  # Pattern-mixture, tipping-point, selection model
├── 11_model_sensitivity.R   # Alternative likelihoods, horseshoe, restricted covariates
├── 12_falsification_evalue.R # E-values, falsification tests, negative control
├── 13_subgroups_causal_forest.R   # Bayesian subgroup interactions, causal forest
├── 14_stenosis_exploratory.R      # Concomitant stenosis subgroup
├── 15_tables_figures.R      # Summary tables from model results
├── 16_covariate_influence.R # Projpred variable selection
├── 17_eld_approach_comparison.R   # Interlaminar vs transforaminal
├── 18_learning_curve.R      # Operating time and ODI vs case number
├── 19_odi_spider_plot.R     # ODI domain profiles
├── 20_publication_figures.R # Main and supplementary figures
├── 21_frequentist_tmle.R    # TMLE with SuperLearner cross-validation
├── 23b_update_figures_2_3.R # Regenerate primary figures from brms results
├── 24_causal_diagram.R      # DAG visualization
├── 25_update_causal_outputs.R     # Causal inference summary tables
└── run_all.R                # Master pipeline runner
```

## Analysis pipeline

The pipeline runs sequentially from scripts 00 through 25:

| Scripts | Stage | Description |
|---------|-------|-------------|
| 00 | Configuration | Paths, parameters, NI margins, seed |
| 01 | Data preparation | NORspine import, cohort assembly, covariate coding |
| 02-03 | Descriptive | Table 1, propensity score balance |
| 04-08 | Primary analysis | brms Bayesian models: primary, secondary, perioperative, descriptive |
| 09-11 | Sensitivity | Prior, missing data, and model specification sensitivity |
| 12 | Causal diagnostics | E-value, falsification tests, negative control |
| 13 | Heterogeneity | Subgroup interactions, causal forests |
| 14 | Exploratory | Concomitant stenosis subgroup |
| 15-20 | Outputs | Tables, figures, publication formatting |
| 21 | Cross-validation | Frequentist TMLE with SuperLearner |
| 23b-25 | Integration | Update figures and causal outputs |

## Software

- **R 4.5**: brms (cmdstanr backend), grf, projpred, EValue, tmle3, sl3, ggplot2

## Data availability

Individual patient data can be requested through the NORspine registry application
process. All analysis scripts are provided in this repository for transparency and
reproducibility.

## License

Analysis code is provided for transparency and reproducibility. Please cite the
associated publication if reusing any part of this work.
