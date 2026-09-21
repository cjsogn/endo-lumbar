# Statistical review guide

Start with the primary effect definition, cohort construction and missing-data handling. The files below follow their dependency order. Keep the package together so shared functions and model specifications remain available.

| Priority | Files in `scripts/analysis/` | Review focus |
| --- | --- | --- |
| 1 | `00_config.R` | Covariates, priors, calendar interactions, contrast direction, g-computation and non-inferiority rule |
| 2 | `01_import_registry.R`, `01_prepare_data.R` | Registry coding, eligibility, outcomes, imputation and calendar basis |
| 3 | `02_primary_jobs.R`, `02_fit_primary.R`, `02_summarize_primary.R` | Primary ODI and perioperative likelihoods, target populations, posterior summaries and stricter margins |
| 4 | `03_model_jobs.R`, `04_run_calendar_models.R` | Secondary specifications, common calendar adjustment and endpoint-specific scales |
| 5 | `07_tmle_matched_checks.R`, `09_model_and_population_audit.R` | Learners, observed versus standardisation populations, pooling, censoring weights and targeting diagnostics |
| 6 | `05_primary_derived_sensitivities.R`, `06_primary_model_sensitivities.R`, `08_selection_model.R` | Confounding and missingness sensitivity, priors, likelihoods and adjustment sets |
| 7 | `11_predictive_diagnostics.R`, `12_finalize_result_tables.R` | Sampling checks, predictive assessment, influence diagnostics and summaries |
| 8 | `13_descriptive_analyses.R`, `14_exploratory_analyses.R` | Clinical coding, subgroups, causal forest and bootstrap weighting |

`run_all.R` and `04_model_scheduler.py` manage execution. They are dependencies and lower-priority files for methodological review. `scripts/check_release.py` checks package contents and does not analyse registry data.

## Questions to assess

- Do the registry filters and baseline coding identify the stated target population?
- Are model families, priors, calendar terms and interactions appropriate for the outcomes and available information?
- Are contrast direction, rescaling, non-inferiority margins and the length-of-stay estimand implemented consistently?
- Are missing-data assumptions and standardisation populations stated accurately, particularly at twelve months?
- Does every outcome-imputation analysis include all datasets and pool uncertainty correctly?
- Do the TMLE and Bayesian comparisons distinguish population differences from differences in assumptions?
- Are sensitivity analyses and predictive diagnostics sufficient to assess the statistical conclusions?
- Are exploratory analyses and departures from the published SAP identified accurately?

## Materials and execution

Read [README.md](README.md), [DEPENDENCIES.md](DEPENDENCIES.md) and [SAP_DEVIATIONS.md](SAP_DEVIATIONS.md). The published SAP is available through [OSF](https://osf.io/82qra/).

Code review does not require patient data. Reproducing estimates requires the authorised frozen registry export and a private output directory. This package contains no registry data, fitted objects or numerical results. Execution instructions are in the README. Selected stages require outputs from preceding stages.
