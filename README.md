# ENDO-LUMBAR analysis scripts

Analysis code for a target trial emulation comparing full-endoscopic (ELD) and microsurgical (MSD) lumbar discectomy during institutional adoption, using data from the Norwegian Registry for Spine Surgery (NORspine). The statistical analysis plan is published on OSF: <https://osf.io/82qra/>.

The scripts cover data preparation, Bayesian g-computation, targeted maximum likelihood estimation (TMLE), sensitivity analyses and exploratory analyses. Registry data are not included. Access to the frozen registry extract requires authorisation from NORspine.

## Requirements

R 4.5.3, CmdStan 2.38.0 with a working C++ toolchain, and Python 3 (standard library only, used for job scheduling). The analyses were run with the following R package versions.

| Package | Version |
| --- | --- |
| brms | 2.23.0 |
| posterior | 1.7.0 |
| haven | 2.5.5 |
| cmdstanr | 0.9.0 |
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |
| tibble | 3.3.1 |
| tableone | 0.13.2 |
| loo | 2.9.0 |
| tmle | 2.1.1 |
| SuperLearner | 2.0.40 |
| mice | 3.19.0 |
| glmnet | 4.1.10 |
| ranger | 0.18.0 |
| xgboost | 3.2.1.1 |
| grf | 2.6.1 |

Other versions may work. Seeds are fixed in `00_config.R`, but new sampling runs are not guaranteed to be bitwise identical across software versions, compilers or platforms.

## Running the analysis

```sh
export ENDO_LUMBAR_RAW_DATA="/absolute/path/to/authorised_registry_export.sav"
export ENDO_WORK_DIR="/absolute/path/to/private_analysis_workspace"
Rscript scripts/analysis/run_all.R
```

Both paths must be outside the repository. Data, fitted models, numerical results, diagnostics and logs are written under `ENDO_WORK_DIR`. These outputs contain individual patient records and must remain in an authorised study workspace. Cohort cutoffs refer to the frozen study extract, not a live registry query.

To run selected stages, append their exact filenames as listed in `run_all.R`. Each stage needs the outputs of the stages before it, and a failed stage stops execution. Existing `brms` fits are reused only when the model code, data and priors are unchanged (`file_refit = "on_change"`).

The scheduler runs Bayesian fits on up to 14 CPU cores in independent processes, with four chains per fit.

### Follow-up opportunity sensitivity

`16_followup_eligibility.R` runs only when a reference date and its basis are supplied:

```sh
export ENDO_FOLLOWUP_REFERENCE="2026-03-03"
export ENDO_FOLLOWUP_REFERENCE_BASIS="SPSS file creation date"
Rscript scripts/analysis/run_all.R 16_followup_eligibility.R
```

The SPSS file header records creation on 3 March 2026, the date NORspine supplied the extract. The routine excludes operations performed after the date three calendar months before the reference date. It reuses the primary posterior only if no observed ODI outcome is excluded and every model-input row matches the retained fit. Otherwise it stops, because a refit would be required.

## Scripts

| Source in `scripts/analysis/` | Purpose |
| --- | --- |
| `00_config.R` | Shared covariates, priors, effect definitions, seed and diagnostic functions |
| `01_import_registry.R`, `01_prepare_data.R` | Registry coding, eligibility, baseline imputation, standardisation and calendar basis |
| `02_primary_jobs.R`, `02_fit_primary.R`, `02_summarize_primary.R` | Primary ODI and perioperative specifications, fitting and posterior summaries |
| `03_model_jobs.R`, `04_run_calendar_models.R` | Secondary outcomes, operating time and negative-control models |
| `05_primary_derived_sensitivities.R` | E-values, pattern-mixture shifts and tipping-point calculations |
| `06_primary_model_sensitivities.R` | Alternative priors, likelihoods, populations and adjustment sets |
| `07_tmle_matched_checks.R` | Observed-case, outcome-imputed and censoring-weighted TMLE |
| `08_selection_model.R` | Joint outcome and response selection-model sensitivity |
| `09_model_and_population_audit.R` | Matched populations, targeting diagnostics, overlap and field availability |
| `11_predictive_diagnostics.R` | Posterior predictive statistics, PSIS-LOO and prior/posterior dispersion |
| `12_finalize_result_tables.R` | Endpoint summaries, conditional-association checks and sampler diagnostics |
| `13_descriptive_analyses.R` | Baseline characteristics, institutional case sequence, inferred corridor and recorded reoperations |
| `14_exploratory_analyses.R` | Subgroup interactions, causal forest and propensity weighting |
| `15_discharge_descriptives.R` | Planned versus actual discharge and historical MSD discharge rates |
| `16_followup_eligibility.R` | Primary ODI standardisation after restricting follow-up opportunity |
| `04_model_scheduler.py`, `run_all.R` | Parallel scheduling and pipeline execution |

Each script states its inputs and outputs in its header.

## Effect direction

The primary contrast is MSD minus ELD in Oswestry Disability Index (ODI) points, so positive values favour ELD. Non-inferiority requires a posterior probability above 0.95 that the contrast exceeds the registered margin of −7 points. Margins of 5 and 3 points are evaluated as sensitivity analyses. Other contrasts use the endpoint-specific direction and scale set where each model is summarised.

## License

MIT. See [LICENSE](LICENSE).
