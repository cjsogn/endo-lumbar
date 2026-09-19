# ENDO-LUMBAR analysis code

Statistical code for the Neurospine major revision comparing full-endoscopic and microsurgical lumbar discectomy during institutional adoption.

This release contains analysis source and execution instructions. Patient data, fitted models, numerical results, figures, plotting scripts and manuscript generators are excluded. The NORspine export requires authorised registry access. It is not included or downloaded by this code.

The published statistical analysis plan is registered on [OSF](https://osf.io/82qra/). This code implements the revised analysis described in the manuscript and its reported deviations from that plan. It does not represent every historical exploratory analysis as an updated analysis. Earlier code remains in Git history.

## Running the analyses

Use R, Python 3 and CmdStan with a working C++ toolchain. Install the R packages listed in `DEPENDENCIES.md`. No alternative backend or reduced learner library is selected automatically when a dependency is missing.

Use the study's frozen SPSS registry extract. The cohort cutoffs in the import stage do not define a live registry query. Set the private export and a private output directory, then run:

```sh
export ENDO_LUMBAR_RAW_DATA="/absolute/path/to/authorised_registry_export.sav"
export ENDO_WORK_DIR="/absolute/path/to/private_analysis_workspace"
Rscript scripts/analysis/run_all.R
```

Both locations should be outside the repository. The runner rejects output directories inside the checkout, including paths through symlinks. It uses the supplied EQ-5D index scores without rescoring. Their valuation tariff is not documented in the extract.

The entry point works from any current directory. To run selected stages, append their exact filenames as listed in `run_all.R`. Dependencies produced by earlier stages must already exist. A failed stage stops the run. Existing `brms` fits are reused only through `file_refit="on_change"`, which checks the model code, data and priors. Selection models and TMLE checks are rerun.

Generated data, models, numerical summaries, diagnostics and logs are written beneath `ENDO_WORK_DIR`. Treat this entire directory as protected study material. Numerical output can contain individual records and identifiers even when the final manuscript uses aggregate summaries.

## Analysis sequence

| Source | Function |
| --- | --- |
| `01_import_registry.R`, `01_prepare_data.R` | Import the SPSS fields, define the overlap and date-eligible cohorts, apply baseline median/mode imputation, standardise covariates and construct the calendar basis |
| `02_primary_jobs.R`, `02_fit_primary.R`, `02_summarize_primary.R` | Fit and summarise primary ODI and perioperative models |
| `03_model_jobs.R`, `04_run_calendar_models.R` | Fit secondary outcomes and negative controls |
| `05_primary_derived_sensitivities.R` | Margin-referenced E-values, pattern-mixture shifts and tipping-point calculations |
| `06_primary_model_sensitivities.R` | Alternative priors, likelihoods and adjustment sets |
| `07_tmle_matched_checks.R` | TMLE checks, outcome imputation and censoring weighting |
| `08_selection_model.R` | Selection-model sensitivity grid |
| `09_model_and_population_audit.R` | Estimator population checks, targeting diagnostics and overlap summaries |
| `11_predictive_diagnostics.R` | Numerical posterior predictive and PSIS-LOO diagnostics |
| `12_finalize_result_tables.R` | Endpoint summaries, conditional-association checks and sampler diagnostics |
| `13_descriptive_analyses.R` | Baseline characteristics, institutional case sequence, inferred corridor and recorded reoperations |

`00_config.R` contains the covariate names, priors, contrast definitions and seed. `04_model_scheduler.py` runs independent R processes with up to 14 CPU cores, allocated as 4 + 4 + 4 + 2. Each Bayesian fit retains four chains. Compilation is sequential to limit memory use.

## Specifications and interpretation

- Primary ODI uses the 27 baseline covariates, two natural cubic spline terms for calendar time (`ns(df=2)`) and both treatment-by-spline interactions. Other adjusted outcomes use the same baseline covariates and calendar spline main effects. Primary-model sensitivity analyses retain the interactions, including their explicitly specified reduced adjustment sets.
- Continuous baseline covariates are standardised in the relevant cohort. The twelve-month models retain the calendar basis constructed from the full cohort and restandardise the continuous baseline covariates within the date-eligible cohort.
- The primary contrast is MSD minus ELD in ODI units. Positive contrasts favour ELD. Non-inferiority is assessed against the registered 7-point margin, with 5-point and 3-point sensitivity margins. Other contrasts use the endpoint-specific direction and scale in the code. Length of stay uses a cumulative-logit model and a log odds ratio for shorter stay.
- Zero-inflated beta and binary models fit observed outcomes and standardise predictions over their stated target cohorts. Gaussian models include latent missing outcomes through `mi()`. Baseline median/mode imputation is distinct from the multiple-imputation sensitivity analyses.
- TMLE uses the calendar basis, without explicit treatment-by-spline product columns. The observed-case library contains GLM, glmnet, ranger and xgboost. Imputation and censoring-weighted checks use the specified smaller GLM/glmnet library. `cvQinit=TRUE` and learner cross-validation do not constitute full cross-fitting of both nuisance models.
- Each imputation check fits all 20 imputed datasets and combines estimates and variances with `mice::pool.scalar`. No imputation is selected as the final estimate. The censoring-weighted check uses its specified probability floor and weight truncation.
- Twelve-month effectiveness is exploratory. Case order describes the registry-recorded institutional sequence, not a complete surgeon-specific learning curve. Corridor is inferred from approach codes. The same-level reoperation flag is a proxy and cannot adjudicate residual versus recurrent disc.

This release reproduces the current numerical analyses. Historical submitted-versus-revised comparison columns refer to previously saved submitted analyses and are not reconstructed by the default revised pipeline. Figure rendering, manuscript formatting and journal response preparation remain outside this repository.

## Release checks

Run `python3 scripts/check_release.py` to verify the explicit file allowlist and excluded file types. Run it with `--staged` before committing. `RELEASE_FILES.txt` lists the permitted release files. These checks supplement review of the files being published. They are not a substitute for data-access controls.
