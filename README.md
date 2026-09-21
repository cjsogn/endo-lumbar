# ENDO-LUMBAR analysis code

Statistical analysis of full-endoscopic versus microsurgical lumbar discectomy during institutional adoption, using the Norwegian Registry for Spine Surgery (NORspine).

This package contains data preparation, Bayesian g-computation, targeted maximum likelihood estimation (TMLE), sensitivity analyses and exploratory analyses. The registry extract requires authorised access and is not included or downloaded by the code.

The [published SAP on OSF](https://osf.io/82qra/) is the prespecification reference. [SAP_DEVIATIONS.md](SAP_DEVIATIONS.md) distinguishes registered methods, implementation departures and additional analyses. Not every registered component was completed.

## Getting started

Read [REVIEW_GUIDE.md](REVIEW_GUIDE.md) for the recommended order of statistical review. Software requirements and recorded versions are in [DEPENDENCIES.md](DEPENDENCIES.md).

Use R, Python 3 and CmdStan with a working C++ toolchain. Set the location of the frozen SPSS export and a separate output directory, then run:

```sh
export ENDO_LUMBAR_RAW_DATA="/absolute/path/to/authorised_registry_export.sav"
export ENDO_WORK_DIR="/absolute/path/to/private_analysis_workspace"
Rscript scripts/analysis/run_all.R
```

Both locations should be outside the repository. The runner rejects output directories inside the checkout, including paths through symlinks. Cohort cutoffs refer to the frozen study extract, not a live registry query. The supplied EQ-5D index values are used without rescoring. Their valuation tariff is not documented in the extract.

The entry point works from any working directory. To run selected stages, append their exact filenames as listed in `run_all.R`. Their inputs must already have been produced by preceding stages. A failed stage stops execution. Existing `brms` fits are reused only through `file_refit="on_change"`, which checks the model code, data and priors. Selection models and TMLE checks are rerun.

Data, fitted models, numerical summaries, diagnostics and logs are written under `ENDO_WORK_DIR`. These outputs can include individual records and identifiers and must remain in an authorised study workspace.

## Code map

| Source in `scripts/analysis/` | Purpose |
| --- | --- |
| `00_config.R` | Shared covariates, priors, effect definitions, seed and diagnostic functions |
| `01_import_registry.R`, `01_prepare_data.R` | Registry coding, eligibility, baseline imputation, standardisation and calendar basis |
| `02_primary_jobs.R`, `02_fit_primary.R`, `02_summarize_primary.R` | Primary ODI and perioperative specifications, fitting and posterior summaries |
| `03_model_jobs.R`, `04_run_calendar_models.R` | Secondary outcomes, operating time and negative-control models |
| `05_primary_derived_sensitivities.R` | E-values, pattern-mixture shifts and tipping-point calculations |
| `06_primary_model_sensitivities.R` | Alternative priors, likelihoods, populations and adjustment sets |
| `07_tmle_matched_checks.R` | Observed-case, outcome-imputed and censoring-weighted TMLE |
| `08_selection_model.R` | Joint outcome/response selection-model sensitivity |
| `09_model_and_population_audit.R` | Matched populations, targeting diagnostics, overlap and field availability |
| `11_predictive_diagnostics.R` | Posterior predictive statistics, PSIS-LOO and prior/posterior dispersion |
| `12_finalize_result_tables.R` | Endpoint summaries, conditional-association checks and sampler diagnostics |
| `13_descriptive_analyses.R` | Baseline characteristics, institutional case sequence, inferred corridor and recorded reoperations |
| `14_exploratory_analyses.R` | Subgroup interactions, causal forest and propensity weighting |
| `04_model_scheduler.py`, `run_all.R` | Parallel scheduling and pipeline execution |

The scheduler uses up to 14 CPU cores in independent processes, allocated as 4 + 4 + 4 + 2. Each Bayesian fit retains four chains. Compilation is sequential to limit memory use. The forest uses 14 threads and the weighting bootstrap uses 14 local socket workers.

## Models and effect definitions

Primary ODI, adjusted secondary and perioperative outcomes, operating time and negative controls use the same 27 baseline covariates, two natural cubic calendar terms (`ns(df=2)`) and both treatment-by-spline interactions. Alternative primary models retain the calendar interactions, including explicitly reduced baseline adjustment sets.

Continuous baseline covariates are standardised in the relevant cohort. Twelve-month models retain the calendar basis constructed from the full cohort and restandardise continuous baseline covariates within the date-eligible cohort.

The primary contrast is MSD minus ELD in ODI units. Positive values favour ELD. Non-inferiority requires a posterior probability greater than 0.95 that the contrast exceeds the negative margin. The registered ODI margin is 7 points, with additional sensitivity margins of 5 and 3 points. Other contrasts use the endpoint-specific direction and scale in the code.

Bounded and binary models fit observed outcomes and standardise predictions over their stated target cohorts. Gaussian models include latent missing outcomes through `mi()`. Deterministic baseline median/mode imputation is separate from outcome multiple imputation.

Length of stay uses a cumulative-logit model. Its effect summary is the conditional log odds ratio for shorter stay averaged over the calendar dates of procedures with observed stay. Exponentiation gives a geometric mean conditional odds ratio, not a marginal common odds ratio or a difference in mean nights.

## TMLE and missing outcomes

Observed-case TMLE covers 17 effectiveness outcomes, day surgery, postoperative stay and patient-reported complications. All 17 effectiveness outcomes and day surgery also have outcome multiple-imputation checks. All eight twelve-month effectiveness outcomes have censoring-weighted checks. A day-surgery model without calendar adjustment is an additional diagnostic, stored under the output key `day_surgery_original`.

TMLE includes the calendar basis without explicit treatment-by-spline product columns. The observed-case learner library contains GLM, glmnet, ranger and xgboost. Imputation and censoring-weighted checks use GLM and glmnet. `cvQinit=TRUE` and learner cross-validation do not constitute full cross-fitting of both nuisance models. Required learners are not replaced automatically if fitting fails.

Outcome imputation uses 20 datasets, 10 iterations and seed 20260204. Every dataset is analysed, and estimates and variances are combined with Rubin pooling through `mice::pool.scalar`. No single imputation supplies the final estimate.

Censoring weights use a common model of twelve-month ODI response as a proxy for follow-up attendance. They are stabilised by treatment arm, use a probability floor of 0.01 and are capped at the pooled 99th percentile. This model does not separately describe item-specific nonresponse.

## Exploratory analyses

Four subgroup models examine age, baseline ODI severity, leg symptom duration and number of operated levels. They use zero-inflated beta likelihoods, Normal(0,1) mean-coefficient priors, the common calendar basis and both calendar interactions. The multilevel indicator is represented once. Every retained posterior draw contributes to standardisation over the full subgroup. Differences between subgroup-average ODI effects are distinct from conditional model-scale interaction coefficients. They do not support separate formal subgroup non-inferiority decisions.

The causal forest uses baseline predictors and calendar basis columns, with treatment supplied separately. It uses honesty, 4,000 trees, tuning of all supported parameters, and out-of-bag nuisance estimates and conditional predictions. Outputs include average-score inference, calibration, linear projection, split importance and conditional prediction summaries. Split importance is not a causal modifier ranking.

Propensity weighting uses a logistic treatment model with all baseline covariates and calendar basis columns. Stabilised inverse-probability weights are capped at the pooled 99th percentile. Overlap weights are also evaluated. Outcome means use observed ODI only. Percentile intervals use 2,000 procedure bootstrap samples, refitting the propensity model and capping in every sample. Individually seeded resamples are processed in batches. Failed fits or nonfinite estimates prevent final interval output.

The weighting and forest analyses do not adjust for outcome nonresponse. Capping and overlap weighting also change the weighted population. Effective sample sizes, weight distributions and balance are saved for full and observed-outcome cohorts.

## Interpretation and reproducibility

Twelve-month effectiveness analyses are exploratory. Case order describes the registry-recorded institutional sequence and does not identify individual surgeons' experience. Corridor is inferred from approach codes. The same-level reoperation flag is a proxy and cannot adjudicate residual versus recurrent disc.

Each script states its inputs and outputs. Existing numerical output keys are retained so downstream stages can locate their inputs. The runner records installed package versions. Seeds are fixed, but exact binary equality across software versions, compilers and platforms is not assumed.

Run `python3 scripts/check_release.py` to verify the file allowlist and source checks. Use `--staged` to check the Git index before committing. `RELEASE_FILES.txt` lists the permitted contents.
