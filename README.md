# ENDO-LUMBAR: Endoscopic vs Microsurgical Lumbar Discectomy

Analysis code for:

**Endoscopic vs Microsurgical Lumbar Discectomy During Technique Adoption: A Target Trial Emulation**

## Study overview

Registry-based cohort study emulating a target trial, using data from the
Norwegian Registry for Spine Surgery (NORspine). The study compares
endoscopic lumbar discectomy (ELD, n = 124) with microsurgical discectomy
(MSD, n = 297) at Oslo University Hospital during the introduction of the
endoscopic technique (October 2023 through December 2025).

The primary outcome is the Oswestry Disability Index (ODI) at 3 months,
assessed using Bayesian G-computation with a pre-specified non-inferiority
margin of 7 ODI points.

Pre-registered on OSF Registries: <https://osf.io/82qra/>

## Repository contents

- `scripts/analysis/` the statistical analysis pipeline
- `tables/` aggregate result tables backing every number in the publication
- `figures/main/` the three manuscript figures

The repository holds the analysis record. Figure-drawing code and the scripts
that assemble the manuscript document are not included, since they produce no
estimate: every number they render is read from `tables/`, which is published
here in full.

No individual patient data are included, and none can be. NORspine data are
governed by a registry data-sharing agreement. Everything under `data/`,
`results/` and `models/` is patient-level or derived model state and is
excluded by `.gitignore`; the published tables are group-level summaries and
posterior summaries only, and surgery dates appear at month granularity at
finest.

```
scripts/analysis/
├── 00_config.R                    Shared configuration, paths, parameters, priors
├── 01_data_preparation.R          NORspine import, cohort assembly, covariate coding
├── 02_descriptive_table1.R        Baseline characteristics (Table 1)
├── 03_propensity_balance.R        Propensity score estimation, covariate balance
├── 04_primary_analysis.R          Primary outcome: brms ZIB, G-computation ATE
├── 05_mcmc_diagnostics.R          Convergence checks (Rhat, ESS, divergences)
├── 06_secondary_effectiveness.R   Tier 2 outcomes (ODI 12m, NRS, EQ-5D, RTW, ...)
├── 07_perioperative_superiority.R Tier 3 outcomes (day surgery, LOS, complications)
├── 08_descriptive_tier4.R         Operating time, negative control outcome
├── 09_prior_sensitivity.R         Skeptical, reference, diffuse prior comparison
├── 10_missing_data_sensitivity.R  Pattern-mixture, tipping-point, Gaussian mi()
├── 11_model_sensitivity.R         Alternative likelihoods, horseshoe, restricted covariates
├── 12_falsification_evalue.R      E-values, falsification tests, negative control
├── 13_subgroups_causal_forest.R   Bayesian subgroup interactions, causal forest
├── 14_eld_approach_comparison.R   Interlaminar vs transforaminal
├── 15_learning_curve.R            Operating time and ODI vs cumulative case number
├── 16_frequentist_tmle.R          TMLE with SuperLearner (cross-fitted initial Q)
├── 17_refit_12m_clean_eligibility.R  12-month models under date-based eligibility
├── 19_supplement_additions.R      Surgical-level sensitivity, multiplicity, attrition
└── run_all.R                      Master pipeline runner
```

## Analysis pipeline

The pipeline runs sequentially from `00_config.R` through `16_frequentist_tmle.R`.

| Scripts | Stage            | Description                                                           |
|---------|------------------|-----------------------------------------------------------------------|
| 00      | Configuration    | Paths, parameters, NI margins, priors, seed                           |
| 01      | Data preparation | NORspine import, cohort assembly, covariate coding                    |
| 02-03   | Descriptive      | Table 1, propensity score balance                                     |
| 04      | Primary analysis | Bayesian ZIB regression and G-computation for ODI at 3 months         |
| 05      | Diagnostics      | MCMC convergence (Rhat, ESS, divergent transitions)                   |
| 06-07   | Secondary        | Effectiveness and perioperative outcomes (Tier 2 and Tier 3)          |
| 08      | Descriptive      | Operating time and negative control (Tier 4)                          |
| 09-11   | Sensitivity      | Prior, missing data, and model specification sensitivity              |
| 12      | Causal diagnostics | E-value, falsification tests, negative control outcome              |
| 13      | Heterogeneity    | Subgroup interactions and causal forest                               |
| 14-15   | Exploratory      | Endoscopic approach comparison, learning curve                        |
| 16      | Cross-validation | Frequentist TMLE with SuperLearner ensemble                           |

## Software requirements

- **R** 4.5 or later
- **R packages**: `brms`, `cmdstanr`, `rms`, `loo`, `bayesplot`, `posterior`,
  `grf`, `projpred`, `EValue`, `tableone`, `cobalt`, `tidyverse`, `haven`,
  `here`, `scales`, `gt`, `gtsummary`, `patchwork`
- **Stan**: Install `cmdstanr` and run `cmdstanr::install_cmdstan()` for the
  recommended backend.

## Running the pipeline

From the repository root:

```bash
Rscript scripts/analysis/run_all.R
```

All paths are resolved from the repository root via the `here` package, so
the pipeline works from any clone. The only file the pipeline cannot locate
automatically is the raw NORspine SPSS export (see below).

### Raw data

Individual patient data from NORspine cannot be shared publicly. To run
the full pipeline, request access through the NORspine registry and set
the path to the SPSS export before running:

```bash
export ENDO_LUMBAR_RAW_DATA="/path/to/norspine_export.sav"
Rscript scripts/analysis/run_all.R
```

Alternatively, edit `paths$data_raw` directly in `scripts/analysis/00_config.R`.

## Computational notes

- Default MCMC settings: 4 chains, 2000 iterations per chain (1000 warmup),
  seed 20260204. A sensitivity check with 4000 iterations per chain left
  every non-inferiority conclusion unchanged; point estimates for the
  patient-reported outcomes agreed to within 0.09 ODI points. The
  length-of-stay row of `tables/comparison_2k_vs_4k.csv` is not valid: the
  comparison script treats the cumulative ordinal model as continuous.
  The reported length-of-stay effect comes from `07_perioperative_superiority.R`,
  which handles the ordinal model correctly.
- Full pipeline run time on a 14-core workstation: approximately 45-60 minutes.

## Data availability

Individual patient data can be requested through the NORspine registry
application process:
<https://www.kvalitetsregistre.no/registers/norsk-ryggkirurgiregister>.
All analysis scripts are provided here for transparency and reproducibility.

## AI-assisted development

Claude Code (Anthropic; Claude Opus 4.6 for the original analysis and Claude
Opus 5 for the revision) was used as a coding and editing assistant. Its role
covered drafting, debugging and optimising the analysis scripts, auditing the
manuscript against the analysis output, and language editing. No analysis was
accepted on the assistant's word: every reported estimate was regenerated from
the registry data by scripts the authors read and ran. The study question,
design, pre-registered protocol, clinical interpretation and all conclusions
are the authors' own, and the authors take full responsibility for the analysis
and findings. The manuscript Methods section carries the full statement,
including representative prompts.

## License

This analysis code is released under the MIT License. See `LICENSE` for
details. Please cite the associated publication if reusing any part of
this work.
