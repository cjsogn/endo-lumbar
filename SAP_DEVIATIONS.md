# Relationship to the published SAP

The [published OSF statistical analysis plan](https://osf.io/82qra/) is the prespecification reference. This document records the methodological status of the implemented analyses. It is not a replacement SAP or a claim that all registered components were completed.

## Registered components and implementation departures

| Topic | Implemented analysis and departure |
| --- | --- |
| Primary outcome and margin | Three-month ODI retains the registered 7-point margin and posterior non-inferiority threshold greater than 0.95. |
| Likelihood and priors | The registered primary likelihood was Gaussian. The primary implementation uses zero-inflated beta for bounded ODI and NRS. Treatment-prior standard deviations of 0.5, 1 and 2 on the logit scale are not equivalent to the registered 4, 10 and 25 on the ODI scale. A Gaussian sensitivity is included. Exact priors are in `00_config.R` and the model builders. |
| Missing outcomes | Bounded and binary models use observed-outcome likelihoods and standardise over the stated cohort. This differs from the registered latent-outcome approach. Gaussian models retain latent outcomes through `mi()`. |
| Calendar time | The registered spline and treatment interaction are included in all adjusted Bayesian mean models. The basis has two natural cubic terms, corresponding to three knots. Empirical knot positions were not registered. Decision chronology is stated below. |
| Pain scale | The registry NRS scale is 0 to 10. The registered 10-point margin on a 0 to 100 scale is proportionally rescaled to 1 point. |
| Binary priors | Endpoint-specific Normal priors are implementation choices. The published SAP supplies Gaussian-model priors and does not specify a separate Cauchy prior for binary covariates. |
| Study period | The cohort begins in October 2023, when both techniques were available, rather than January 2022. The supplied extract ends in December 2025 rather than February 2026. |
| Baseline ODI eligibility | Procedures with missing baseline ODI were retained with deterministic fill, although baseline ODI was required by the SAP. |
| Operated level | Main models omit operated level. A sensitivity adds the source level fields. |
| Surgical side | A dedicated operative-side field has not been established in the extract, so this registered covariate is omitted. Complication fields referring to wrong side do not supply operative laterality. |
| Baseline missingness | Median/mode fill is used. The registered joint-covariate missingness sensitivity was not completed. Outcome imputation does not propagate uncertainty in filled baseline values. |
| Projection predictive selection | The registered projection-predictive analysis was not completed. The full adjustment set and restricted-set sensitivity do not replace it. |
| Separate stenosis population | This package addresses disc herniation. The separately registered non-disc stenosis comparison and across-indication interaction are outside its implemented scope. |
| Complication testing | The rate rule is applied to the pooled observed rate. The SAP did not specify pooled versus arm-specific adjudication. |

## Additional analyses and implementation limits

TMLE is an additional, post-hoc estimator check rather than a method specified in the published SAP. Its outcome multiple-imputation checks also differ from the registered Bayesian missing-outcome strategy. Every imputation contributes to pooling. The twelve-month censoring-weighted checks use ODI response as a common follow-up proxy, without separate item-specific response models.

The 5-point and 3-point ODI margins and matched-population Bayesian/TMLE comparisons are additional analyses. They do not replace the registered 7-point primary margin. Stricter-margin probabilities are calculated from the same primary posterior.

Subgroup and causal-forest methods were registered as exploratory. Propensity weighting is post hoc. The weighting implementation was reconstructed from the retained method description because the earlier generating script was unavailable. Its variance calculation uses a procedure bootstrap with propensity estimation and capping repeated in each sample. It is not an exact rerun of the unavailable code.

The fixed-nu Student-t sensitivity uses the registered nu of 5. The transformed-beta sensitivity applies the inverse transformation when returning to ODI units. These specifications include corrections identified during analysis checking.

Excluding January to December 2022 cannot assess early endoscopic adoption because that period predates endoscopic availability and contributes no included procedures. Descriptive case-order analyses do not estimate time-varying adjusted treatment effects. The observed-outcome Bayesian check standardises the primary fitted model over responders rather than fitting a separate outcome model. A Docker environment is not supplied.

## Decision chronology

Treatment identities were known to the author making analytic decisions. The author reports that initial departures arose after data receipt during diagnostic assessment and other analytic considerations, before reviewing comparative treatment results. This is a retrospective account. Outcome information was available during diagnostics, and there is no contemporaneous item-by-item decision log establishing masking.

Calendar terms were initially omitted because of perceived absence of a learning-curve effect and concern about model complexity. That rationale does not exclude secular confounding. Their inclusion, the common calendar-interaction specification, stricter margins and matched-population investigations were decided after peer review. The common interaction specification for secondary models was fixed before those refits, with earlier estimates already known. It was not selected separately by endpoint. Calendar updates to exploratory models and the weighting reconstruction were also specified before their refits. These decisions are not retrospectively described as preregistered.
