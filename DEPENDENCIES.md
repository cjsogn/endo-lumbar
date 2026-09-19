# Dependencies

Versions installed during verification of this release. These are recorded versions, not a claim that other versions are incompatible. R and CmdStan versions are recorded below.

| R package | Version |
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

Install the packages explicitly before running the pipeline. CmdStan and cmdstanr must be configured separately. Python 3 uses only its standard library. The runner records installed package versions in the private output directory. Seeds are fixed in `00_config.R`, but bitwise equality of new sampling runs is not guaranteed across software versions, compilers or platforms.
R version: 4.5.3 (2026-03-11). CmdStan version: 2.38.0.
