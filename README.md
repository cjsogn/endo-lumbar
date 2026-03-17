# ENDO-LUMBAR: Endoscopic vs Microsurgical Lumbar Discectomy

Analysis code and manuscript materials for:

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
ENDO_LUMBAR/
├── scripts/
│   ├── analysis/              # Core analysis pipeline (34 scripts)
│   │   ├── 00_config.R        # Shared configuration and parameters
│   │   ├── 01-21_*.R          # R analysis pipeline
│   │   ├── 22_pymc_bayesian.py # Primary PyMC Bayesian G-computation
│   │   ├── 22a-d_*.R/py       # Data export and model refitting
│   │   ├── 23-25_*.R/py       # Output assembly and figure updates
│   │   └── run_all.R          # Master pipeline runner
│   └── manuscript/            # Document generation scripts (12 scripts)
│       └── generate_*.py/R    # Reproducible .docx generation from CSVs
├── manuscript/
│   ├── main/                  # Main text, tables, and figure legends
│   ├── supplement/            # eAppendix, eTables
│   └── submission/            # Cover letter, checklists, RIS references
├── figures/
│   ├── main/                  # Figures 1-3
│   └── supplement/            # eFigures 1-22
├── tables/                    # Analysis output CSVs (39 files)
├── results/                   # Model objects (.rds, 17 files)
├── diagnostics/               # MCMC diagnostic plots
├── complete_case/             # Complete-case sensitivity tables
├── supplementary_followup_eligible/  # Follow-up eligible sensitivity
│   ├── scripts/
│   ├── results/
│   ├── data/
│   └── models/
└── data/                      # Not included (see data/README.md)
```

## Analysis pipeline

The analysis pipeline runs sequentially from scripts 00 through 25. Key stages:

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
| 22-22d | PyMC models | Primary Bayesian G-computation (Python) |
| 23-25 | Integration | Merge PyMC results, update figures and outputs |

## Software

- **R 4.5**: brms, grf, projpred, EValue, tmle3, sl3, ggplot2
- **Python 3.12**: PyMC 5.26.1, ArviZ 0.22.0, NumPy, pandas
- **Document generation**: python-docx

## Data availability

Individual patient data can be requested through the NORspine registry application
process. See `data/README.md` for details. All analysis scripts are provided in
this repository.

## License

Analysis code is provided for transparency and reproducibility. Please cite the
associated publication if reusing any part of this work.
