# Small & Large Airway Cytokine Timecourse Analysis

Analysis code for studying cytokine inflammatory responses in small airway epithelial (SAE) and large airway epithelial (LAE) cells exposed to biomass smoke particulates across two timepoints (4 h and 144 h post-exposure).

## Overview

This repository contains modular R code for:

- MSD cytokine data loading and preprocessing
- LLOD-based imputation
- Timepoint-aware (`_4`, `_144`) WSTC assembly and mixed-effects model screening across airway type, hormone treatment, sex, and timepoint
- Cytokine dotplots, barplots, and histograms (faceted by `CELLTYPE × TIMEPOINT`)
- Variance partitioning and PCA across timepoints
- LDH / TEER analysis scaffolding

## Project structure

```text
.
├── .gitignore
├── 01_timecourse_analysis.Rmd
├── README.md
├── _load_all.R
├── config.R
├── functions_data.R
├── functions_analysis.R
├── functions_plots.R
├── R/
├── analysis/
├── data/
│   ├── raw/
│   └── derived/
└── output/
    ├── figures/
    ├── tables/
    └── reports/
```
