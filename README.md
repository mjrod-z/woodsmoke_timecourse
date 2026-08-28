# Small & Large Airway Cytokine Analysis

Analysis code for studying cytokine inflammatory responses in small airway epithelial (SAE) and large airway epithelial (LAE) cells exposed to biomass smoke particulates.

## Overview

This repository contains modular R code for:

- MSD cytokine data loading and preprocessing
- LLOD-based imputation
- Mixed-effects model screening across airway type, hormone treatment, and sex
- Cytokine dotplots, barplots, and histograms
- Smoker variance partition analysis
- LDH / TEER analysis scaffolding

## Project structure

```text
.
├── 01_sala_analysis.Rmd
├── README.md
├── _load_all.R
├── config.R
├── functions_data.R
├── functions_analysis.R
├── functions_plots.R
├── data/
│   ├── raw/
│   └── derived/
└── output/
    ├── figures/
    ├── tables/
    └── reports/