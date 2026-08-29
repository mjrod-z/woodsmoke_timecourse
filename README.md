# Small & Large Airway Cytokine Timecourse Analysis

Analysis code for studying cytokine inflammatory responses in small airway epithelial (SAE) and large airway epithelial (LAE) cells exposed to biomass smoke particulates across two timepoints (4 h and 144 h post-exposure).

## Overview

This repository contains modular R code for:

- MSD cytokine data loading and preprocessing
- LLOD-based imputation
- Timepoint-aware (`_4`, `_144`) WSTC assembly and mixed-effects model screening across airway type, hormone treatment, sex, and timepoint
- DiD-annotated timecourse spaghetti plots and DiD heatmaps relative to PBS baseline
- Cytokine dotplots, barplots, and histograms (faceted by `CELLTYPE × TIMEPOINT`)
- Variance partitioning and PCA across timepoints
- LDH / TEER analysis scaffolding

## DiD timecourse interpretation

Timecourse inference now includes a difference-in-differences (DiD) mixed-model workflow for cytokines:
\[
[(\text{Exposure} - \text{PBS\_Control})_{144h}] - [(\text{Exposure} - \text{PBS\_Control})_{4h}]
\]
using `log2(Value + pseudocount)` and donor random intercepts.
Calls use BH-FDR (`q < 0.05`) plus direction:
- **Increase over time vs PBS**: `q < 0.05` and estimate > 0
- **Decrease over time vs PBS**: `q < 0.05` and estimate < 0
- **No significant change over time vs PBS**: otherwise

Timepoint visualizations use this DiD table for cytokine ranking, annotations, and heatmap values, so plots do not present within-exposure 4h-vs-144h statistics as the main inference.

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
