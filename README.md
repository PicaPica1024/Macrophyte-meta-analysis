# Reproducibility files

This folder contains the data and R code accompanying the manuscript **“Context dependence of interactions among global change stressors on submerged plants: roles of null models, experimental conditions, and plant traits.”**

## Files

- `Analysis code.R`: complete data cleaning, effect-size calculation, meta-analysis, moderator analysis, model selection, publication-bias analysis, sensitivity analysis, and figure-generation workflow.
- `data.xlsx`: observation-level meta-analysis dataset read by the R script.
- `plant traits.xlsx`: species-level trait classifications and supporting references.

## Running the analysis

1. Place `Analysis code.R`, `data.xlsx`, and `plant traits.xlsx` in the same directory.
2. Set the R working directory to that directory.
3. Run `Analysis code.R` from top to bottom.

The script creates an `analysis outputs` subdirectory and writes the generated figures and tables there. It reproduces the data-derived manuscript Figures 2–5; Figure 2 is also saved as separate chord-diagram and Venn-diagram panels for inspection. The conceptual Figure 1 is not generated from the analysis data. The analyses were developed in R 4.4.1.

Required packages are `meta`, `metafor`, `ggplot2`, `tidyr`, `readxl`, `clubSandwich`, `dmetar`, `MuMIn`, `dplyr`, `cowplot`, `patchwork`, `grid`, `circlize`, `gridExtra`, `eulerr`, and `png`.

## Key analytical settings

- Additive, multiplicative, and dominance null-model predictions are compared using the common Hedges’ *d* standardization described in the manuscript.
- Non-estimable effect sizes, non-positive sampling variances, and observations with `|d| > 30` are excluded.
- Multilevel models include study and observation nested within study as random effects and use study-clustered CR2 inference.
- The primary within-study sampling correlation is 0.6; sensitivity analyses repeat the overall models using correlations of 0, 0.3, 0.6, and 0.9 and use leave-one-study-out analyses.
- Model selection uses maximum likelihood, AICc, and retains models with `Delta AICc < 3`.
