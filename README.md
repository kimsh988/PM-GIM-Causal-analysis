# Long-term Particulate Matter Exposure and Gastric Intestinal Metaplasia

This repository contains the R code used for the analyses in the study examining the causal effect of long-term particulate matter exposure on gastric intestinal metaplasia.

## Study

The study evaluates hypothetical interventions on annual average PM2.5 and PM10 concentrations using the parametric g-formula.

## Code

The repository contains two R scripts:

### `PM-GIM_gformula.R`

This is the main analysis script and contains the code for:

* Data preprocessing
* Parametric g-formula analyses
* Sensitivity analyses
* Subgroup analyses
* Generation of figures and plots

### `gform_plot_functions.R`

This script contains custom functions used to generate plots from the parametric g-formula results.

The plotting functions in `gform_plot_functions.R` are sourced and used within `PM-GIM_gformula.R`.

## Data Availability

Individual-level data are not publicly available due to institutional and privacy restrictions.

## Software

Analyses were conducted in R.
