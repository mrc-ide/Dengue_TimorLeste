# 1. Seroprevalence Analysis in Timor-Leste

This repository contains R scripts to process and analyse a dengue IgG seroprevalence dataset from Timor-Leste using a Bayesian hierarchical model developed in **RStan**.

## Workflow

- **Dataset:** `Dataset_TL` with variables: EA (village number), IgG dengue result (positive/negative), age, household number, number of people tested from the same household.  
- **Model:** Run `1_run_hierarchical_model` to load `Dataset_TL`, clean and tidy the data and it will use the model described in `Hierarchical_model.stan` file to estimate the EA-level FOI. 

## Package versions

The analysis was run using the following R package versions:

| Package    | Version |
|------------|---------|
| R          | 4.3.2   |
| tidyverse  | 2.0.0   |
| rstan      | 2.32.6  |
| ggplot2    | 3.5.1   |
| ggpubr     | 0.6.0   |
| binom      | 1.1.1.1 |
| Hmisc      | 5.2.2   |
| gridExtra  | 2.3     |

To download these packages please follow these steps:
install.packages(c("tidyverse", "rstan", "ggplot2", "ggpubr", "binom", "Hmisc", "gridExtra", "purr", "truncnorm"))

# 2. Estimating force-of-infection across Timor-Leste

## Workflow

- **Dataset:** 'spatial_model_data/FOI_catalytic_model_posteriors.xlsx' (pre-run estimates from the sero-catalytic model of this repository)  
- **Optional Step 1:** Run the 3_setup_spatial_model.R script which creates the mesh and spde needed in step 4.
- **Step 2:** Run the 4_run_spatial_model.R script which runs the spatial model on sampled FOI estimates from the serocatalytic model to estimate FOI for the whole of Timor-Leste using a spatial model implemented in **R-INLA**.


## Package versions

The analysis was run using the following R package versions:

| Package    | Version  |
|------------|----------|
| INLA       | 24.06.27 |
| INLAutils  | 0.0.6    |
| sf         | 1.0-19   |

To download R-INLA and INLAutils please follow these steps: 
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
install.packages('devtools')
devtools::install_github('timcdlucas/INLAutils')

