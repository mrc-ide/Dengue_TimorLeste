# 1. Seroprevalence Analysis in Timor-Leste

This repository contains R scripts to process and analyse a dengue IgG seroprevalence dataset from Timor-Leste using a Bayesian model developed in **RStan**.

## Workflow

- **Dataset:** `Dataset_TL` with variables: EA (village number), IgG dengue result (positive/negative), age, household number, number of people tested from the same household.  
- **Step 1:** Run `1_prod_database.R` to load `Dataset_TL`, clean and tidy the data.  
- **Step 2:** Run `2_run_model.R` to use `model_TL.stan` with RStan for FOI estimation. 

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
install.packages(c("tidyverse", "rstan", "ggplot2", "ggpubr", "binom", "Hmisc", "gridExtra"))

# 2. Estimating force-of-infection across Timor-Leste

The 'spatial_model' folder contains R scripts to process the FOI estimates from the serocatalytic model (in step 1, and see 'FOI_model' folder) and to estimate FOI for the whole of Timor-Leste using a spatial model implemented in **R-INLA**.

## Workflow

- Run the spatial_model_script.R which reads in 'FOI_catalytic_model_posteriors.xlsx' (pre-run estimates from the sero-catalytic model in step 1 of this repository) and runs the spatial model on these sampled FOI estimates.


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

