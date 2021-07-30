# Reproducible datasets and code for:
## Early life stepping stones to sex role evolution: Sex-specific offspring survival dynamics of two sympatric coucal species differing in mating system
#### Ignas Safari, Luke J. Eberhart-Hertel, and Wolfgang Goymann
#### *In prep* 

In this repository you can find all the raw data and code needed to reproduce our investigation of sex-specific survival of offspring black coucals (_Centropus grillii_) and white-browed coucals (_Centropus superciliosus_) monitored over a 15-year period in south-western Tanzania.

#### Repository Contents
[**`R/`**](https://github.com/leberhartphillips/coucal_demography/tree/main/R)
  - [`project/`](https://github.com/leberhartphillips/coucal_demography/blob/main/R/project) folder contains scripts for housekeeping code (loading libraries and custom plotting themes)
  - [`wrangle/`](https://github.com/leberhartphillips/coucal_demography/blob/main/R/wrangle) folder contains code used to wrangle our raw field data into the format used in the models presented in this investigation
  - [`function/`](https://github.com/leberhartphillips/coucal_demography/blob/main/R/function) folder contains custom functions used in our analysis (i.e., referenced in the `analysis_X.R` scripts)
  - [`analysis/`](https://github.com/leberhartphillips/coucal_demography/blob/main/R/analysis) folder contains code used for modeling and exporting output
  - [`results/`](https://github.com/leberhartphillips/coucal_demography/blob/main/R/results) folder contains code used to process the model output into figures, tables, and inferential statistics
  
[**`data/`**](https://github.com/leberhartphillips/coucal_demography/tree/main/data)
  - `status_dat_all.csv` are the survival data from radio tracked offspring of both species (used in the [`analysis_hazard_boot.R`](https://github.com/leberhartphillips/coucal_demography/blob/main/R/analysis/analysis_hazard_boot.R) script)
  - `egg_data.csv` are the egg data of both species (used in the [`analysis_egg_survival.R`](https://github.com/leberhartphillips/coucal_demography/blob/main/R/analysis/analysis_egg_survival.R),  [`analysis_clutch_size.R`](https://github.com/leberhartphillips/coucal_demography/blob/main/R/analysis/analysis_clutch_size.R), and [`analysis_hatching_sex_ratio.R`](https://github.com/leberhartphillips/coucal_demography/blob/main/R/analysis/analysis_hatching_sex_ratio.R) scripts)
  - `mating_dat.csv` are the egg survival data of both species (used in the [`analysis_immigration_sex_ratio.R`](https://github.com/leberhartphillips/coucal_demography/blob/main/R/analysis/analysis_immigration_sex_ratio.R) and [`analysis_mating_system.R`](https://github.com/leberhartphillips/coucal_demography/blob/main/R/analysis/analysis_mating_system.R) scripts)
  - `fledge_dat.csv` are the fledging data of both species (used in the [`analysis_fledge_age.R`](https://github.com/leberhartphillips/coucal_demography/blob/main/R/analysis/analysis_fledge_age.R) script)
  - `flight_dat.csv` are the flight data of both species (used in the [`analysis_age_first_flight.R`](https://github.com/leberhartphillips/coucal_demography/blob/main/R/analysis/analysis_age_first_flight.R) script)
  
[**`output/`**](https://github.com/leberhartphillips/coucal_demography/tree/main/output)
  - `offpring_survival_bootstrap.rds` 1000 smoothed splines of the hazard function of offspring survival since hatch
  - `BC_treat_sensitivity_analysis_JSR.rds` sensitivities of treatment life table for black coucal
  - `BC_Mprime_sensitivity_analysis_male_JSR.rds` sensitivities of male-based M-prime life table for black coucal
  - `BC_Mprime_sensitivity_analysis_female_JSR.rds` sensitivities of female-based M-prime life table for black coucal
  - `WBC_treat_sensitivity_analysis_JSR.rds` sensitivities of treatment life table for white-browed coucal
  - `WBC_Mprime_sensitivity_analysis_male_JSR.rds` sensitivities of male-based M-prime life table for white-browed coucal
  - `WBC_Mprime_sensitivity_analysis_female_JSR.rds` sensitivities of female-based M-prime life table for white-browed coucal
  - `LTRE_coucal_JSR.rds` life table response experiment contributions for both species and all demographic parameters
  
[**`products/`**](https://github.com/leberhartphillips/coucal_demography/tree/main/products)