# load packages
source("R/project/project_libraries.R")

# load functions
function.sources = list.files(path = "R/functions",
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
try (sapply(function.sources, source), silent = TRUE)

# load wrangled data
source("R/wrangle/wrangle_juvenile_status_data.R")

data.sources = list.files(path = "cooked_data", 
                          pattern="*ch.rds$", full.names = TRUE, 
                          ignore.case = TRUE)
sapply(data.sources, load, .GlobalEnv)

# load parameter distributions
source("R/analysis/analysis_egg_survival.R")
source("R/analysis/analysis_clutch_size.R")
source("R/analysis/analysis_hatching_sex_ratio.R")
source("R/analysis/analysis_immigration_sex_ratio.R")
source("R/analysis/analysis_mating_system.R")
source("R/analysis/analysis_fledge_age.R")
source("R/analysis/analysis_age_first_flight.R")

# consolidate
parameter_distributions <- 
  bind_rows(coucal_egg_survival,
          coucal_clutch_size,
          coucal_HSR,
          coucal_ISR,
          coucal_mating_system,
          coucal_fledge_age,
          coucal_flight_age) %>% 
  dplyr::select(trait, species, sex, mean, sd)

# filter to species of interest
BC_all_dat <- 
  status_dat_all %>% 
  filter(species == "BC")

WBC_all_dat <- 
  status_dat_all %>% 
  filter(species == "WBC")

# load WBC bootstrapped sex specific adult survival estimates
load("../coucal_ASR/coucal_ASR/output/wrangled/White-browed_Coucal_bootstrap_result_clean.rds")

WBC_adult_surival_boot_out <- 
  WBC_boot_out$survival_rates_boot %>% 
  filter(stage == "adult" & rate == "survival")

#### Black Coucals ----
niter = 1000
set.seed(14)

k_dist_BC = c(pull(filter(coucal_clutch_size, species == "BC"), mean), 
              pull(filter(coucal_clutch_size, species == "BC"), sd))

HSR_dist_BC = c(pull(filter(coucal_HSR, species == "BC"), mean), 
                pull(filter(coucal_HSR, species == "BC"), sd))

h_dist_BC = c(pull(filter(coucal_mating_system, species == "BC"), mean), 
              pull(filter(coucal_mating_system, species == "BC"), sd))

egg_surv_dist_BC = c(pull(filter(coucal_egg_survival, species == "BC"), mean),
                     pull(filter(coucal_egg_survival, species == "BC"), sd))

ISR_dist_BC = c(pull(filter(coucal_ISR, species == "BC"), mean), 
                pull(filter(coucal_ISR, species == "BC"), sd))

fledge_age_distF_BC = c(pull(filter(coucal_fledge_age, species == "BC" & sex == "F"), mean), 
                        pull(filter(coucal_fledge_age, species == "BC" & sex == "F"), sd))

fledge_age_distM_BC = c(pull(filter(coucal_fledge_age, species == "BC" & sex == "M"), mean), 
                        pull(filter(coucal_fledge_age, species == "BC" & sex == "M"), sd))

flight_age_distF_BC = c(pull(filter(coucal_flight_age, species == "BC" & sex == "F"), mean), 
                        pull(filter(coucal_flight_age, species == "BC" & sex == "F"), sd))

flight_age_distM_BC = c(pull(filter(coucal_flight_age, species == "BC" & sex == "M"), mean), 
                        pull(filter(coucal_flight_age, species == "BC" & sex == "M"), sd))

BC_hazard_ASR_bootstrap_result_w_WBC_ad_surv_stoc <-
  pbsapply(1:niter, run_bootstrap_juv_hazd_ad_surv_ASR,
    offspring = BC_all_dat,
    adult_surival_boot_out = WBC_adult_surival_boot_out,
    immigrant_pop_size = 100,
    bootstrap_name = "BC_boot_w_WBC_ad_surv_stoc",
    species = "BC",
    iter_add = 1,
    prefix_number = "boot_w_WBC_ad_surv_stoc",
    max_time = 70,
    niter = niter,
    alpha_value = 1.4,
    k_dist = k_dist_BC,
    HSR_dist = HSR_dist_BC,
    h_dist = h_dist_BC,
    egg_surv_dist = egg_surv_dist_BC,
    ISR_dist = ISR_dist_BC,
    fledge_age_distF = fledge_age_distF_BC,
    fledge_age_distM = fledge_age_distM_BC,
    flight_age_distF = flight_age_distF_BC,
    flight_age_distM = flight_age_distM_BC)

#### White-browed Coucals ----
niter = 1000
set.seed(14)

k_dist_WBC = c(pull(filter(coucal_clutch_size, species == "WBC"), mean), 
               pull(filter(coucal_clutch_size, species == "WBC"), sd))

HSR_dist_WBC = c(pull(filter(coucal_HSR, species == "WBC"), mean), 
                 pull(filter(coucal_HSR, species == "WBC"), sd))

h_dist_WBC = c(pull(filter(coucal_mating_system, species == "WBC"), mean), 
               pull(filter(coucal_mating_system, species == "WBC"), sd))

egg_surv_dist_WBC = c(pull(filter(coucal_egg_survival, species == "WBC"), mean),
                      pull(filter(coucal_egg_survival, species == "WBC"), sd))

ISR_dist_WBC = c(pull(filter(coucal_ISR, species == "WBC"), mean), 
                 pull(filter(coucal_ISR, species == "WBC"), sd))

fledge_age_distF_WBC = c(pull(filter(coucal_fledge_age, species == "WBC" & sex == "F"), mean), 
                         pull(filter(coucal_fledge_age, species == "WBC" & sex == "F"), sd))

fledge_age_distM_WBC = c(pull(filter(coucal_fledge_age, species == "WBC" & sex == "M"), mean), 
                         pull(filter(coucal_fledge_age, species == "WBC" & sex == "M"), sd))

flight_age_distF_WBC = c(pull(filter(coucal_flight_age, species == "WBC" & sex == "F"), mean), 
                         pull(filter(coucal_flight_age, species == "WBC" & sex == "F"), sd))

flight_age_distM_WBC = c(pull(filter(coucal_flight_age, species == "WBC" & sex == "M"), mean), 
                         pull(filter(coucal_flight_age, species == "WBC" & sex == "M"), sd))

WBC_hazard_ASR_bootstrap_result_w_WBC_ad_surv_stoc <-
  pbsapply(1:niter, run_bootstrap_juv_hazd_ad_surv_ASR,
           offspring = WBC_all_dat,
           adult_surival_boot_out = WBC_adult_surival_boot_out,
           immigrant_pop_size = 100,
           bootstrap_name = "WBC_boot_w_WBC_ad_surv_stoc",
           species = "WBC",
           iter_add = 1,
           prefix_number = "boot_w_WBC_ad_surv_stoc",
           max_time = 70,
           niter = niter,
           alpha_value = 1.4,
           k_dist = k_dist_WBC,
           HSR_dist = HSR_dist_WBC,
           h_dist = h_dist_WBC,
           egg_surv_dist = egg_surv_dist_WBC,
           ISR_dist = ISR_dist_WBC,
           fledge_age_distF = fledge_age_distF_WBC,
           fledge_age_distM = fledge_age_distM_WBC,
           flight_age_distF = flight_age_distF_WBC,
           flight_age_distM = flight_age_distM_WBC)

#### save and wrangle bootstrap results ----
# save model output
saveRDS(object = BC_hazard_ASR_bootstrap_result_w_WBC_ad_surv_stoc, 
        file = "output/bootstraps/hazard/cooked/BC_hazard_ASR_bootstrap_result_w_WBC_ad_surv_stoc.rds")

saveRDS(object = WBC_hazard_ASR_bootstrap_result_w_WBC_ad_surv_stoc, 
        file = "output/bootstraps/hazard/cooked/WBC_hazard_ASR_bootstrap_result_w_WBC_ad_surv_stoc.rds")

# clean up the output from the bootstrap procedure and save as rds
BC_hazard_rate_boot_tidy <- 
  hazard_boot_out_wrangle(species = "BC", niter = 1000, 
                          output_dir = "output/bootstraps/hazard/cooked/",
                          rds_file = "_hazard_ASR_bootstrap_result_w_WBC_ad_surv_stoc")

WBC_hazard_rate_boot_tidy <- 
  hazard_boot_out_wrangle(species = "WBC", niter = 1000, 
                          output_dir = "output/bootstraps/hazard/cooked/",
                          rds_file = "_hazard_ASR_bootstrap_result_w_WBC_ad_surv_stoc")

# save model output
saveRDS(object = BC_hazard_rate_boot_tidy, 
        file = "output/bootstraps/hazard/cooked/BC_haz_sur_ASR_boot_tidy_stoc.rds")

saveRDS(object = WBC_hazard_rate_boot_tidy, 
        file = "output/bootstraps/hazard/cooked/WBC_haz_sur_ASR_boot_tidy_stoc.rds")

#### Bootstrap run ----
niter = 1000
set.seed(14)

# run bootstrap procedure on Black Coucals
BC_hazard_ASR_bootstrap_result_50_ISR <-
  pbsapply(1:niter, run_bootstrap_hazard_ASR,
           offspring = BC_all_dat,
           k = 4,
           HSR = 0.4955,
           h = 1/2.9,
           egg_survival = 0.32,
           adult_survival_rate = 0.3,
           ISR = 0.5,
           immigrant_pop_size = 100,
           fledge_age = 15,
           flight_age = 36,
           bootstrap_name = "BC_boot_test",
           species = "BC",
           iter_add = 1,
           prefix_number = "boot_test_",
           max_time = 70)

# save model output
saveRDS(object = BC_hazard_ASR_bootstrap_result_50_ISR, 
        file = "output/bootstraps/hazard/cooked/BC_hazard_ASR_bootstrap_result_50_ISR.rds")