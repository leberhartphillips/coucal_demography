# load packages
source("R/project/project_libraries.R")

# load functions
function.sources = list.files(path = "R/functions",
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
sapply(function.sources, source)

# clean up the output from the bootstrap procedure and save as rds
BC_hazard_rate_boot_tidy <- 
  hazard_boot_out_wrangle(species = "BC", niter = 1000, 
                          output_dir = "output/bootstraps/hazard/cooked/",
                          rds_file = "_hazard_ASR_bootstrap_result_one")

# clean up the output from the bootstrap procedure and save as rds
WBC_hazard_rate_boot_tidy <- 
  hazard_boot_out_wrangle(species = "WBC", niter = 1000, 
                          output_dir = "output/bootstraps/hazard/cooked/",
                          rds_file = "_hazard_ASR_bootstrap_result_one")

#### Run BC JSR Bootstrap ----
BC_F_rates <- 
  BC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  filter(sex == "Female" & species == "BC")

BC_M_rates <- 
  BC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  filter(sex == "Male" & species == "BC")

BC_h = 1/2.9
BC_HSR = 0.4955
BC_k = 4
BC_ISR = 0.738
BC_egg_survival = 0.32

BC_JSR_run <- 
  JSR_bootstrap(niter = 1000, 
                F_rates = BC_F_rates,
                M_rates = BC_M_rates, 
                immigrant_pop_size = 100, 
                ISR = BC_ISR, h = BC_h, k = BC_k, species = "BC",
                HSR = BC_HSR, egg_survival = BC_egg_survival)

# save model output
saveRDS(object = BC_JSR_run, 
        file = "output/bootstraps/hazard/cooked/BC_hazard_JSR_bootstrap_result.rds")

#### Run WBC JSR Bootstrap ----
WBC_F_rates <- 
  WBC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  filter(sex == "Female" & species == "WBC")

WBC_M_rates <- 
  WBC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  filter(sex == "Male" & species == "WBC")

WBC_h = 1/1.1
WBC_HSR = 0.5198
WBC_k = 4
WBC_ISR = 0.524
WBC_egg_survival = 0.18

WBC_JSR_run <- 
  JSR_bootstrap(niter = 1000, 
                F_rates = WBC_F_rates,
                M_rates = WBC_M_rates, 
                immigrant_pop_size = 100, 
                ISR = WBC_ISR, h = WBC_h, k = WBC_k, species = "WBC",
                HSR = WBC_HSR, egg_survival = WBC_egg_survival)

# save model output
saveRDS(object = WBC_JSR_run, 
        file = "output/bootstraps/hazard/cooked/WBC_hazard_JSR_bootstrap_result.rds")