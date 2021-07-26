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
                          rds_file = "_hazard_ASR_bootstrap_result_w_WBC_ad_surv_stoc")

# clean up the output from the bootstrap procedure and save as rds
WBC_hazard_rate_boot_tidy <- 
  hazard_boot_out_wrangle(species = "WBC", niter = 1000, 
                          output_dir = "output/bootstraps/hazard/cooked/",
                          rds_file = "_hazard_ASR_bootstrap_result_w_WBC_ad_surv_stoc")

#### Run BC JSR Bootstrap ----
BC_F_rates <- 
  BC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  filter(sex == "Female" & species == "BC")

BC_M_rates <- 
  BC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  filter(sex == "Male" & species == "BC")

# niter = 1000
# F_rates = BC_F_rates
# M_rates = BC_M_rates
# immigrant_pop_size = 100
# ISR_rates = filter(BC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "ISR")
# h_rates = filter(BC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "h")
# k_rates = filter(BC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "k")
# species = "BC"
# HSR_rates = filter(BC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "HSR")
# egg_survival_rates = filter(BC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "egg")

BC_JSR_run_stoch <- 
  JSR_bootstrap_stoch(niter = 1000, 
                      F_rates = BC_F_rates,
                      M_rates = BC_M_rates, 
                      immigrant_pop_size = 100, 
                      ISR_rates = filter(BC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "ISR"), 
                      h_rates = filter(BC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "h"), 
                      k_rates = filter(BC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "k"), 
                      species = "BC",
                      HSR_rates = filter(BC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "HSR"), 
                      egg_survival_rates = filter(BC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "egg"))

HSR_rates_5050 <- 
  data.frame(sex = NA,
             value = 0.5,
             stage = "HSR",
             rate = "fecundity",
             iter = c(1:1000),
             species = "BC")

BC_JSR_run_5050_HSR <- 
  JSR_bootstrap_stoch(niter = 1000, 
                      F_rates = BC_F_rates,
                      M_rates = BC_M_rates, 
                      immigrant_pop_size = 100, 
                      ISR_rates = filter(BC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "ISR"), 
                      h_rates = filter(BC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "h"), 
                      k_rates = filter(BC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "k"), 
                      species = "BC",
                      HSR_rates = HSR_rates_5050, 
                      egg_survival_rates = filter(BC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "egg"))

# save model output
saveRDS(object = BC_JSR_run, 
        file = "output/bootstraps/hazard/cooked/BC_hazard_JSR_bootstrap_result_stoch.rds")

#### Run WBC JSR Bootstrap ----
WBC_F_rates <- 
  WBC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  filter(sex == "Female" & species == "WBC")

WBC_M_rates <- 
  WBC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  filter(sex == "Male" & species == "WBC")

# WBC_h = 1/1.1
# WBC_HSR = 0.5198
# WBC_k = 4
# WBC_ISR = 0.524
# WBC_egg_survival = 0.18

WBC_JSR_run <- 
  JSR_bootstrap_stoch(niter = 1000, 
                      F_rates = WBC_F_rates,
                      M_rates = WBC_M_rates, 
                      immigrant_pop_size = 100, 
                      ISR_rates = filter(WBC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "ISR"), 
                      h_rates = filter(WBC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "h"), 
                      k_rates = filter(WBC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "k"), 
                      species = "WBC",
                      HSR_rates = filter(WBC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "HSR"), 
                      egg_survival_rates = filter(WBC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "egg"))

HSR_rates_5050 <- 
  data.frame(sex = NA,
             value = 0.5,
             stage = "HSR",
             rate = "fecundity",
             iter = c(1:1000),
             species = "WBC")

WBC_JSR_run_5050_HSR <- 
  JSR_bootstrap_stoch(niter = 1000, 
                      F_rates = WBC_F_rates,
                      M_rates = WBC_M_rates, 
                      immigrant_pop_size = 100, 
                      ISR_rates = filter(WBC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "ISR"), 
                      h_rates = filter(WBC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "h"), 
                      k_rates = filter(WBC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "k"), 
                      species = "WBC",
                      HSR_rates = HSR_rates_5050, 
                      egg_survival_rates = filter(WBC_hazard_rate_boot_tidy$vital_rate_ests_boot, stage == "egg"))
# save model output
saveRDS(object = WBC_JSR_run, 
        file = "output/bootstraps/hazard/cooked/WBC_hazard_JSR_bootstrap_result_stoch.rds")
