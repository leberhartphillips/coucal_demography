# load packages
source("R/project/project_libraries.R")

# load functions
function.sources = list.files(path = "R/functions",
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
sapply(function.sources, source)

# load parameter distributions
source("R/analysis/analysis_egg_survival.R")
source("R/analysis/analysis_clutch_size.R")
source("R/analysis/analysis_hatching_sex_ratio.R")
source("R/analysis/analysis_immigration_sex_ratio.R")
source("R/analysis/analysis_mating_system.R")
source("R/analysis/analysis_age_first_flight.R")
source("R/analysis/analysis_fledge_age.R")

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

# clean up the output from the bootstrap procedure
BC_hazard_rate_boot_tidy <- 
  hazard_boot_out_wrangle(species = "BC", niter = 1000, 
                          output_dir = "output/bootstraps/hazard/cooked/",
                          rds_file = "_hazard_ASR_bootstrap_result_w_WBC_ad_surv_stoc")

WBC_hazard_rate_boot_tidy <- 
  hazard_boot_out_wrangle(species = "WBC", niter = 1000, 
                          output_dir = "output/bootstraps/hazard/cooked/",
                          rds_file = "_hazard_ASR_bootstrap_result_w_WBC_ad_surv_stoc")

hazard_rate_boot_tidy <- 
  bind_rows(BC_hazard_rate_boot_tidy, WBC_hazard_rate_boot_tidy)

# extract age-specific survival rates from the bootstrap procedure and summarize by age
survival_rates_boot_summary <- 
  bind_rows(BC_hazard_rate_boot_tidy$hazard_rates_boot,
            WBC_hazard_rate_boot_tidy$hazard_rates_boot) %>% 
  Rmisc::summarySE(.,
                   measurevar = "fit",
                   groupvars = c("species", "sex", "age"),
                   conf.interval = 0.95) %>% 
  mutate(fit = 1 - fit) %>% 
  dplyr::select(species, sex, age, fit)

#### Arguments ----
species_name = "BC"
imm_N_base = 100
ISR_rate = pull(filter(parameter_distributions, 
                       species == species_name & trait == "immigration_sex_ratio"), 
                mean)
k_rate = pull(filter(parameter_distributions, 
                     species == species_name & trait == "clutch_size"), 
              mean)
h_rate = (1 / pull(filter(parameter_distributions, 
                          species == species_name & trait == "mating_system"), 
                   mean))
HSR_rate = pull(filter(parameter_distributions, 
                       species == species_name & trait == "hatching_sex_ratio"), 
                mean)
egg_S_rate = pull(filter(parameter_distributions, 
                         species == species_name & trait == "egg_survival"), 
                  mean)
JSR_base = 
  bind_rows(BC_JSR_out,
            WBC_JSR_out) %>%
  mutate(species = factor(species, levels = c("BC", "WBC"))) %>% 
  dplyr::group_by(species, age) %>%
  dplyr::summarise(ucl_JSR = stats::quantile(JSR, (1 - CI)/2, na.rm = TRUE),
                   lcl_JSR = stats::quantile(JSR, 1 - (1 - CI)/2, na.rm = TRUE),
                   avg_JSR = mean(JSR, na.rm = TRUE),
                   med_JSR = median(JSR, na.rm = TRUE),
                   max_JSR = max(JSR, na.rm = TRUE),
                   min_JSR = min(JSR, na.rm = TRUE)) %>% 
  filter(age == (length(unique(survival_rates$age)) - 2) & species == species_name) %>% 
  pull(avg_JSR)

#### analysis ----
BC_VR_treat <-
  make_mprime_and_treat_summary_JSR(survival_rates_boot_summary = survival_rates_boot_summary,
                                    h_rate = h_rate,
                                    HSR_rate = HSR_rate,
                                    k_rate = k_rate,
                                    ISR_rate = ISR_rate, 
                                    imm_N_base = imm_N_base,
                                    egg_S_rate = egg_S_rate,
                                    species_name = "BC", 
                                    mprime = FALSE)

BC_VR_mprime_male <- 
  make_mprime_and_treat_summary_JSR(survival_rates_boot_summary = survival_rates_boot_summary,
                                    h_rate = h_rate,
                                    HSR_rate = HSR_rate,
                                    k_rate = k_rate,
                                    ISR_rate = ISR_rate, 
                                    imm_N_base = imm_N_base,
                                    egg_S_rate = egg_S_rate,
                                    species_name = "BC", 
                                    base_sex = "Male")

BC_VR_mprime_female <- 
  make_mprime_and_treat_summary_JSR(survival_rates_boot_summary = survival_rates_boot_summary,
                                    h_rate = h_rate,
                                    HSR_rate = HSR_rate,
                                    k_rate = k_rate,
                                    ISR_rate = ISR_rate, 
                                    imm_N_base = imm_N_base,
                                    egg_S_rate = egg_S_rate,
                                    species_name = "BC", 
                                    base_sex = "Female")

BC_treatment_JSR_analysis <- 
  JSR_bootstrap_deter(VR_dataframe = BC_VR_treat,
                      immigrant_pop_size = 100,
                      species = "BC")

BC_M_prime_JSR_analysis_male <- 
  JSR_bootstrap_deter(VR_dataframe = BC_VR_mprime_male,
                      immigrant_pop_size = 100,
                      species = "BC")

BC_M_prime_JSR_analysis_female <- 
  JSR_bootstrap_deter(VR_dataframe = BC_VR_mprime_female,
                      immigrant_pop_size = 100,
                      species = "BC")

BC_JSR_treat <- 
  BC_treatment_JSR_analysis[nrow(BC_treatment_JSR_analysis), "JSR"] 
BC_JSR_treat

BC_JSR_mprime_male <- 
  BC_M_prime_JSR_analysis_male[nrow(BC_M_prime_JSR_analysis_male), "JSR"] 
BC_JSR_mprime_male

BC_JSR_mprime_female <- 
  BC_M_prime_JSR_analysis_female[nrow(BC_M_prime_JSR_analysis_female), "JSR"] 
BC_JSR_mprime_female 

BC_treat_sensitivity_analysis <- 
  sensitivity_analysis_JSR(VR_dataframe = BC_VR_treat, 
                           JSR_base = BC_JSR_treat)

saveRDS(object = BC_treat_sensitivity_analysis, 
        file = "output/sensitivity_analysis/BC_treat_sensitivity_analysis_JSR.rds")

BC_Mprime_sensitivity_analysis_male <- 
  sensitivity_analysis_JSR(VR_dataframe = BC_VR_mprime_male, 
                           JSR_base = BC_JSR_mprime_male)

saveRDS(object = BC_Mprime_sensitivity_analysis_male, 
        file = "output/sensitivity_analysis/BC_Mprime_sensitivity_analysis_male_JSR.rds")

BC_Mprime_sensitivity_analysis_female <- 
  sensitivity_analysis_JSR(VR_dataframe = BC_VR_mprime_female, 
                           JSR_base = BC_JSR_mprime_female)

saveRDS(object = BC_Mprime_sensitivity_analysis_female, 
        file = "output/sensitivity_analysis/BC_Mprime_sensitivity_analysis_female_JSR.rds")

BC_LTRE_male <- 
  LTRE_analysis_JSR(Mprime_sens = BC_Mprime_sensitivity_analysis_male, 
                    VR_dataframe = BC_VR_treat,
                    base_sex = "male")

BC_LTRE_female <- 
  LTRE_analysis_JSR(Mprime_sens = BC_Mprime_sensitivity_analysis_female, 
                    VR_dataframe = BC_VR_treat,
                    base_sex = "female")

# WBC_VR_treat <- 
#   make_mprime_matrix_JSR(survival_rates_boot_summary = survival_rates_boot_summary,
#                          h_rate = h_rate,
#                          HSR_rate = HSR_rate,
#                          k_rate = k_rate,
#                          ISR_rate = ISR_rate, 
#                          imm_N_base = imm_N_base,
#                          egg_S_rate = egg_S_rate,
#                          species_name = "WBC", 
#                          mprime = FALSE)
# 
# WBC_VR_mprime_male <- 
#   make_mprime_matrix_JSR(survival_rates_boot_summary = survival_rates_boot_summary,
#                          h_rate = h_rate,
#                          HSR_rate = HSR_rate,
#                          k_rate = k_rate,
#                          ISR_rate = ISR_rate, 
#                          imm_N_base = imm_N_base,
#                          egg_S_rate = egg_S_rate,
#                          species_name = "WBC", 
#                          base_sex = "Male")
# 
# WBC_VR_mprime_female <- 
#   make_mprime_matrix_JSR(survival_rates_boot_summary = survival_rates_boot_summary,
#                          h_rate = h_rate,
#                          HSR_rate = HSR_rate,
#                          k_rate = k_rate,
#                          ISR_rate = ISR_rate, 
#                          imm_N_base = imm_N_base,
#                          egg_S_rate = egg_S_rate,
#                          species_name = "WBC", 
#                          base_sex = "Female")



# WBC_treatment_JSR_analysis <- 
#   JSR_bootstrap_deter(VR_dataframe = WBC_VR_treat,
#                       immigrant_pop_size = 100,
#                       species = "WBC")
# 
# WBC_ASR_treat <- 
#   WBC_treatment_JSR_analysis[nrow(WBC_treatment_JSR_analysis), "JSR"] 
# 
# WBC_ASR_treat
# 
# WBC_M_prime_ASR_analysis_male <- 
#   matrix_ASR(M = WBC_M_prime_matrix_male, 
#              h = WBC_VR_mprime_male$h, 
#              k = WBC_VR_mprime_male$k,
#              HSR = WBC_VR_mprime_male$HSR,
#              ISR = WBC_VR_mprime_male$ISR,
#              immigrant_pop_size = 100,
#              iterations = 100,
#              num_boot = 1,
#              iter_add = 1,
#              species = "BC")
# 
# WBC_ASR_mprime_male <- WBC_M_prime_ASR_analysis_male$ASR
# WBC_ASR_mprime_male
# 
# WBC_M_prime_ASR_analysis_female <- 
#   matrix_ASR(M = WBC_M_prime_matrix_female, 
#              h = WBC_VR_mprime_female$h, 
#              k = WBC_VR_mprime_female$k,
#              HSR = WBC_VR_mprime_female$HSR,
#              ISR = WBC_VR_mprime_female$ISR,
#              immigrant_pop_size = 100,
#              iterations = 100,
#              num_boot = 1,
#              iter_add = 1,
#              species = "BC")
# 
# WBC_ASR_mprime_female <- WBC_M_prime_ASR_analysis_female$ASR
# WBC_ASR_mprime_female
# 
# BC_lambda_treat <- BC_treatment_ASR_analysis$lambda
# BC_lambda_treat
# 
# BC_lambda_mprime_male <- BC_M_prime_ASR_analysis_male$lambda
# BC_lambda_mprime_male
# 
# BC_lambda_mprime_female <- BC_M_prime_ASR_analysis_female$lambda
# BC_lambda_mprime_female
# 
# WBC_lambda_treat <- WBC_treatment_ASR_analysis$lambda
# WBC_lambda_treat
# 
# WBC_lambda_mprime_male <- WBC_M_prime_ASR_analysis_male$lambda
# WBC_lambda_mprime_male
# 
# WBC_lambda_mprime_female <- WBC_M_prime_ASR_analysis_female$lambda
# WBC_lambda_mprime_female


# WBC_treat_sensitivity_analysis <- 
#   sensitivity_analysis(vital_rate_summary = WBC_VR_treat, 
#                        matrix_str = matrix_structure, 
#                        h = WBC_VR_treat$h, 
#                        k = WBC_VR_treat$k, 
#                        HSR = WBC_VR_treat$HSR, 
#                        niter = 1000, 
#                        ASR = WBC_ASR_treat,
#                        lambda = WBC_lambda_treat, 
#                        ISR = WBC_VR_treat$ISR, immigrant_pop_size = 100)
# 
# WBC_Mprime_sensitivity_analysis_male <- 
#   sensitivity_analysis(vital_rate_summary = WBC_VR_mprime_male, 
#                        matrix_str = matrix_structure, 
#                        h = WBC_VR_mprime_male$h, 
#                        k = WBC_VR_mprime_male$k, 
#                        HSR = WBC_VR_mprime_male$HSR, 
#                        niter = 1000, 
#                        ASR = WBC_ASR_mprime_male,
#                        lambda = WBC_lambda_mprime_male, 
#                        ISR = WBC_VR_mprime_male$ISR, immigrant_pop_size = 100)
# 
# WBC_Mprime_sensitivity_analysis_female <- 
#   sensitivity_analysis(vital_rate_summary = WBC_VR_mprime_female, 
#                        matrix_str = matrix_structure, 
#                        h = WBC_VR_mprime_female$h, 
#                        k = WBC_VR_mprime_female$k, 
#                        HSR = WBC_VR_mprime_female$HSR, 
#                        niter = 1000, 
#                        ASR = WBC_ASR_mprime_female,
#                        lambda = WBC_lambda_mprime_female, 
#                        ISR = WBC_VR_mprime_female$ISR, immigrant_pop_size = 100)
# ```
#conduct the LTRE comparing the two matrices
# ```{r}

WBC_LTRE_male <- 
  LTRE_analysis(Mprime_sens = WBC_Mprime_sensitivity_analysis_male, 
                matrix_str = matrix_str, 
                vital_rates = WBC_VR_treat,
                species_name = "WBC",
                sex = "male")

WBC_LTRE_female <- 
  LTRE_analysis(Mprime_sens = WBC_Mprime_sensitivity_analysis_female, 
                matrix_str = matrix_str, 
                vital_rates = WBC_VR_treat,
                species_name = "WBC",
                sex = "female")

LTRE_coucal_male_ASR <- 
  bind_rows(BC_LTRE_male$LTRE_ASR, WBC_LTRE_male$LTRE_ASR) %>% 
  mutate(sex = "male")

LTRE_coucal_female_ASR <- 
  bind_rows(BC_LTRE_female$LTRE_ASR, WBC_LTRE_female$LTRE_ASR) %>% 
  mutate(sex = "female")

LTRE_coucal_ASR <- rbind(LTRE_coucal_male_ASR, LTRE_coucal_female_ASR)

LTRE_coucal_ASR$parameter <- 
  factor(LTRE_coucal_ASR$parameter, 
         levels = c("Hatching sex ratio",
                    "Nestling survival",
                    "Groundling survival",
                    "Fledgling survival",
                    "Adult survival",
                    "Mating system",
                    "Immigrant sex ratio"))