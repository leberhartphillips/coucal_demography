# load packages
source("R/project/project_libraries.R")

# load functions
function.sources = list.files(path = "R/functions",
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
sapply(function.sources, source)


BC_h = 1/2.9
BC_HSR = 0.4955
BC_k = 4
BC_ISR = 0.738

WBC_h = 1/1.1
WBC_HSR = 0.5198
WBC_k = 4
WBC_ISR = 0.524

niter = 1000

# load output
BC_hazard_rate_boot <- 
  readRDS("output/bootstraps/hazard/cooked/BC_hazard_ASR_bootstrap_result_one.rds")

# load output
WBC_hazard_rate_boot <- 
  readRDS("output/bootstraps/hazard/cooked/WBC_hazard_ASR_bootstrap_result_one.rds")

# load output
BC_hazard_rate_boot_tidy <- 
  readRDS("output/bootstraps/hazard/cooked/BC_haz_sur_ASR_boot_tidy.rds")

# load output
WBC_hazard_rate_boot_tidy <- 
  readRDS("output/bootstraps/hazard/cooked/WBC_haz_sur_ASR_boot_tidy.rds")

BC_hazard_rate_boot_tidy$vital_rate_ests_boot$iter <- as.factor(BC_hazard_rate_boot_tidy$vital_rate_ests_boot$iter)
WBC_hazard_rate_boot_tidy$vital_rate_ests_boot$iter <- as.factor(WBC_hazard_rate_boot_tidy$vital_rate_ests_boot$iter)

# summarize the vital rates
survival_rates_boot_summary <-
  bind_rows(BC_hazard_rate_boot_tidy$vital_rate_ests_boot, WBC_hazard_rate_boot_tidy$vital_rate_ests_boot) %>% 
  mutate(vital_rate = paste(sex, stage, rate, sep = "_")) %>% 
  Rmisc::summarySE(.,
                   measurevar = "value",
                   groupvars = c("vital_rate", "species"),
                   conf.interval = 0.95) %>% 
  arrange(species, vital_rate)


BC_VR_treat <- 
  make_treat_matrix(survival_rates_boot_summary = survival_rates_boot_summary,
                    h = BC_h,
                    HSR = BC_HSR,
                    k = BC_k,
                    ISR = BC_ISR,
                    species = "BC")

BC_VR_mprime_male <- 
  make_mprime_matrix(survival_rates_boot_summary = survival_rates_boot_summary,
                     h = BC_h,
                     HSR = BC_HSR,
                     k = BC_k,
                     ISR = BC_ISR,
                     species = "BC",
                     sex = "male")

BC_VR_mprime_female <- 
  make_mprime_matrix(survival_rates_boot_summary = survival_rates_boot_summary,
                     h = BC_h,
                     HSR = BC_HSR,
                     k = BC_k,
                     ISR = BC_ISR,
                     species = "BC",
                     sex = "female")
WBC_VR_treat <- 
  make_treat_matrix(survival_rates_boot_summary = survival_rates_boot_summary,
                    h = WBC_h,
                    HSR = WBC_HSR,
                    k = WBC_k,
                    ISR = WBC_ISR,
                    species = "WBC")

WBC_VR_mprime_male <- 
  make_mprime_matrix(survival_rates_boot_summary = survival_rates_boot_summary,
                     h = WBC_h,
                     HSR = WBC_HSR,
                     k = WBC_k,
                     ISR = WBC_ISR,
                     species = "WBC",
                     sex = "male")

WBC_VR_mprime_female <- 
  make_mprime_matrix(survival_rates_boot_summary = survival_rates_boot_summary,
                     h = WBC_h,
                     HSR = WBC_HSR,
                     k = WBC_k,
                     ISR = WBC_ISR,
                     species = "WBC",
                     sex = "female")

BC_treatment_matrix <- coucal_matrix(BC_VR_treat)
BC_M_prime_matrix_male <- coucal_matrix(BC_VR_mprime_male)
BC_M_prime_matrix_female <- coucal_matrix(BC_VR_mprime_female)

WBC_treatment_matrix <- coucal_matrix(WBC_VR_treat)
WBC_M_prime_matrix_male <- coucal_matrix(WBC_VR_mprime_male)
WBC_M_prime_matrix_female <- coucal_matrix(WBC_VR_mprime_female)

BC_treatment_ASR_analysis <- 
  matrix_ASR(M = BC_treatment_matrix,  
             h = BC_VR_treat$h, 
             k = BC_VR_treat$k,
             HSR = BC_VR_treat$HSR,
             ISR = BC_VR_treat$ISR, 
             immigrant_pop_size = 100,
             iterations = 100, 
             num_boot = 1, 
             iter_add = 1,
             species = "BC")

BC_ASR_treat <- BC_treatment_ASR_analysis$ASR
BC_ASR_treat

WBC_treatment_ASR_analysis <- 
  matrix_ASR(M = WBC_treatment_matrix,  
             h = WBC_VR_treat$h, 
             k = WBC_VR_treat$k,
             HSR = WBC_VR_treat$HSR,
             ISR = WBC_VR_treat$ISR, 
             immigrant_pop_size = 100,
             iterations = 100, 
             num_boot = 1, 
             iter_add = 1,
             species = "BC")

WBC_ASR_treat <- WBC_treatment_ASR_analysis$ASR
WBC_ASR_treat

BC_M_prime_ASR_analysis_male <- 
  matrix_ASR(M = BC_M_prime_matrix_male, 
             h = BC_VR_mprime_male$h, 
             k = BC_VR_mprime_male$k,
             HSR = BC_VR_mprime_male$HSR,
             ISR = BC_VR_mprime_male$ISR,
             immigrant_pop_size = 100,
             iterations = 100,
             num_boot = 1,
             iter_add = 1,
             species = "BC")

BC_ASR_mprime_male <- BC_M_prime_ASR_analysis_male$ASR
BC_ASR_mprime_male

BC_M_prime_ASR_analysis_female <- 
  matrix_ASR(M = BC_M_prime_matrix_female, 
             h = BC_VR_mprime_female$h, 
             k = BC_VR_mprime_female$k,
             HSR = BC_VR_mprime_female$HSR,
             ISR = BC_VR_mprime_female$ISR,
             immigrant_pop_size = 100,
             iterations = 100,
             num_boot = 1,
             iter_add = 1,
             species = "BC")

BC_ASR_mprime_female <- BC_M_prime_ASR_analysis_female$ASR
BC_ASR_mprime_female

WBC_M_prime_ASR_analysis_male <- 
  matrix_ASR(M = WBC_M_prime_matrix_male, 
             h = WBC_VR_mprime_male$h, 
             k = WBC_VR_mprime_male$k,
             HSR = WBC_VR_mprime_male$HSR,
             ISR = WBC_VR_mprime_male$ISR,
             immigrant_pop_size = 100,
             iterations = 100,
             num_boot = 1,
             iter_add = 1,
             species = "BC")

WBC_ASR_mprime_male <- WBC_M_prime_ASR_analysis_male$ASR
WBC_ASR_mprime_male

WBC_M_prime_ASR_analysis_female <- 
  matrix_ASR(M = WBC_M_prime_matrix_female, 
             h = WBC_VR_mprime_female$h, 
             k = WBC_VR_mprime_female$k,
             HSR = WBC_VR_mprime_female$HSR,
             ISR = WBC_VR_mprime_female$ISR,
             immigrant_pop_size = 100,
             iterations = 100,
             num_boot = 1,
             iter_add = 1,
             species = "BC")

WBC_ASR_mprime_female <- WBC_M_prime_ASR_analysis_female$ASR
WBC_ASR_mprime_female

BC_lambda_treat <- BC_treatment_ASR_analysis$lambda
BC_lambda_treat

BC_lambda_mprime_male <- BC_M_prime_ASR_analysis_male$lambda
BC_lambda_mprime_male

BC_lambda_mprime_female <- BC_M_prime_ASR_analysis_female$lambda
BC_lambda_mprime_female

WBC_lambda_treat <- WBC_treatment_ASR_analysis$lambda
WBC_lambda_treat

WBC_lambda_mprime_male <- WBC_M_prime_ASR_analysis_male$lambda
WBC_lambda_mprime_male

WBC_lambda_mprime_female <- WBC_M_prime_ASR_analysis_female$lambda
WBC_lambda_mprime_female

BC_treat_sensitivity_analysis <- 
  sensitivity_analysis(vital_rate_summary = BC_VR_treat, 
                       matrix_str = matrix_structure, 
                       h = BC_VR_treat$h, 
                       k = BC_VR_treat$k, 
                       HSR = BC_VR_treat$HSR, 
                       niter = 1000, 
                       ASR = BC_ASR_treat,
                       lambda = BC_lambda_treat, 
                       ISR = BC_VR_treat$ISR, immigrant_pop_size = 100)

BC_Mprime_sensitivity_analysis_male <- 
  sensitivity_analysis(vital_rate_summary = BC_VR_mprime_male, 
                       matrix_str = matrix_structure, 
                       h = BC_VR_mprime_male$h, 
                       k = BC_VR_mprime_male$k, 
                       HSR = BC_VR_mprime_male$HSR, 
                       niter = 1000, 
                       ASR = BC_ASR_mprime_male,
                       lambda = BC_lambda_mprime_male, 
                       ISR = BC_VR_mprime_male$ISR, immigrant_pop_size = 100)

BC_Mprime_sensitivity_analysis_female <- 
  sensitivity_analysis(vital_rate_summary = BC_VR_mprime_female, 
                       matrix_str = matrix_structure, 
                       h = BC_VR_mprime_female$h, 
                       k = BC_VR_mprime_female$k, 
                       HSR = BC_VR_mprime_female$HSR, 
                       niter = 1000, 
                       ASR = BC_ASR_mprime_female,
                       lambda = BC_lambda_mprime_female, 
                       ISR = BC_VR_mprime_female$ISR, immigrant_pop_size = 100)

WBC_treat_sensitivity_analysis <- 
  sensitivity_analysis(vital_rate_summary = WBC_VR_treat, 
                       matrix_str = matrix_structure, 
                       h = WBC_VR_treat$h, 
                       k = WBC_VR_treat$k, 
                       HSR = WBC_VR_treat$HSR, 
                       niter = 1000, 
                       ASR = WBC_ASR_treat,
                       lambda = WBC_lambda_treat, 
                       ISR = WBC_VR_treat$ISR, immigrant_pop_size = 100)

WBC_Mprime_sensitivity_analysis_male <- 
  sensitivity_analysis(vital_rate_summary = WBC_VR_mprime_male, 
                       matrix_str = matrix_structure, 
                       h = WBC_VR_mprime_male$h, 
                       k = WBC_VR_mprime_male$k, 
                       HSR = WBC_VR_mprime_male$HSR, 
                       niter = 1000, 
                       ASR = WBC_ASR_mprime_male,
                       lambda = WBC_lambda_mprime_male, 
                       ISR = WBC_VR_mprime_male$ISR, immigrant_pop_size = 100)

WBC_Mprime_sensitivity_analysis_female <- 
  sensitivity_analysis(vital_rate_summary = WBC_VR_mprime_female, 
                       matrix_str = matrix_structure, 
                       h = WBC_VR_mprime_female$h, 
                       k = WBC_VR_mprime_female$k, 
                       HSR = WBC_VR_mprime_female$HSR, 
                       niter = 1000, 
                       ASR = WBC_ASR_mprime_female,
                       lambda = WBC_lambda_mprime_female, 
                       ISR = WBC_VR_mprime_female$ISR, immigrant_pop_size = 100)
# ```
#conduct the LTRE comparing the two matrices
# ```{r}
BC_LTRE_male <- 
  LTRE_analysis(Mprime_sens = BC_Mprime_sensitivity_analysis_male, 
                matrix_str = matrix_str, 
                vital_rates = BC_VR_treat,
                species_name = "BC",
                sex = "male")

BC_LTRE_female <- 
  LTRE_analysis(Mprime_sens = BC_Mprime_sensitivity_analysis_female, 
                matrix_str = matrix_str, 
                vital_rates = BC_VR_treat,
                species_name = "BC",
                sex = "female")

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