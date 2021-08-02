# load packages
source("R/project/project_libraries.R")

# load functions
function.sources = list.files(path = "R/functions",
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
sapply(function.sources, source)

# # load tidy output
# WBC_hazard_rate_boot_tidy <-
#   readRDS("output/bootstraps/hazard/cooked/WBC_haz_sur_ASR_boot_tidy_stoc_no_imm.rds")

# LTRE analysis
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
