# load packages
source("scripts/01_libraries.R")

# load functions
function.sources = list.files(path = "scripts", 
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
sapply(function.sources, source, .GlobalEnv)

# load capture histories
data.sources = list.files(path = "cooked_data", 
                          pattern="*ch.rds$", full.names = TRUE, 
                          ignore.case = TRUE)
sapply(data.sources, load, .GlobalEnv)

niter = 250
set.seed(87)

# run bootstrap procedure on White-browed Coucals
WBC_survival_ASR_bootstrap_result_four <-
  pbsapply(1:niter, run_bootstrap_survival_ASR, 
           adult = White_browed_Coucal_adult_CJS_ch,
           fledgling = White_browed_Coucal_fledgling_Burnham_ch,
           nestling = White_browed_Coucal_nestling_Known_ch,
           k = 4, 
           HSR = 0.5198, 
           h = 1/1.1, 
           egg_survival = 0.18, 
           ISR = 0.524, 
           immigrant_pop_size = 100, 
           flight_age = 32,
           first_year = 2005,
           bootstrap_name = "WBC_boot_four",
           species = "WBC",
           iter_add = 4,
           prefix_number = "boot_four_")

save(WBC_survival_ASR_bootstrap_result_four,
     file = "output/bootstraps/WBC_survival_ASR_bootstrap_result_four.rds")