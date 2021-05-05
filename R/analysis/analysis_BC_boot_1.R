# load packages
source("scripts/01_libraries.R")

# load functions
function.sources = list.files(path = "scripts",
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
sapply(function.sources, source)

# load capture histories
data.sources = list.files(path = "cooked_data", 
                          pattern="*ch.rds$", full.names = TRUE, 
                          ignore.case = TRUE)
sapply(data.sources, load, .GlobalEnv)

niter = 250
set.seed(14)

# run bootstrap procedure on Black Coucals
BC_survival_ASR_bootstrap_result_one <-
  pbsapply(1:niter, run_bootstrap_survival_ASR,
           adult = Black_Coucal_adult_CJS_ch,
           fledgling = Black_Coucal_fledgling_Burnham_ch,
           nestling = Black_Coucal_nestling_Known_ch,
           k = 4,
           HSR = 0.4955,
           h = 1/2.9,
           egg_survival = 0.32,
           ISR = 0.738,
           immigrant_pop_size = 100,
           flight_age = 36,
           first_year = 2001,
           bootstrap_name = "BC_boot_one",
           species = "BC",
           iter_add = 1,
           prefix_number = "boot_one_")

save(BC_survival_ASR_bootstrap_result_one,
     file = "output/bootstraps/BC_survival_ASR_bootstrap_result_one.rds")
