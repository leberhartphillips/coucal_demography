# load libraries
source("R/project/project_libraries.R")

# load functions
function.sources = list.files(path = "R/functions", 
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
sapply(function.sources, source, .GlobalEnv)

# load capture histories
data.sources = list.files(path = "data/cooked", 
                          pattern="*ch.rds$", full.names = TRUE, 
                          ignore.case = TRUE)
sapply(data.sources, load, .GlobalEnv)

# assign capture histories
niter = 1000

# run bootstrap procedure on Black Coucals
BC_known_fate_single_mod_boot_nest <-
  pbsapply(1:niter, run_known_fate_boot, 
           ch_input = Black_Coucal_nestling_Known_ch,
           bootstrap_name = "BC_boot_nest_known_fate",
           species = "BC",
           iter_add = 1,
           prefix_number = "BC_boot_nest_known_fate",
           stage_name = "nest")

# run bootstrap procedure on Black Coucals
BC_known_fate_single_mod_boot_fled <-
  pbsapply(1:niter, run_known_fate_boot, 
           ch_input = Black_Coucal_fledgling_Burnham_ch,
           bootstrap_name = "BC_boot_fled_known_fate",
           species = "BC",
           iter_add = 1,
           prefix_number = "BC_boot_fled_known_fate",
           stage_name = "fled")

# save model output
saveRDS(object = BC_known_fate_single_mod_boot_nest, 
        file = "output/bootstraps/single_models/raw/BC_known_fate_single_mod_boot_nest.rds")
saveRDS(object = BC_known_fate_single_mod_boot_fled, 
        file = "output/bootstraps/single_models/raw/BC_known_fate_single_mod_boot_fled.rds")