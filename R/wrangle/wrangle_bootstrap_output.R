# load packages
source("scripts/01_libraries.R")

# load functions
function.sources = list.files(path = "scripts", 
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
sapply(function.sources, source, .GlobalEnv)

# clean up the output from the bootstrap procedure and save as rds
BC_boot_out <- 
  boot_out_wrangle(species = "BC", niter = 250) %>% 
  AIC_table_wrangle()

str(BC_boot_out$ASR_boot)
sum(is.na(BC_boot_out$ASR_boot$ASR_estimate))

WBC_boot_out <- 
  boot_out_wrangle(species = "WBC", niter = 250) %>% 
  AIC_table_wrangle()

# save cleaned bootstrap output
save(BC_boot_out, 
     file = "output/wrangled/Black_Coucal_bootstrap_result_clean.rds")
save(WBC_boot_out, 
     file = "output/wrangled/White-browed_Coucal_bootstrap_result_clean.rds")
