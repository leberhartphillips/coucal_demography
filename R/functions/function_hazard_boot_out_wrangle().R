# single_model_boot_out_wrangle() loads all the output rds files from the 
# run_known_fate_boot() output and consolidates them into a single object for 
# plotting, etc.
hazard_boot_out_wrangle <- 
  function(species, niter, output_dir, rds_file){
    
    # specify directory string
    boot_path <- paste0(output_dir, species, rds_file, ".rds")
    
    # load each of the four bootstrap output files
    bootstrap_out <- readRDS(file = boot_path)
    
    # bind all the fledgling daily survival rates output together
    hazard_rates_boot <-
      do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                            function(x) bootstrap_out["offspring_hazard_rates", x]$offspring_hazard_rates)) %>%
      arrange(iter, sex, age)
    
    # bind all the ASR estimates output together
    ASR_ests_boot <-
      do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                            function(x) bootstrap_out["ASR_SSD", x]$ASR_SSD$ASR)) %>% 
      as.data.frame() %>% 
      mutate(species = species)
    
    # bind all the ASR estimates output together
    vital_rate_ests_boot <-
      do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                            function(x) bootstrap_out["coucal_vital_rates", x]$coucal_vital_rates))
    
    # tidy all the bound data together into a mega list of results
    boot_out_tidy <- 
      list(hazard_rates_boot = hazard_rates_boot,
           vital_rate_ests_boot = vital_rate_ests_boot,
           ASR_ests_boot = ASR_ests_boot)
    
    boot_out_tidy
  }