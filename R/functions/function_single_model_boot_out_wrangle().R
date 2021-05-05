# single_model_boot_out_wrangle() loads all the output rds files from the 
# run_known_fate_boot() output and consolidates them into a single object for 
# plotting, etc.
single_model_boot_out_wrangle <- 
  function(species, stage, niter, output_dir, rds_file){
    
    # specify directory string
    boot_path <- paste0(output_dir, species, rds_file, stage, ".rds")
    
    # load each of the four bootstrap output files
    bootstrap_out <- readRDS(file = boot_path)

    # bind all the fledgling daily survival rates output together
    reals_known_fate_boot <-
      do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                            function(x) bootstrap_out["known_fate_model_output", x]$known_fate_model_output$reals)) %>%
      arrange(iter, sex, age)
    
    # bind all the fledgling beta estimates output together
    betas_known_fate_boot <-
      do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                            function(x) bootstrap_out["known_fate_model_output", x]$known_fate_model_output$betas))
    
    # tidy all the bound data together into a mega list of results
    boot_out_tidy <- 
      list(reals_known_fate_boot = reals_known_fate_boot,
           betas_known_fate_boot = betas_known_fate_boot)
    
    boot_out_tidy
  }