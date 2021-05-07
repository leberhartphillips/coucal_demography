source("R/project/project_libraries.R")

# bootstrap_hazard_data() randomly samples one sibling per nest from the 
# offspring survival dataset
bootstrap_hazard_data <- 
  function(offspring, num_boot, species, iter_add, alpha_value, max_time) {
  
    
    # set attempt to 0 at start of each loop
    attempt <- 0
    
    max_time <- max_time
    
    # sample a new offspring dataset containing only one nest 
    # member per draw
    offspring_boot <- 
      offspring %>% 
      group_by(nest_ID) %>% 
      sample_n(1)
    
    # store simulated estimates only if peak >= 1 and <= 10 and it's less than
    # 100 attempts
    while( (exists("hzd_curve_try_M") == FALSE | exists("hzd_curve_try_F") == FALSE) && attempt <= max_time) {
      
      time_vector <- seq(0, max_time - attempt, 1)
      
      # next attempt
      attempt <- attempt + 1
      
      try(
        hzd_ss_try_M <- sshzd(Surv(exit, event, entry) ~ exit, 
                            data = filter(offspring_boot, sex == "M"), 
                            alpha = alpha_value),
        silent = TRUE
      )
      
      try(
        hzd_ss_try_F <- sshzd(Surv(exit, event, entry) ~ exit, 
                            data = filter(offspring_boot, sex == "F"), 
                            alpha = alpha_value),
        silent = TRUE
      )
      
      # simulate an estimate
      try(
        hzd_curve_try_M <- 
          hzdcurve.sshzd(object = hzd_ss_try_M, time = time_vector, se = TRUE),
        silent = TRUE
      )
      try(
        hzd_curve_try_F <- 
          hzdcurve.sshzd(object = hzd_ss_try_F, time = time_vector, se = TRUE),
        silent = TRUE
      )
      
    }
    
    # make a list of these two datasets, which will be used in the next function
    out <- list(offspring_boot = offspring_boot, 
                iter = num_boot + ((iter_add - 1) * niter),
                species = species,
                time_vector = time_vector)
  }