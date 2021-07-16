source("R/project/project_libraries.R")

# sex_diff_survival() takes the list of results from the bootstrap proceedure
# and calculates the sex difference in survival at each life stage for each
# iteration

sex_diff_hazard <- function(boot_out_list, niter) {
  
  # make an empty datarame to store the results
  sex_diff_surv_output <- data.frame(Adult = niter,
                                     Fledgling = niter,
                                     Groundling = niter,
                                     Nestling = niter,
                                     ISR = niter,
                                     HSR = niter,
                                     h = niter)
  
  # for loop to go through each iteration and calculate the difference between 
  # female and male survival rates for each stage.
  for(i in 1:niter){
    Adult <- 
      pull(filter(boot_out_list$vital_rate_ests_boot[which(boot_out_list$vital_rate_ests_boot$iter == i), ], 
                  stage == "adult" & rate == "survival" & sex == "Female"), value) -
      pull(filter(boot_out_list$vital_rate_ests_boot[which(boot_out_list$vital_rate_ests_boot$iter == i), ], 
                  stage == "adult" & rate == "survival" & sex == "Male"), value)
    Fledgling <- 
      pull(filter(boot_out_list$vital_rate_ests_boot[which(boot_out_list$vital_rate_ests_boot$iter == i), ], 
                  stage == "fledgling" & rate == "survival" & sex == "Female"), value) -
      pull(filter(boot_out_list$vital_rate_ests_boot[which(boot_out_list$vital_rate_ests_boot$iter == i), ], 
                  stage == "fledgling" & rate == "survival" & sex == "Male"), value)
    Groundling <- 
      pull(filter(boot_out_list$vital_rate_ests_boot[which(boot_out_list$vital_rate_ests_boot$iter == i), ], 
                  stage == "groundling" & rate == "survival" & sex == "Female"), value) -
      pull(filter(boot_out_list$vital_rate_ests_boot[which(boot_out_list$vital_rate_ests_boot$iter == i), ], 
                  stage == "groundling" & rate == "survival" & sex == "Male"), value)
    Nestling <- 
      pull(filter(boot_out_list$vital_rate_ests_boot[which(boot_out_list$vital_rate_ests_boot$iter == i), ], 
                  stage == "nestling" & rate == "survival" & sex == "Female"), value) -
      pull(filter(boot_out_list$vital_rate_ests_boot[which(boot_out_list$vital_rate_ests_boot$iter == i), ], 
                  stage == "nestling" & rate == "survival" & sex == "Male"), value)
    ISR <- 
      0.5 - pull(filter(boot_out_list$vital_rate_ests_boot[which(boot_out_list$vital_rate_ests_boot$iter == i), ], 
                  stage == "ISR"), value)
    HSR <- 
      0.5 - pull(filter(boot_out_list$vital_rate_ests_boot[which(boot_out_list$vital_rate_ests_boot$iter == i), ], 
                  stage == "HSR"), value)
    h <- 
      1 - pull(filter(boot_out_list$vital_rate_ests_boot[which(boot_out_list$vital_rate_ests_boot$iter == i), ], 
                  stage == "h"), value)
    
    sex_diff_surv_output[i, 1] <- ifelse(length(Adult) > 0, Adult, NA)
    sex_diff_surv_output[i, 2] <- ifelse(length(Fledgling) > 0, Fledgling, NA)
    sex_diff_surv_output[i, 3] <- ifelse(length(Groundling) > 0, Groundling, NA)
    sex_diff_surv_output[i, 4] <- ifelse(length(Nestling) > 0, Nestling, NA)
    sex_diff_surv_output[i, 5] <- ifelse(length(ISR) > 0, ISR, NA)
    sex_diff_surv_output[i, 6] <- ifelse(length(HSR) > 0, HSR, NA)
    sex_diff_surv_output[i, 7] <- ifelse(length(h) > 0, h, NA)
    
    
  }
  
  # restructure the output and lable columns
  sex_diff_surv_output <- reshape2::melt(data = sex_diff_surv_output)
  colnames(sex_diff_surv_output) <- c("stage", "difference")
  
  # return the output
  sex_diff_surv_output
}
