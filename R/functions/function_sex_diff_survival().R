source("R/project/project_libraries.R")

# sex_diff_survival() takes the list of results from the bootstrap proceedure
# and calculates the sex difference in survival at each life stage for each
# iteration

sex_diff_survival <- function(boot_out_list, niter) {
  
  # make an empty datarame to store the results
  sex_diff_surv_output <- data.frame(Adult = numeric(niter * 4),
                                     Fledgling = numeric(niter * 4),
                                     Groundling = numeric(niter * 4),
                                     Nestling = numeric(niter * 4))
  
  # for loop to go through each iteration and calculate the differece between 
  # female and male survival rates for each stage.
  for(i in 1:(niter * 4)){
    Adult <- 
      boot_out_list$survival_rates_boot[which(boot_out_list$survival_rates_boot$iter == i), 2][9] -
      boot_out_list$survival_rates_boot[which(boot_out_list$survival_rates_boot$iter == i), 2][8]
    Fledgling <- 
      boot_out_list$survival_rates_boot[which(boot_out_list$survival_rates_boot$iter == i), 2][7] -
      boot_out_list$survival_rates_boot[which(boot_out_list$survival_rates_boot$iter == i), 2][6]
    Groundling <- 
      boot_out_list$survival_rates_boot[which(boot_out_list$survival_rates_boot$iter == i), 2][5] -
      boot_out_list$survival_rates_boot[which(boot_out_list$survival_rates_boot$iter == i), 2][4]
    Nestling <- 
      boot_out_list$survival_rates_boot[which(boot_out_list$survival_rates_boot$iter == i), 2][3] -
      boot_out_list$survival_rates_boot[which(boot_out_list$survival_rates_boot$iter == i), 2][2]
    
    sex_diff_surv_output[i, 1] <- Adult
    sex_diff_surv_output[i, 2] <- Fledgling
    sex_diff_surv_output[i, 3] <- Groundling
    sex_diff_surv_output[i, 4] <- Nestling
    
  }
  
  # restructure the output and lable columns
  sex_diff_surv_output <- reshape2::melt(data = sex_diff_surv_output)
  colnames(sex_diff_surv_output) <- c("stage", "difference")
  
  # return the output
  sex_diff_surv_output
}