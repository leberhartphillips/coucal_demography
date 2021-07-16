source("R/project/project_libraries.R")
source("R/functions/function_bootstrap_hazard_data().R")
source("R/functions/function_bootstrap_juv_hazd_ad_surv_ASR().R")

# run_bootstrap_survival_ASR() initiates the bootstrap_data() and 
# bootstrap_survival_ASR() functions in sequence.

run_bootstrap_juv_hazd_ad_surv_ASR <- function(num_boot, 
                                               offspring,
                                               k_dist, 
                                               HSR_dist, 
                                               h_dist, 
                                               egg_surv_dist, 
                                               flight_age_distF,
                                               flight_age_distM,
                                               ISR_dist, 
                                               immigrant_pop_size, 
                                               fledge_age_distF, 
                                               fledge_age_distM,
                                               bootstrap_name, 
                                               adult_surival_boot_out,
                                               species, 
                                               iter_add, 
                                               prefix_number,
                                               alpha_value = 1.4,
                                               max_time = 70,
                                               niter)
{
  # run the sampling function and specify the datasets
  bootstrap_data_list <- 
    bootstrap_hazard_data(offspring = offspring, 
                          num_boot = num_boot,
                          species = species, 
                          iter_add = iter_add,
                          alpha_value = alpha_value,
                          max_time = max_time,
                          niter = niter)
  
  # run the survival analysis and ASR deduction on the sampled data
  result <- 
    bootstrap_juv_hazd_ad_surv_ASR(coucal_boot_list = bootstrap_data_list, 
                                   num_boot = num_boot, 
                                   egg_surv_dist = egg_surv_dist, 
                                   adult_surival_boot_out = adult_surival_boot_out,
                                   ISR_dist = ISR_dist, 
                                   immigrant_pop_size = immigrant_pop_size, 
                                   HSR_dist = HSR_dist, 
                                   h_dist = h_dist, 
                                   k_dist = k_dist, 
                                   fledge_age_distF = fledge_age_distF,
                                   flight_age_distF = flight_age_distF, 
                                   fledge_age_distM = fledge_age_distM,
                                   flight_age_distM = flight_age_distM, 
                                   bootstrap_name = bootstrap_name,
                                   species = species,
                                   iter_add = iter_add,
                                   prefix_number = prefix_number)
  result
}
