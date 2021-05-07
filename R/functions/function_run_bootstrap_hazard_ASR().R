source("R/project/project_libraries.R")
source("R/functions/function_bootstrap_hazard_data().R")
source("R/functions/function_bootstrap_survival_ASR().R")

# run_bootstrap_survival_ASR() initiates the bootstrap_data() and 
# bootstrap_survival_ASR() functions in sequence.

run_bootstrap_hazard_ASR <- function(num_boot, offspring,
                                     k, HSR, h, egg_survival, flight_age,
                                     ISR, immigrant_pop_size, fledge_age, 
                                     bootstrap_name, adult_survival_rate,
                                     species, iter_add, prefix_number,
                                     alpha_value = 1.4,
                                     max_time = 70)
{
  # run the sampling function and specify the datasets
  bootstrap_data_list <- 
    bootstrap_hazard_data(offspring = offspring, 
                          num_boot = num_boot,
                          species = species, 
                          iter_add = iter_add,
                          alpha_value = alpha_value,
                          max_time = max_time)
  
  # run the survival analysis and ASR deduction on the sampled data
  result <- bootstrap_hazard_ASR(coucal_boot_list = bootstrap_data_list, 
                                 num_boot = num_boot, 
                                 egg_survival = egg_survival, 
                                 adult_survival_rate = adult_survival_rate,
                                 ISR = ISR, 
                                 immigrant_pop_size = immigrant_pop_size, 
                                 HSR = HSR, h = h, k = k, 
                                 fledge_age = fledge_age,
                                 flight_age = flight_age, 
                                 bootstrap_name = bootstrap_name,
                                 species = species,
                                 iter_add = iter_add,
                                 prefix_number = prefix_number)
  result
}