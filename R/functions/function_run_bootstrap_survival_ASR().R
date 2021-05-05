source("scripts/01_libraries.R")
source("scripts/02_func_bootstrap_data().R")
source("scripts/02_func_bootstrap_survival_ASR().R")

# run_bootstrap_survival_ASR() initiates the bootstrap_data() and 
# bootstrap_survival_ASR() functions in sequence.

run_bootstrap_survival_ASR <- function(num_boot, adult, fledgling, nestling,
                                       k, HSR, h, egg_survival, flight_age,
                                       ISR, immigrant_pop_size, 
                                       first_year, bootstrap_name,
                                       species, iter_add, prefix_number)
{
  # run the sampling function and specify the datasets
  bootstrap_data_list <- bootstrap_data(fledgling = fledgling, 
                                        nestling = nestling, 
                                        adult = adult,
                                        num_boot = num_boot,
                                        species = species, 
                                        iter_add = iter_add)
  
  # run the survival analysis and ASR deduction on the sampled data
  result <- bootstrap_survival_ASR(coucal_boot_list = bootstrap_data_list, 
                                   num_boot = num_boot, 
                                   egg_survival = egg_survival, 
                                   ISR = ISR, 
                                   immigrant_pop_size = immigrant_pop_size, 
                                   HSR = HSR, h = h, k = k, 
                                   flight_age = flight_age, 
                                   first_year = first_year,
                                   bootstrap_name = bootstrap_name,
                                   species = species,
                                   iter_add = iter_add,
                                   prefix_number = prefix_number)
  result
}