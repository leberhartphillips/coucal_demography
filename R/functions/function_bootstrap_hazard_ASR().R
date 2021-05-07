source("R/project/project_libraries.R")
source("R/functions/function_matrix_ASR().R")
source("R/functions/function_coucal_matrix().R")

# bootstrap_hazard_ASR() extracts the stage and sex-specific hazard rates from the 
# bootstrapped sample and feeds them into the matrix model to estimates the ASR.
# coucal_boot_list :  the output list from bootstrap_data()
# ISR : the immigration sex ratio (i.e., the observed adult sex ratio during the
# breeding season)
# immigrant_pop_size : the number of immigrants entering the population at the 
# start of each season
# default is a vector with 10 individuals in each stage)
# h : harem size index
# k : clutch size
# HSR : hatching sex ratio
# num_boot : bootstrap number in the loop (do not specify)
# species : coucal species name
# iter_add : addition to bootstrap number in the loop (do not specify)
# bootstrap_name : ID of bootstrap iteration (do not specify)
# flight_age : age at first flight in days
# first_year : the first year of the study period (e.g., 2001)
# egg_survival : hatching probability

bootstrap_hazard_ASR <- 
  function(coucal_boot_list, num_boot, iter_add,
           egg_survival,
           adult_survival_rate,
           ISR,
           immigrant_pop_size = 100,
           HSR = 0.5,
           h,
           k,
           fledge_age,
           flight_age,
           first_year,
           bootstrap_name,
           species,
           prefix_number,
           alpha_value = 1.4) {
    
    # specify the bootstrapped data samples (from the previous function)
    offspring_data <- coucal_boot_list[["offspring_boot"]]
    
    # clean up capture histories
    offspring_data <-    
      offspring_data %>% 
      ungroup() %>% 
      as.data.frame()
    
    # fit smoothed spline of hazard function for either sex
    M_haz_ss <- sshzd(Surv(exit, event, entry) ~ exit, 
                      data = filter(offspring_data, sex == "M"), 
                      alpha = alpha_value)

    F_haz_ss <- sshzd(Surv(exit, event, entry) ~ exit, 
                      data = filter(offspring_data, sex == "F"), 
                      alpha = alpha_value)

    haz_ss_function <- list(Male_haz_ss = M_haz_ss,
                            Female_haz_ss = F_haz_ss)
    
    # extract fitted estimates from the spline function
    M_haz_ss_curve <- 
      hzdcurve.sshzd(object = M_haz_ss, time = coucal_boot_list[["time_vector"]], se = TRUE)
    
    F_haz_ss_curve <- 
      hzdcurve.sshzd(object = F_haz_ss, time = coucal_boot_list[["time_vector"]], se = TRUE)
    
    haz_ss_curve <- 
      expand.grid(species = species, 
                  age = coucal_boot_list[["time_vector"]],
                  sex = c("Male", "Female")) %>% 
      mutate(fit = c(M_haz_ss_curve$fit, F_haz_ss_curve$fit),
             se = c(M_haz_ss_curve$se, F_haz_ss_curve$se)) %>% 
      mutate(estimate = 1 - fit,
             upper = 1 - fit * exp(1.96 * se),
             lower = 1 - fit / exp(1.96 * se),
             iter = num_boot)
    
    # transform the daily nestling survival (DCS) to apparent fledgling success
    # by calculating the product of all DCS estimates:
    coucal_nestling_survival <-
      haz_ss_curve %>% 
      filter(age <= fledge_age) %>% 
      group_by(sex) %>% 
      dplyr::summarise(value = prod(estimate), .groups = 'drop') %>%
      mutate(stage = "nestling",
             rate = "survival")
    
    coucal_groundling_survival <-
      haz_ss_curve %>% 
      filter(age < flight_age & age > fledge_age) %>% 
      group_by(sex) %>% 
      dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
      mutate(stage = "groundling",
             rate = "survival")
    
    coucal_fledgling_survival <-
      haz_ss_curve %>% 
      filter(age >= flight_age) %>% 
      group_by(sex) %>% 
      dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
      mutate(stage = "fledgling",
             rate = "survival")
    
    coucal_adult_survival <-
      expand.grid(sex = c("Male", "Female"),
                  value = adult_survival_rate,
                  stage = c("adult"),
                  rate = c("survival"))
    
    adult_F_immigrants <- immigrant_pop_size * (1 - ISR)
    adult_M_immigrants <- immigrant_pop_size * ISR
    
    coucal_adult_immigration <- 
      data.frame(sex = c("Female", "Male"),
                 value = c(adult_F_immigrants, adult_M_immigrants),
                 stage = c("adult"),
                 rate = c("immigration"))
    
    coucal_egg_survival <- 
      data.frame(sex = NA,
                 value = egg_survival,
                 stage = c("egg"),
                 rate = c("survival"))
    
    # Bind the juvenile and adult dataframe with the nestlings
    coucal_vital_rates <- 
      bind_rows(coucal_egg_survival,
                coucal_nestling_survival,
                coucal_groundling_survival,
                coucal_fledgling_survival,
                coucal_adult_survival,
                coucal_adult_immigration) %>% 
      as.data.frame() %>% 
      mutate(iter = num_boot, #+ ((iter_add - 1) * niter),
             species = species) %>% 
      arrange(sex, stage, rate)
    
    # Create a list of demographic rates from the survival analyses above
    demographic_rates <- list(Egg_survival = coucal_vital_rates[11, 2],
                              F_Nestling_survival = coucal_vital_rates[5, 2],
                              F_Groundling_survival = coucal_vital_rates[4, 2],
                              F_Fledgling_survival = coucal_vital_rates[3, 2],
                              F_Adult_survival = coucal_vital_rates[2, 2],
                              F_Adult_immigration = coucal_vital_rates[1, 2],
                              M_Nestling_survival = coucal_vital_rates[10, 2],
                              M_Groundling_survival = coucal_vital_rates[9, 2],
                              M_Fledgling_survival = coucal_vital_rates[8, 2],
                              M_Adult_survival = coucal_vital_rates[7, 2],
                              M_Adult_immigration = coucal_vital_rates[6, 2],
                              
                              # Define hatching sex ratio
                              HSR = HSR,
                              
                              # Define the mating system (h), and clutch size (k)
                              h = h,
                              k = k)
    
    # Build matrix based on rates specified in the list above
    demographic_matrix <- coucal_matrix(demographic_rates)
    
    # populate sex-specific adult immigration matrix
    immigration_vector <- matrix(data = c(0, demographic_rates$F_Adult_immigration, 
                                          0, demographic_rates$M_Adult_immigration), 
                                 nrow = 4, ncol = 1)
    
    # Determine the ASR at the stable stage distribution
    ASR_SSD <- matrix_ASR(M = demographic_matrix,
                          h = demographic_rates$h,
                          HSR = demographic_rates$HSR, iterations = 100,
                          num_boot = num_boot,
                          species = species,
                          immigrant_pop_size = immigrant_pop_size,
                          ISR = ISR, iter_add = 1)
    
    # Extract ASR
    ASR_estimate <- ASR_SSD$ASR
    
    # make a list of all the results from this iteration
    bootstrap_results_list <- 
      list(#offspring_hazard_function = haz_ss_function, 
           offspring_hazard_rates = haz_ss_curve,
           coucal_vital_rates = coucal_vital_rates,
           ASR_SSD = ASR_SSD,
           bootstrapped_data = coucal_boot_list)
    
    # extract file directory
    file_directory <- paste0(getwd(), "/output/bootstraps/hazard/raw/", bootstrap_name, "_", num_boot, ".Rds")
    
    # save the results of the iteration
    save(bootstrap_results_list, file = file_directory)
    
    bootstrap_results_list
  }