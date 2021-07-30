# make_mprime_and_treat_summary_JSR() makes the M-prime life table (i.e., a 
# table with all rates halfway between the treatment table and the unbiased 
# table (either unbiased with all female rates == male rates (when base_sex is
# "Male") or unbiased with all males rates == female rates (when base_sex is
# "Female"))

make_mprime_and_treat_summary_JSR <- 
  function(survival_rates_boot_summary, species_name, 
           h_rate, k_rate, HSR_rate, ISR_rate, egg_S_rate, imm_N_base = 100,
           base_sex = "Male", mprime = TRUE){
    
    survival_rates_LTRE_mprime <-
      survival_rates_boot_summary %>% 
      filter(species == species_name) %>% 
      pivot_wider(names_from = sex, values_from = fit) %>% 
      mutate(mprime = (Male + Female)/2) %>% 
      pivot_longer(cols = c(Male, Female), names_to = "sex", values_to = "obs") %>% 
      arrange(sex, age) %>% 
      dplyr::rename(parameter = age) %>% 
      mutate(parameter = paste(sex, "age", parameter, sep = "_")) %>% 
      dplyr::select(species, parameter, obs, mprime)
    
    other_rates_LTRE_mprime <- 
      data.frame(species = species_name,
                 parameter = c("h", "k", "HSR", "ISR", "egg_S", "imm_N"),
                 obs = c(h_rate, k_rate, HSR_rate, ISR_rate, egg_S_rate, imm_N_base)) %>% 
      mutate(mprime = ifelse(parameter == "h", (obs + 1) / 2,
                             ifelse(parameter == "k", obs,
                                    ifelse(parameter == "HSR", (obs + 0.5) / 2,
                                           ifelse(parameter == "ISR", (obs + 0.5) / 2,
                                                  ifelse(parameter == "egg_S", obs,
                                                         ifelse(parameter == "imm_N", obs, 0)))))))
    
    LTRE_mprime <- 
      bind_rows(survival_rates_LTRE_mprime, other_rates_LTRE_mprime)
    
    # make M-prime life table (i.e., if "Male", then the female vital rates are
    # set to halfway between male and female rates, and male rates stay as observed.
    # If "Female", then male rates are set to halfway between female and male 
    # rates, and female rates stay as observed)
    if(mprime == TRUE){
      if(base_sex == "Male"){
        LTRE_mprime <- 
          LTRE_mprime %>% 
          mutate(obs = ifelse(str_detect(parameter, "Male", negate = TRUE), mprime, obs)) %>% 
          dplyr::select(!mprime)
      }
      else{
        LTRE_mprime <- 
          LTRE_mprime %>% 
          mutate(obs = ifelse(str_detect(parameter, "Female", negate = TRUE), mprime, obs)) %>% 
          dplyr::select(!mprime)
      }
    }
    
    # make treatment life table (i.e., all the observed sex-specific rates)
    else{
      LTRE_mprime <- 
        LTRE_mprime %>% 
        dplyr::select(!mprime)
    }
    
  }
