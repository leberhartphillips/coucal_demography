source("R/project/project_libraries.R")

# JSR_bootstrap takes the bootstrapped hazard functions and simulates the
# age-specific population dynamics of juveniles after hatching while also
# taking into consideration the hatching sex ratio and mating function
JSR_bootstrap_deter <- 
  function(VR_dataframe,
           immigrant_pop_size = 100,
           species){
    
    surv_rates <- 
      filter(VR_dataframe, str_detect(parameter, "_age_")) %>% 
      mutate(sex = ifelse(str_detect(parameter, "Female"), "Female", "Male"),
             age = readr::parse_number(parameter))
    
    ISR_rate <- filter(VR_dataframe, parameter == "ISR") %>% pull(obs)
    
    h_rate <- filter(VR_dataframe, parameter == "h") %>% pull(obs)
    
    HSR_rate <- filter(VR_dataframe, parameter == "HSR") %>% pull(obs)
    
    k_rate <- filter(VR_dataframe, parameter == "k") %>% pull(obs)
    
    egg_survival_rate  <- filter(VR_dataframe, parameter == "egg_S") %>% pull(obs)

    n <- length(unique(surv_rates$age))
    
    F_juv_pop_storage_boot <- 
      matrix(, nrow = n, ncol = 1)
    
    M_juv_pop_storage_boot <- 
      matrix(, nrow = n, ncol = 1)
    
    N_F_immigrants <- 
      immigrant_pop_size * (1 - ISR_rate)
    
    N_M_immigrants <- 
      immigrant_pop_size * ISR_rate
    
    # Female freq-dep fecundity of Female chicks
    F_F_eggs <- 
      ((k_rate * N_M_immigrants) / (N_M_immigrants + (N_F_immigrants / h_rate)) * (1 - HSR_rate))
    
    # Female freq-dep fecundity of Male chicks
    F_M_eggs <- 
      ((k_rate * N_M_immigrants) / (N_M_immigrants + (N_F_immigrants / h_rate)) * HSR_rate)
    
    # Male freq-dep fecundity of Female chicks
    M_F_eggs <- 
      ((k_rate * N_F_immigrants) / (N_M_immigrants + (N_F_immigrants / h_rate)) * (1 - HSR_rate))
    
    # Male freq-dep fecundity of Male chicks
    M_M_eggs <- 
      ((k_rate * N_F_immigrants) / (N_M_immigrants + (N_F_immigrants / h_rate)) * HSR_rate)
    
    N_F_eggs <- F_F_eggs * N_F_immigrants + M_F_eggs * N_M_immigrants
    N_M_eggs <- F_M_eggs * N_F_immigrants + M_M_eggs * N_M_immigrants
    
    F_juv_pop_storage_boot[1, 1] <- N_F_eggs * egg_survival_rate
    M_juv_pop_storage_boot[1, 1] <- N_M_eggs * egg_survival_rate
    
    for(i in 2:n){
      F_juv_pop_storage_boot[i, 1] <- 
        F_juv_pop_storage_boot[i - 1, 1] * surv_rates[which(surv_rates$sex == "Female"), ][["obs"]][i - 1]
      M_juv_pop_storage_boot[i, 1] <- 
        M_juv_pop_storage_boot[i - 1, 1] * surv_rates[which(surv_rates$sex == "Male"), ][["obs"]][i - 1]
    }
    
    F_juv_pop_storage_boot_clean <-
      data.frame(F_juv_pop_storage_boot) %>% 
      dplyr::rename(pop_size = F_juv_pop_storage_boot) %>% 
      mutate(age = c(0:(n-1))) %>% 
      mutate(sex = "F", 
             species = species)
    
    M_juv_pop_storage_boot_clean <-
      data.frame(M_juv_pop_storage_boot) %>% 
      dplyr::rename(pop_size = M_juv_pop_storage_boot) %>% 
      mutate(age = c(0:(n-1))) %>% 
      mutate(sex = "M",
             species = species)
    
    # calc JSR as the proportion of the juvenile population that is male
    out <- 
      left_join(F_juv_pop_storage_boot_clean, 
                M_juv_pop_storage_boot_clean, 
                by = c("species", "age"))  %>% 
      mutate(JSR = pop_size.y / (pop_size.y + pop_size.x))
    
    out
    
  }
