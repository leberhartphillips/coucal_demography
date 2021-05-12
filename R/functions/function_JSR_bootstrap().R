source("R/project/project_libraries.R")

# JSR_bootstrap takes the bootstrapped hazard functions and simulates the
# age-specific population dynamics of juveniles after hatching while also
# taking into consideration the hatching sex ratio and mating function
JSR_bootstrap <- 
  function(niter,
           F_rates, M_rates,
           immigrant_pop_size, 
           ISR, h, HSR, k, egg_survival, species){
    
    F_juv_pop_storage_boot <- 
      matrix(, nrow = max(F_rates$age), ncol = niter)
    
    M_juv_pop_storage_boot <- 
      matrix(, nrow = max(F_rates$age), ncol = niter)
    
    N_F_immigrants <- 
      immigrant_pop_size * (1 - ISR)
    
    N_M_immigrants <- 
      immigrant_pop_size * ISR
    
    # Female freq-dep fecundity of Female chicks
    F_F_eggs <- 
      ((k * N_M_immigrants) / (N_M_immigrants + (N_F_immigrants / h)) * (1 - HSR))
    
    # Female freq-dep fecundity of Male chicks
    F_M_eggs <- 
      ((k * N_M_immigrants) / (N_M_immigrants + (N_F_immigrants / h)) * HSR)
    
    # Male freq-dep fecundity of Female chicks
    M_F_eggs <- 
      ((k * N_F_immigrants) / (N_M_immigrants + (N_F_immigrants / h)) * (1 - HSR))
    
    # Male freq-dep fecundity of Male chicks
    M_M_eggs <- 
      ((k * N_F_immigrants) / (N_M_immigrants + (N_F_immigrants / h)) * HSR)
    
    N_F_eggs <- F_F_eggs * N_F_immigrants + M_F_eggs * N_M_immigrants
    N_M_eggs <- F_M_eggs * N_F_immigrants + M_M_eggs * N_M_immigrants
    
    for(j in 1:niter){
      F_juv_pop_storage_boot[1, j] <- N_F_eggs * egg_survival
      M_juv_pop_storage_boot[1, j] <- N_M_eggs * egg_survival
      for(i in 2:max(F_rates$age)){
        F_juv_pop_storage_boot[i, j] <- 
          F_juv_pop_storage_boot[i - 1, j] * F_rates[which(F_rates$iter == j), ][["estimate"]][i - 1]
        M_juv_pop_storage_boot[i, j] <- 
          M_juv_pop_storage_boot[i - 1, j] * M_rates[which(M_rates$iter == j), ][["estimate"]][i - 1]
      }
    }
    
    F_juv_pop_storage_boot_clean <-
      data.frame(F_juv_pop_storage_boot) %>% 
      mutate(age = c(0:(max(F_rates$age)-1))) %>% 
      pivot_longer(!age, names_to = "iteration", values_to = "n") %>% 
      mutate(sex = "F", 
             species = species) %>% 
      arrange(iteration)
    
    M_juv_pop_storage_boot_clean <-
      data.frame(M_juv_pop_storage_boot) %>% 
      mutate(age = c(0:(max(F_rates$age)-1))) %>% 
      pivot_longer(!age, names_to = "iteration", values_to = "n") %>% 
      mutate(sex = "M",
             species = species) %>% 
      arrange(iteration)
    
    out <- 
      left_join(F_juv_pop_storage_boot_clean, 
                M_juv_pop_storage_boot_clean, 
                by = c("species", "age", "iteration"))  %>% 
      mutate(JSR = n.y / (n.y + n.x))
    
  }