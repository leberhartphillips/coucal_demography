# sensitivity_analysis() takes the vital rate summary of the bootstrap procedure
# and conducts a perturbation analysis on each rate to assess how a proportional
# change in a given vital rate changes the ASR

sensitivity_analysis_JSR <-
  function(VR_dataframe, JSR_base){
    
    surv_rates <- 
      filter(VR_dataframe, str_detect(parameter, "_age_")) %>% 
      mutate(sex = ifelse(str_detect(parameter, "Female"), "Female", "Male"),
             age = readr::parse_number(parameter))
    
    ISR_rate <- filter(VR_dataframe, parameter == "ISR") %>% pull(obs)
    
    h_rate <- filter(VR_dataframe, parameter == "h") %>% pull(obs)
    
    HSR_rate <- filter(VR_dataframe, parameter == "HSR") %>% pull(obs)
    
    k_rate <- filter(VR_dataframe, parameter == "k") %>% pull(obs)
    
    egg_survival_rate  <- filter(VR_dataframe, parameter == "egg_S") %>% pull(obs)
    
    imm_N_base <- filter(VR_dataframe, parameter == "imm_N") %>% pull(obs)
    
    JSR_pert_results <-
      data.frame(parameter = c(unique(paste0("M_age_", survival_rates$age)),
                               unique(paste0("F_age_", survival_rates$age)),
                               "h", "k", "HSR", "ISR", "egg_S", "imm_N"),
                 sensitivities = numeric((length(unique(survival_rates$age)) * 2) + 6),
                 elasticities = numeric((length(unique(survival_rates$age)) * 2) + 6))
    
    # specify how many sex-specific survival rates there are 
    n <- length(unique(survival_rates$age)) * 2
    
    # create vectors of perturbations to test on parameters of the survival process
    vr_nums <- seq(0, 1, 0.01) # proportional changes in survival, HSR, ISR, egg_S (i.e., between 0 an 1)
    h_nums <- seq(0, 2, 0.02) # proportional changes in h index (i.e., between 0 and 2)
    k_nums <- seq(3, 5, 0.02) # proportional changes in k (i.e, between 3 and 5)
    imm_N_nums <- seq(10, 1010, 10) # proportional changes in imm_N (i.e, between 10 and 1010)
    
    # create empty dataframes to store the perturbation results for JSR
    vr_pert_JSR <- matrix(numeric(n * length(vr_nums)),
                          ncol = n, dimnames = list(vr_nums, 
                                                    JSR_pert_results$parameter[1:(nrow(JSR_pert_results)-6)]))
    
    h_pert_JSR <- matrix(numeric(length(h_nums)),
                         ncol = 1, dimnames = list(h_nums, "h"))
    
    k_pert_JSR <- matrix(numeric(length(k_nums)),
                         ncol = 1, dimnames = list(k_nums, "k"))
    
    HSR_pert_JSR <- matrix(numeric(length(vr_nums)),
                           ncol = 1, dimnames = list(vr_nums, "HSR"))
    
    ISR_pert_JSR <- matrix(numeric(length(vr_nums)),
                           ncol = 1, dimnames = list(vr_nums, "ISR"))
    
    egg_S_pert_JSR <- matrix(numeric(length(vr_nums)),
                             ncol = 1, dimnames = list(vr_nums, "egg_S"))
    
    imm_N_pert_JSR <- matrix(numeric(length(imm_N_nums)),
                             ncol = 1, dimnames = list(imm_N_nums, "imm_N"))
    
    vr <- 
      filter(survival_rates, species == species_name) %>% 
      dplyr::select(species, sex, age, fit)
    
    ##### perturbation of survival rates ####
    for (g in 1:n){ # pick a column (i.e., a sex-specific age)
      
      vr2 <- vr # reset the vital rates to the original
      
      F_juv_pop_storage_boot <- 
        matrix(, nrow = n/2, ncol = 1)
      
      M_juv_pop_storage_boot <- 
        matrix(, nrow = n/2, ncol = 1)
      
      N_F_immigrants <- 
        imm_N_base * (1 - ISR_rate)
      
      N_M_immigrants <- 
        imm_N_base * ISR_rate
      
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
      
      F_juv_pop_storage_boot[1, 1] <- N_F_eggs * egg_S_rate
      M_juv_pop_storage_boot[1, 1] <- N_M_eggs * egg_S_rate
      
      for (i in 1:length(vr_nums)){ # pick a perturbation level
        
        vr2[g, "fit"] <- vr_nums[i] # specify the vital rate with the new perturbation level
        
        for (j in 2:(n/2)) { # push the population through the survival process
          
          F_juv_pop_storage_boot[j, 1] <- 
            F_juv_pop_storage_boot[j - 1, 1] * vr2[which(vr2$sex == "Female"), ][["fit"]][j - 1]
          M_juv_pop_storage_boot[j, 1] <- 
            M_juv_pop_storage_boot[j - 1, 1] * vr2[which(vr2$sex == "Male"), ][["fit"]][j - 1]
        }
        
        F_juv_pop_storage_boot_clean <-
          data.frame(F_juv_pop_storage_boot) %>% 
          dplyr::rename(pop_size = F_juv_pop_storage_boot) %>% 
          mutate(age = c(0:((n/2)-1))) %>% 
          mutate(sex = "F", 
                 species = species_name)
        
        M_juv_pop_storage_boot_clean <-
          data.frame(M_juv_pop_storage_boot) %>% 
          dplyr::rename(pop_size = M_juv_pop_storage_boot) %>% 
          mutate(age = c(0:((n/2)-1))) %>% 
          mutate(sex = "M", 
                 species = species_name)
        
        # calc JSR as the proportion of the juvenile population that is male
        out <- 
          left_join(F_juv_pop_storage_boot_clean, 
                    M_juv_pop_storage_boot_clean, 
                    by = c("species", "age"))  %>% 
          mutate(JSR = pop_size.y / (pop_size.y + pop_size.x))
        
        # extract and store the final JSR value
        vr_pert_JSR[i, g] <- out[(n/2), "JSR"]
      }
      
      # get the spline function of JSR
      spl_JSR <- 
        smooth.spline(vr_pert_JSR[which(!is.na(vr_pert_JSR[, g])), g] ~ 
                        names(vr_pert_JSR[which(!is.na(vr_pert_JSR[, g])), g]))
      
      # estimate the slope of the tangent of the spline at the vital rate
      JSR_pert_results[g, 2] <- predict(spl_JSR, x = vr[g, "fit"], deriv = 1)$y
      
      # re-scale sensitivity into elasticity
      JSR_pert_results[g, 3] <- vr[g, "fit"] / JSR_base * JSR_pert_results[g, 2]
      
    }
    
    ##### perturbation of egg survival rate ####
    for (g in 1:1){ 
      
      for (i in 1:length(vr_nums)){ # pick a perturbation level
        
        vr2 <- vr # reset the vital rates to the original
        
        F_juv_pop_storage_boot <- 
          matrix(, nrow = n/2, ncol = 1)
        
        M_juv_pop_storage_boot <- 
          matrix(, nrow = n/2, ncol = 1)
        
        N_F_immigrants <- 
          imm_N_base * (1 - ISR_rate)
        
        N_M_immigrants <- 
          imm_N_base * ISR_rate
        
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
        
        # specify the egg survival with the new perturbation level
        F_juv_pop_storage_boot[1, 1] <- N_F_eggs * vr_nums[i]
        M_juv_pop_storage_boot[1, 1] <- N_M_eggs * vr_nums[i]
        
        # vr2[g, "fit"] <- vr_nums[i] # specify the vital rate with the new perturbation level
        
        for (j in 2:(n/2)) { # push the population through the survival process
          
          F_juv_pop_storage_boot[j, 1] <- 
            F_juv_pop_storage_boot[j - 1, 1] * vr2[which(vr2$sex == "Female"), ][["fit"]][j - 1]
          M_juv_pop_storage_boot[j, 1] <- 
            M_juv_pop_storage_boot[j - 1, 1] * vr2[which(vr2$sex == "Male"), ][["fit"]][j - 1]
        }
        
        F_juv_pop_storage_boot_clean <-
          data.frame(F_juv_pop_storage_boot) %>% 
          dplyr::rename(pop_size = F_juv_pop_storage_boot) %>% 
          mutate(age = c(0:((n/2)-1))) %>% 
          mutate(sex = "F", 
                 species = species_name)
        
        M_juv_pop_storage_boot_clean <-
          data.frame(M_juv_pop_storage_boot) %>% 
          dplyr::rename(pop_size = M_juv_pop_storage_boot) %>% 
          mutate(age = c(0:((n/2)-1))) %>% 
          mutate(sex = "M", 
                 species = species_name)
        
        # calc JSR as the proportion of the juvenile population that is male
        out <- 
          left_join(F_juv_pop_storage_boot_clean, 
                    M_juv_pop_storage_boot_clean, 
                    by = c("species", "age"))  %>% 
          mutate(JSR = pop_size.y / (pop_size.y + pop_size.x))
        
        # extract and store the final JSR value
        egg_S_pert_JSR[i, g] <- out[(n/2), "JSR"]
      }
      
      # get the spline function of JSR
      spl_JSR <- 
        smooth.spline(egg_S_pert_JSR[which(!is.na(egg_S_pert_JSR[, g])), g] ~ 
                        names(egg_S_pert_JSR[which(!is.na(egg_S_pert_JSR[, g])), g]))
      
      # estimate the slope of the tangent of the spline at the vital rate
      JSR_pert_results[which(JSR_pert_results$parameter == "egg_S"), 2] <- 
        predict(spl_JSR, x = egg_S_rate, deriv = 1)$y
      
      # re-scale sensitivity into elasticity
      JSR_pert_results[which(JSR_pert_results$parameter == "egg_S"), 3] <- 
        egg_S_rate / JSR_base * JSR_pert_results[which(JSR_pert_results$parameter == "egg_S"), 2]
      
    }
    
    ##### perturbation of hatching sex ratio ####
    for (g in 1:1){ 
      
      for (i in 1:length(vr_nums)){ # pick a perturbation level
        
        vr2 <- vr # reset the vital rates to the original
        
        F_juv_pop_storage_boot <- 
          matrix(, nrow = n/2, ncol = 1)
        
        M_juv_pop_storage_boot <- 
          matrix(, nrow = n/2, ncol = 1)
        
        N_F_immigrants <- 
          imm_N_base * (1 - ISR_rate)
        
        N_M_immigrants <- 
          imm_N_base * ISR_rate
        
        # Female freq-dep fecundity of Female chicks
        F_F_eggs <- 
          ((k_rate * N_M_immigrants) / (N_M_immigrants + (N_F_immigrants / h_rate)) * (1 - vr_nums[i]))
        
        # Female freq-dep fecundity of Male chicks
        F_M_eggs <- 
          ((k_rate * N_M_immigrants) / (N_M_immigrants + (N_F_immigrants / h_rate)) * vr_nums[i])
        
        # Male freq-dep fecundity of Female chicks
        M_F_eggs <- 
          ((k_rate * N_F_immigrants) / (N_M_immigrants + (N_F_immigrants / h_rate)) * (1 - vr_nums[i]))
        
        # Male freq-dep fecundity of Male chicks
        M_M_eggs <- 
          ((k_rate * N_F_immigrants) / (N_M_immigrants + (N_F_immigrants / h_rate)) * vr_nums[i])
        
        N_F_eggs <- F_F_eggs * N_F_immigrants + M_F_eggs * N_M_immigrants
        N_M_eggs <- F_M_eggs * N_F_immigrants + M_M_eggs * N_M_immigrants
        
        # specify the egg survival with the new perturbation level
        F_juv_pop_storage_boot[1, 1] <- N_F_eggs * egg_S_rate
        M_juv_pop_storage_boot[1, 1] <- N_M_eggs * egg_S_rate
        
        for (j in 2:(n/2)) { # push the population through the survival process
          
          F_juv_pop_storage_boot[j, 1] <- 
            F_juv_pop_storage_boot[j - 1, 1] * vr2[which(vr2$sex == "Female"), ][["fit"]][j - 1]
          M_juv_pop_storage_boot[j, 1] <- 
            M_juv_pop_storage_boot[j - 1, 1] * vr2[which(vr2$sex == "Male"), ][["fit"]][j - 1]
        }
        
        F_juv_pop_storage_boot_clean <-
          data.frame(F_juv_pop_storage_boot) %>% 
          dplyr::rename(pop_size = F_juv_pop_storage_boot) %>% 
          mutate(age = c(0:((n/2)-1))) %>% 
          mutate(sex = "F", 
                 species = species_name)
        
        M_juv_pop_storage_boot_clean <-
          data.frame(M_juv_pop_storage_boot) %>% 
          dplyr::rename(pop_size = M_juv_pop_storage_boot) %>% 
          mutate(age = c(0:((n/2)-1))) %>% 
          mutate(sex = "M", 
                 species = species_name)
        
        # calc JSR as the proportion of the juvenile population that is male
        out <- 
          left_join(F_juv_pop_storage_boot_clean, 
                    M_juv_pop_storage_boot_clean, 
                    by = c("species", "age"))  %>% 
          mutate(JSR = pop_size.y / (pop_size.y + pop_size.x))
        
        # extract and store the final JSR value
        HSR_pert_JSR[i, g] <- out[(n/2), "JSR"]
      }
      
      # get the spline function of JSR
      spl_JSR <- 
        smooth.spline(HSR_pert_JSR[which(!is.na(HSR_pert_JSR[, g])), g] ~ 
                        names(HSR_pert_JSR[which(!is.na(HSR_pert_JSR[, g])), g]))
      
      # estimate the slope of the tangent of the spline at the vital rate
      JSR_pert_results[which(JSR_pert_results$parameter == "HSR"), 2] <- 
        predict(spl_JSR, x = HSR_rate, deriv = 1)$y
      
      # re-scale sensitivity into elasticity
      JSR_pert_results[which(JSR_pert_results$parameter == "HSR"), 3] <- 
        HSR_rate / JSR_base * JSR_pert_results[which(JSR_pert_results$parameter == "HSR"), 2]
      
    }
    
    ##### perturbation of mating system ####
    for (g in 1:1){ 
      
      for (i in 1:length(vr_nums)){ # pick a perturbation level
        
        vr2 <- vr # reset the vital rates to the original
        
        F_juv_pop_storage_boot <- 
          matrix(, nrow = n/2, ncol = 1)
        
        M_juv_pop_storage_boot <- 
          matrix(, nrow = n/2, ncol = 1)
        
        N_F_immigrants <- 
          imm_N_base * (1 - ISR_rate)
        
        N_M_immigrants <- 
          imm_N_base * ISR_rate
        
        # Female freq-dep fecundity of Female chicks
        F_F_eggs <- 
          ((k_rate * N_M_immigrants) / (N_M_immigrants + (N_F_immigrants / vr_nums[i])) * (1 - HSR_rate))
        
        # Female freq-dep fecundity of Male chicks
        F_M_eggs <- 
          ((k_rate * N_M_immigrants) / (N_M_immigrants + (N_F_immigrants / vr_nums[i])) * HSR_rate)
        
        # Male freq-dep fecundity of Female chicks
        M_F_eggs <- 
          ((k_rate * N_F_immigrants) / (N_M_immigrants + (N_F_immigrants / vr_nums[i])) * (1 - HSR_rate))
        
        # Male freq-dep fecundity of Male chicks
        M_M_eggs <- 
          ((k_rate * N_F_immigrants) / (N_M_immigrants + (N_F_immigrants / vr_nums[i])) * HSR_rate)
        
        N_F_eggs <- F_F_eggs * N_F_immigrants + M_F_eggs * N_M_immigrants
        N_M_eggs <- F_M_eggs * N_F_immigrants + M_M_eggs * N_M_immigrants
        
        # specify the egg survival with the new perturbation level
        F_juv_pop_storage_boot[1, 1] <- N_F_eggs * egg_S_rate
        M_juv_pop_storage_boot[1, 1] <- N_M_eggs * egg_S_rate
        
        for (j in 2:(n/2)) { # push the population through the survival process
          
          F_juv_pop_storage_boot[j, 1] <- 
            F_juv_pop_storage_boot[j - 1, 1] * vr2[which(vr2$sex == "Female"), ][["fit"]][j - 1]
          M_juv_pop_storage_boot[j, 1] <- 
            M_juv_pop_storage_boot[j - 1, 1] * vr2[which(vr2$sex == "Male"), ][["fit"]][j - 1]
        }
        
        F_juv_pop_storage_boot_clean <-
          data.frame(F_juv_pop_storage_boot) %>% 
          dplyr::rename(pop_size = F_juv_pop_storage_boot) %>% 
          mutate(age = c(0:((n/2)-1))) %>% 
          mutate(sex = "F", 
                 species = species_name)
        
        M_juv_pop_storage_boot_clean <-
          data.frame(M_juv_pop_storage_boot) %>% 
          dplyr::rename(pop_size = M_juv_pop_storage_boot) %>% 
          mutate(age = c(0:((n/2)-1))) %>% 
          mutate(sex = "M", 
                 species = species_name)
        
        # calc JSR as the proportion of the juvenile population that is male
        out <- 
          left_join(F_juv_pop_storage_boot_clean, 
                    M_juv_pop_storage_boot_clean, 
                    by = c("species", "age"))  %>% 
          mutate(JSR = pop_size.y / (pop_size.y + pop_size.x))
        
        # extract and store the final JSR value
        h_pert_JSR[i, g] <- out[(n/2), "JSR"]
      }
      
      # get the spline function of JSR
      spl_JSR <- 
        smooth.spline(h_pert_JSR[which(!is.na(h_pert_JSR[, g])), g] ~ 
                        names(h_pert_JSR[which(!is.na(h_pert_JSR[, g])), g]))
      
      # estimate the slope of the tangent of the spline at the vital rate
      JSR_pert_results[which(JSR_pert_results$parameter == "h"), 2] <- 
        predict(spl_JSR, x = h_rate, deriv = 1)$y
      
      # re-scale sensitivity into elasticity
      JSR_pert_results[which(JSR_pert_results$parameter == "h"), 3] <- 
        h_rate / JSR_base * JSR_pert_results[which(JSR_pert_results$parameter == "h"), 2]
      
    }
    
    ##### perturbation of clutch size ####
    for (g in 1:1){ 
      
      for (i in 1:length(vr_nums)){ # pick a perturbation level
        
        vr2 <- vr # reset the vital rates to the original
        
        F_juv_pop_storage_boot <- 
          matrix(, nrow = n/2, ncol = 1)
        
        M_juv_pop_storage_boot <- 
          matrix(, nrow = n/2, ncol = 1)
        
        N_F_immigrants <- 
          imm_N_base * (1 - ISR_rate)
        
        N_M_immigrants <- 
          imm_N_base * ISR_rate
        
        # Female freq-dep fecundity of Female chicks
        F_F_eggs <- 
          ((vr_nums[i] * N_M_immigrants) / (N_M_immigrants + (N_F_immigrants / h_rate)) * (1 - HSR_rate))
        
        # Female freq-dep fecundity of Male chicks
        F_M_eggs <- 
          ((vr_nums[i] * N_M_immigrants) / (N_M_immigrants + (N_F_immigrants / h_rate)) * HSR_rate)
        
        # Male freq-dep fecundity of Female chicks
        M_F_eggs <- 
          ((vr_nums[i] * N_F_immigrants) / (N_M_immigrants + (N_F_immigrants / h_rate)) * (1 - HSR_rate))
        
        # Male freq-dep fecundity of Male chicks
        M_M_eggs <- 
          ((vr_nums[i] * N_F_immigrants) / (N_M_immigrants + (N_F_immigrants / h_rate)) * HSR_rate)
        
        N_F_eggs <- F_F_eggs * N_F_immigrants + M_F_eggs * N_M_immigrants
        N_M_eggs <- F_M_eggs * N_F_immigrants + M_M_eggs * N_M_immigrants
        
        # specify the egg survival with the new perturbation level
        F_juv_pop_storage_boot[1, 1] <- N_F_eggs * egg_S_rate
        M_juv_pop_storage_boot[1, 1] <- N_M_eggs * egg_S_rate
        
        for (j in 2:(n/2)) { # push the population through the survival process
          
          F_juv_pop_storage_boot[j, 1] <- 
            F_juv_pop_storage_boot[j - 1, 1] * vr2[which(vr2$sex == "Female"), ][["fit"]][j - 1]
          M_juv_pop_storage_boot[j, 1] <- 
            M_juv_pop_storage_boot[j - 1, 1] * vr2[which(vr2$sex == "Male"), ][["fit"]][j - 1]
        }
        
        F_juv_pop_storage_boot_clean <-
          data.frame(F_juv_pop_storage_boot) %>% 
          dplyr::rename(pop_size = F_juv_pop_storage_boot) %>% 
          mutate(age = c(0:((n/2)-1))) %>% 
          mutate(sex = "F", 
                 species = species_name)
        
        M_juv_pop_storage_boot_clean <-
          data.frame(M_juv_pop_storage_boot) %>% 
          dplyr::rename(pop_size = M_juv_pop_storage_boot) %>% 
          mutate(age = c(0:((n/2)-1))) %>% 
          mutate(sex = "M", 
                 species = species_name)
        
        # calc JSR as the proportion of the juvenile population that is male
        out <- 
          left_join(F_juv_pop_storage_boot_clean, 
                    M_juv_pop_storage_boot_clean, 
                    by = c("species", "age"))  %>% 
          mutate(JSR = pop_size.y / (pop_size.y + pop_size.x))
        
        # extract and store the final JSR value
        k_pert_JSR[i, g] <- out[(n/2), "JSR"]
      }
      
      # get the spline function of JSR
      spl_JSR <- 
        smooth.spline(k_pert_JSR[which(!is.na(k_pert_JSR[, g])), g] ~ 
                        names(k_pert_JSR[which(!is.na(k_pert_JSR[, g])), g]))
      
      # estimate the slope of the tangent of the spline at the vital rate
      JSR_pert_results[which(JSR_pert_results$parameter == "k"), 2] <- 
        predict(spl_JSR, x = k_rate, deriv = 1)$y
      
      # re-scale sensitivity into elasticity
      JSR_pert_results[which(JSR_pert_results$parameter == "k"), 3] <- 
        k_rate / JSR_base * JSR_pert_results[which(JSR_pert_results$parameter == "k"), 2]
      
    }
    
    ##### perturbation of immigration sex ratio ####
    for (g in 1:1){ 
      
      for (i in 1:length(vr_nums)){ # pick a perturbation level
        
        vr2 <- vr # reset the vital rates to the original
        
        F_juv_pop_storage_boot <- 
          matrix(, nrow = n/2, ncol = 1)
        
        M_juv_pop_storage_boot <- 
          matrix(, nrow = n/2, ncol = 1)
        
        N_F_immigrants <- 
          imm_N_base * (1 - vr_nums[i])
        
        N_M_immigrants <- 
          imm_N_base * vr_nums[i]
        
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
        
        # specify the egg survival with the new perturbation level
        F_juv_pop_storage_boot[1, 1] <- N_F_eggs * egg_S_rate
        M_juv_pop_storage_boot[1, 1] <- N_M_eggs * egg_S_rate
        
        for (j in 2:(n/2)) { # push the population through the survival process
          
          F_juv_pop_storage_boot[j, 1] <- 
            F_juv_pop_storage_boot[j - 1, 1] * vr2[which(vr2$sex == "Female"), ][["fit"]][j - 1]
          M_juv_pop_storage_boot[j, 1] <- 
            M_juv_pop_storage_boot[j - 1, 1] * vr2[which(vr2$sex == "Male"), ][["fit"]][j - 1]
        }
        
        F_juv_pop_storage_boot_clean <-
          data.frame(F_juv_pop_storage_boot) %>% 
          dplyr::rename(pop_size = F_juv_pop_storage_boot) %>% 
          mutate(age = c(0:((n/2)-1))) %>% 
          mutate(sex = "F", 
                 species = species_name)
        
        M_juv_pop_storage_boot_clean <-
          data.frame(M_juv_pop_storage_boot) %>% 
          dplyr::rename(pop_size = M_juv_pop_storage_boot) %>% 
          mutate(age = c(0:((n/2)-1))) %>% 
          mutate(sex = "M", 
                 species = species_name)
        
        # calc JSR as the proportion of the juvenile population that is male
        out <- 
          left_join(F_juv_pop_storage_boot_clean, 
                    M_juv_pop_storage_boot_clean, 
                    by = c("species", "age"))  %>% 
          mutate(JSR = pop_size.y / (pop_size.y + pop_size.x))
        
        # extract and store the final JSR value
        ISR_pert_JSR[i, g] <- out[(n/2), "JSR"]
      }
      
      # get the spline function of JSR
      spl_JSR <- 
        smooth.spline(ISR_pert_JSR[which(!is.na(ISR_pert_JSR[, g])), g] ~ 
                        names(ISR_pert_JSR[which(!is.na(ISR_pert_JSR[, g])), g]))
      
      # estimate the slope of the tangent of the spline at the vital rate
      JSR_pert_results[which(JSR_pert_results$parameter == "ISR"), 2] <- 
        predict(spl_JSR, x = ISR_rate, deriv = 1)$y
      
      # re-scale sensitivity into elasticity
      JSR_pert_results[which(JSR_pert_results$parameter == "ISR"), 3] <- 
        ISR_rate / JSR_base * JSR_pert_results[which(JSR_pert_results$parameter == "ISR"), 2]
      
    }
    
    ##### perturbation of immigration population size ####
    for (g in 1:1){ 
      
      for (i in 1:length(vr_nums)){ # pick a perturbation level
        
        vr2 <- vr # reset the vital rates to the original
        
        F_juv_pop_storage_boot <- 
          matrix(, nrow = n/2, ncol = 1)
        
        M_juv_pop_storage_boot <- 
          matrix(, nrow = n/2, ncol = 1)
        
        N_F_immigrants <- 
          vr_nums[i] * (1 - ISR_rate)
        
        N_M_immigrants <- 
          vr_nums[i] * ISR_rate
        
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
        
        # specify the egg survival with the new perturbation level
        F_juv_pop_storage_boot[1, 1] <- N_F_eggs * egg_S_rate
        M_juv_pop_storage_boot[1, 1] <- N_M_eggs * egg_S_rate
        
        for (j in 2:(n/2)) { # push the population through the survival process
          
          F_juv_pop_storage_boot[j, 1] <- 
            F_juv_pop_storage_boot[j - 1, 1] * vr2[which(vr2$sex == "Female"), ][["fit"]][j - 1]
          M_juv_pop_storage_boot[j, 1] <- 
            M_juv_pop_storage_boot[j - 1, 1] * vr2[which(vr2$sex == "Male"), ][["fit"]][j - 1]
        }
        
        F_juv_pop_storage_boot_clean <-
          data.frame(F_juv_pop_storage_boot) %>% 
          dplyr::rename(pop_size = F_juv_pop_storage_boot) %>% 
          mutate(age = c(0:((n/2)-1))) %>% 
          mutate(sex = "F", 
                 species = species_name)
        
        M_juv_pop_storage_boot_clean <-
          data.frame(M_juv_pop_storage_boot) %>% 
          dplyr::rename(pop_size = M_juv_pop_storage_boot) %>% 
          mutate(age = c(0:((n/2)-1))) %>% 
          mutate(sex = "M", 
                 species = species_name)
        
        # calc JSR as the proportion of the juvenile population that is male
        out <- 
          left_join(F_juv_pop_storage_boot_clean, 
                    M_juv_pop_storage_boot_clean, 
                    by = c("species", "age"))  %>% 
          mutate(JSR = pop_size.y / (pop_size.y + pop_size.x))
        
        # extract and store the final JSR value
        imm_N_pert_JSR[i, g] <- out[(n/2), "JSR"]
      }
      
      # get the spline function of JSR
      spl_JSR <- 
        smooth.spline(imm_N_pert_JSR[which(!is.na(imm_N_pert_JSR[, g])), g] ~ 
                        names(imm_N_pert_JSR[which(!is.na(imm_N_pert_JSR[, g])), g]))
      
      # estimate the slope of the tangent of the spline at the vital rate
      JSR_pert_results[which(JSR_pert_results$parameter == "imm_N"), 2] <- 
        predict(spl_JSR, x = imm_N_base, deriv = 1)$y
      
      # re-scale sensitivity into elasticity
      JSR_pert_results[which(JSR_pert_results$parameter == "imm_N"), 3] <- 
        imm_N_base / JSR_base * JSR_pert_results[which(JSR_pert_results$parameter == "imm_N"), 2]
      
    }
  
    #### store all results into a list ----
    result <- list(JSR_pert_results = JSR_pert_results,
                   vr_pert_JSR = vr_pert_JSR,
                   egg_S_pert_JSR = egg_S_pert_JSR,
                   HSR_pert_JSR = HSR_pert_JSR,
                   h_pert_JSR = h_pert_JSR,
                   k_pert_JSR = k_pert_JSR,
                   ISR_pert_JSR = ISR_pert_JSR,
                   imm_N_pert_JSR = imm_N_pert_JSR)

  }
