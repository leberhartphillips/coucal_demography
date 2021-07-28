# load libraries
source("R/project/project_libraries.R")
source("R/project/project_plotting.R")

# load functions
function.sources = list.files(path = "R/functions", 
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
sapply(function.sources, source, .GlobalEnv)

# load capture histories
data.sources = list.files(path = "data/cooked", 
                          pattern="*ch.rds$", full.names = TRUE, 
                          ignore.case = TRUE)
sapply(data.sources, load, .GlobalEnv)  

# load parameter distributions
source("R/analysis/analysis_egg_survival.R")
source("R/analysis/analysis_clutch_size.R")
source("R/analysis/analysis_hatching_sex_ratio.R")
source("R/analysis/analysis_immigration_sex_ratio.R")
source("R/analysis/analysis_mating_system.R")
source("R/analysis/analysis_age_first_flight.R")
source("R/analysis/analysis_fledge_age.R")

# consolidate
parameter_distributions <- 
  bind_rows(coucal_egg_survival,
            coucal_clutch_size,
            coucal_HSR,
            coucal_ISR,
            coucal_mating_system,
            coucal_fledge_age,
            coucal_flight_age) %>% 
  dplyr::select(trait, species, sex, mean, sd)

# load output of JSR stochastic bootstrap
BC_JSR_out <- 
  readRDS("output/bootstraps/hazard/cooked/BC_hazard_JSR_bootstrap_result_stoch_5050_HSR.rds")

WBC_JSR_out <- 
  readRDS("output/bootstraps/hazard/cooked/WBC_hazard_JSR_bootstrap_result_stoch_5050_HSR.rds")

CI <- 0.95

# summarise JSR by species and age (i.e., days since hatch)
JSR_boot <- 
  bind_rows(BC_JSR_out,
            WBC_JSR_out) %>%
  mutate(species = factor(species, levels = c("BC", "WBC"))) %>% 
  dplyr::group_by(species, age) %>%
  dplyr::summarise(ucl_JSR = stats::quantile(JSR, (1 - CI)/2, na.rm = TRUE),
                   lcl_JSR = stats::quantile(JSR, 1 - (1 - CI)/2, na.rm = TRUE),
                   avg_JSR = mean(JSR, na.rm = TRUE),
                   med_JSR = median(JSR, na.rm = TRUE),
                   max_JSR = max(JSR, na.rm = TRUE),
                   min_JSR = min(JSR, na.rm = TRUE))

# clean up the output from the bootstrap procedure
BC_hazard_rate_boot_tidy <- 
  hazard_boot_out_wrangle(species = "BC", niter = 1000, 
                          output_dir = "output/bootstraps/hazard/cooked/",
                          rds_file = "_hazard_ASR_bootstrap_result_w_WBC_ad_surv_stoc")

WBC_hazard_rate_boot_tidy <- 
  hazard_boot_out_wrangle(species = "WBC", niter = 1000, 
                          output_dir = "output/bootstraps/hazard/cooked/",
                          rds_file = "_hazard_ASR_bootstrap_result_w_WBC_ad_surv_stoc")

hazard_rate_boot_tidy <- 
  bind_rows(BC_hazard_rate_boot_tidy, WBC_hazard_rate_boot_tidy)

# extract age-specific survival rates from the bootstrap procedure and summarize by age
survival_rates <- 
  bind_rows(BC_hazard_rate_boot_tidy$hazard_rates_boot,
            WBC_hazard_rate_boot_tidy$hazard_rates_boot) %>% 
  Rmisc::summarySE(.,
                   measurevar = "fit",
                   groupvars = c("species", "sex", "age"),
                   conf.interval = 0.95) %>% 
  mutate(fit = 1 - fit)

# # number of stages in the survival process
# no_ages <- max(JSR_boot$age)
# 
# # Define plover life-stages of the Ceuta snowy plover matrix model
# ages <- 
#   data.frame(text = "age",
#              number = seq(1:no_ages)) %>% 
#   mutate(concat = paste(text, number, sep = "_")) %>% 
#   pull(concat)

#### Arguments ----
species_name = "BC"
imm_N_base = 100
ISR_rate = pull(filter(parameter_distributions, 
                        species == species_name & trait == "immigration_sex_ratio"), 
                 mean)
k_rate = pull(filter(parameter_distributions, 
                        species == species_name & trait == "clutch_size"), 
                 mean)
h_rate = pull(filter(parameter_distributions, 
                        species == species_name & trait == "mating_system"), 
                 mean)
HSR_rate = pull(filter(parameter_distributions, 
                      species == species_name & trait == "hatching_sex_ratio"), 
               mean)
egg_S_rate = pull(filter(parameter_distributions, 
                         species == species_name & trait == "egg_survival"), 
                  mean)
JSR_base = 
  bind_rows(BC_JSR_out,
            WBC_JSR_out) %>%
  mutate(species = factor(species, levels = c("BC", "WBC"))) %>% 
  dplyr::group_by(species, age) %>%
  dplyr::summarise(ucl_JSR = stats::quantile(JSR, (1 - CI)/2, na.rm = TRUE),
                   lcl_JSR = stats::quantile(JSR, 1 - (1 - CI)/2, na.rm = TRUE),
                   avg_JSR = mean(JSR, na.rm = TRUE),
                   med_JSR = median(JSR, na.rm = TRUE),
                   max_JSR = max(JSR, na.rm = TRUE),
                   min_JSR = min(JSR, na.rm = TRUE)) %>% 
  filter(age == (length(unique(survival_rates$age)) - 2) & species == species_name) %>% 
  pull(avg_JSR) 
#### -----

# dataframe to store the perturbation results
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

flight_dat <- 
  data.frame(species = c("BC","WBC"),
             end_nestling = c(13, 14),
             end_nestling_lower = c(12, 13),
             end_nestling_upper = c(13, 15),
             end_groundling = c(36, 32),
             end_groundling_lower = c(34, 29),
             end_groundling_upper = c(38, 35))

JSR_pert_results %>%
  filter(str_detect(parameter, "_age_")) %>% 
  mutate(sex = ifelse(str_detect(parameter, "F_"), "Female", "Male"),
         age = c(0:69),
         abs_sen = abs(sensitivities)) %>% 
  filter(age != 69) %>% 
  ggplot(data = .) +
  geom_line(aes(x = age, y = abs_sen, color = sex)) +
  geom_vline(data = filter(flight_dat, species == "BC"),
             aes(xintercept = end_nestling), 
             linetype = "dashed", alpha = 0.5, color = "grey20") +
  geom_vline(data = filter(flight_dat, species == "BC"),
             aes(xintercept = end_groundling),
             linetype = "dashed", alpha = 0.5, color = "grey20") +
  ylab("Sensitivity of JSR to changes in daily survival") +
  xlab("Age (Days since hatching)") +
  annotate(geom = "text", x = 15/2, 
           y = filter(JSR_pert_results, str_detect(parameter, negate = TRUE, pattern = "_69") & str_detect(parameter, pattern = "_age_")) %>% 
                pull(sensitivities) %>% 
                abs() %>% 
                min() %>% 
                round(3),
           label = "nestling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
  annotate(geom = "text", 
           y = filter(JSR_pert_results, str_detect(parameter, negate = TRUE, pattern = "_69") & str_detect(parameter, pattern = "_age_")) %>% 
             pull(sensitivities) %>% 
             abs() %>%  
             min() %>% 
             round(3),
           x = (36 - 15)/2 + 15,
           label = "groundling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
  annotate(geom = "text", 
           y = filter(JSR_pert_results, str_detect(parameter, negate = TRUE, pattern = "_69") & str_detect(parameter, pattern = "_age_")) %>% 
             pull(sensitivities) %>% 
             abs() %>% 
             min() %>% 
             round(3),
           x = (70 - 36)/2 + 36,
           label = "fledgling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
  scale_x_continuous(limits = c(0, 70), expand = c(0, 0))

JSR_pert_results %>%
  filter(str_detect(parameter, "_age_")) %>% 
  mutate(sex = ifelse(str_detect(parameter, "F_"), "Female", "Male"),
         age = c(0:69)) %>% 
  filter(age != 69) %>% 
  ggplot(data = .) +
  geom_line(aes(x = age, y = elasticities, color = sex)) +
  geom_vline(data = filter(flight_dat, species == "BC"),
             aes(xintercept = end_nestling), 
             linetype = "dashed", alpha = 0.5, color = "grey20") +
  geom_vline(data = filter(flight_dat, species == "BC"),
             aes(xintercept = end_groundling),
             linetype = "dashed", alpha = 0.5, color = "grey20") +
  ylab("Elasticity of JSR to changes in daily survival") +
  xlab("Age (Days since hatching)") #+
  annotate(geom = "text", x = 15/2, 
           y = filter(JSR_pert_results, str_detect(parameter, negate = TRUE, pattern = "_69") & str_detect(parameter, pattern = "_age_")) %>% 
             pull(elasticities) %>% 
             abs() %>% 
             min() %>% 
             round(3),
           label = "nestling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
  annotate(geom = "text", 
           y = filter(JSR_pert_results, str_detect(parameter, negate = TRUE, pattern = "_69") & str_detect(parameter, pattern = "_age_")) %>% 
             pull(elasticities) %>% 
             abs() %>% 
             min() %>% 
             round(3),
           x = (36 - 15)/2 + 15,
           label = "groundling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
  annotate(geom = "text", 
           y = filter(JSR_pert_results, str_detect(parameter, negate = TRUE, pattern = "_69") & str_detect(parameter, pattern = "_age_")) %>% 
             pull(elasticities) %>% 
             abs() %>% 
             min() %>% 
             round(3),
           x = (70 - 36)/2 + 36,
           label = "fledgling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
  scale_x_continuous(limits = c(0, 70), expand = c(0, 0))
