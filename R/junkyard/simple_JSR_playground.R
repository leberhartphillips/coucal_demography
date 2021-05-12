#### one step matrix model ----

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

# clean up the output from the bootstrap procedure and save as rds
BC_hazard_rate_boot_tidy <- 
  hazard_boot_out_wrangle(species = "BC", niter = 1000, 
                          output_dir = "output/bootstraps/hazard/cooked/",
                          rds_file = "_hazard_ASR_bootstrap_result_one")

# clean up the output from the bootstrap procedure and save as rds
WBC_hazard_rate_boot_tidy <- 
  hazard_boot_out_wrangle(species = "WBC", niter = 1000, 
                          output_dir = "output/bootstraps/hazard/cooked/",
                          rds_file = "_hazard_ASR_bootstrap_result_one")

average_rates <- 
  bind_rows(BC_hazard_rate_boot_tidy$hazard_rates_boot,
            WBC_hazard_rate_boot_tidy$hazard_rates_boot) %>% 
  group_by(species, sex, age) %>% 
  dplyr::summarise(med_surv = median(estimate))

BC_F_rates <- 
  BC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  filter(sex == "Female" & species == "BC")

BC_M_rates <- 
  BC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  filter(sex == "Male" & species == "BC")

WBC_F_rates <- 
  WBC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  filter(sex == "Female" & species == "WBC")

WBC_M_rates <- 
  WBC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  filter(sex == "Male" & species == "WBC")


BC_ages <- seq(0, max(average_rates[which(average_rates$species == "BC"), "age"]), 1)
WBC_ages <- seq(0, max(average_rates[which(average_rates$species == "WBC"), "age"]), 1)

niter = 1000

BC_juv_pop_storage <- 
  matrix(numeric(2 * length(BC_ages)),
         ncol = 2, nrow = length(BC_ages), 
         dimnames = list(BC_ages, c("nF", "nM")))

WBC_juv_pop_storage <- 
  matrix(numeric(2 * length(WBC_ages)),
         ncol = 2, nrow = length(WBC_ages), 
         dimnames = list(WBC_ages, c("nF", "nM")))

BC_F_juv_pop_storage_boot <- 
  matrix(, nrow = length(BC_ages), ncol = niter)

BC_M_juv_pop_storage_boot <- 
  matrix(, nrow = length(BC_ages), ncol = niter)

WBC_F_juv_pop_storage_boot <- 
  matrix(, nrow = length(WBC_ages), ncol = niter)

WBC_M_juv_pop_storage_boot <- 
  matrix(, nrow = length(WBC_ages), ncol = niter)

# n_eggs = 100
k = 4
HSR = 0.4955
h = 1/2.9
BC_egg_survival = 0.32
adult_survival_rate = 0.3
ISR = 0.738
immigrant_pop_size = 1000
fledge_age = 15
flight_age = 36

A <- matrix(sapply(matrix_str, eval, vr, NULL), 
            nrow = sqrt(length(matrix_str)), byrow=TRUE, 
            dimnames = list(stages, stages))

N_F_immigrants <- 
  immigrant_pop_size * (1 - ISR)

N_M_immigrants <- 
  immigrant_pop_size * ISR

# # reset the starting stage distribution for simulation (all with 10 individuals)
# m <- rep(10, no_stages) 
# 
# # number of male adults at time t
# M2 <- stage[4, i] * A["M_Adult", "M_Adult"] + I[4, 1]
# 
# # number of female adults at time t
# F2 <- stage[2, i] * A["F_Adult", "F_Adult"] + I[2, 1]

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
juv_pop_storage[1, "nF"] <- N_F_eggs * egg_survival
juv_pop_storage[1, "nM"] <- N_M_eggs * egg_survival

for(i in 2:(max(average_rates$age) + 1)){
  juv_pop_storage[i, "nF"] <- 
    juv_pop_storage[i - 1, "nF"] * average_rates[which(average_rates$sex == "Female"), ][["med_surv"]][i - 1]
  juv_pop_storage[i, "nM"] <- 
    juv_pop_storage[i - 1, "nM"] * average_rates[which(average_rates$sex == "Male"), ][["med_surv"]][i - 1]
}

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

BC_JSR_run <- 
  JSR_bootstrap(niter = 1000, 
                F_rates = BC_F_rates,
                M_rates = BC_M_rates, 
                immigrant_pop_size = 100, 
                ISR = BC_ISR, h = BC_h, k = BC_k, species = "BC",
                HSR = BC_HSR, egg_survival = BC_egg_survival)

WBC_JSR_run <- 
  JSR_bootstrap(niter = 1000, 
                F_rates = WBC_F_rates,
                M_rates = WBC_M_rates, 
                immigrant_pop_size = 100, 
                ISR = WBC_ISR, h = WBC_h, k = WBC_k, species = "WBC",
                HSR = WBC_HSR, egg_survival = WBC_egg_survival)

juv_pop_storage_clean <- 
  juv_pop_storage %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "age") %>% 
  pivot_longer(!age, names_to = "sex", values_to = "n")

BC_juv_pop_storage_clean <- 
  data.frame(BC_JSR_run$Females) %>% 
  mutate(age = c(0:69)) %>% 
  # as.data.frame() %>% 
  # rownames_to_column(var = "age") %>% 
  pivot_longer(!age, names_to = "iteration", values_to = "n") %>% 
  mutate(sex = "F") %>% 
  arrange(iteration)

M_juv_pop_storage_clean <- 
  data.frame(M_juv_pop_storage_boot) %>% 
  mutate(age = c(0:69)) %>% 
  # as.data.frame() %>% 
  # rownames_to_column(var = "age") %>% 
  pivot_longer(!age, names_to = "iteration", values_to = "n") %>% 
  mutate(sex = "M") %>% 
  arrange(iteration)

ggplot(data = bind_rows(F_juv_pop_storage_clean,
                        M_juv_pop_storage_clean)) +
  geom_line(aes(y = n, x = as.numeric(age), group = iteration)) +
  facet_grid(sex ~ .)

JSR <- 
  juv_pop_storage %>% 
  as.data.frame() %>% 
  mutate(JSR = nM / (nM + nF)) %>% 
  rownames_to_column(var = "age")

JSR_boot <- 
  left_join(F_juv_pop_storage_clean, M_juv_pop_storage_clean, 
            by = c("age", "iteration"))  %>% 
  mutate(JSR = n.y / (n.y + n.x))

filter(JSR_boot, age == 69) %>% 
  summary(mean(JSR))

ggplot(data = bind_rows(BC_JSR_run, WBC_JSR_run)) +
  geom_line(aes(y = JSR, x = as.numeric(age), group = iteration),
            alpha = 0.05) +
  geom_hline(yintercept = 0.5, color = "white") +
  facet_grid(species ~ .)

ggplot(data = out) +
  geom_line(aes(y = JSR, x = as.numeric(age), group = iteration),
            alpha = 0.05) +
  geom_hline(yintercept = 0.5, color = "white")
