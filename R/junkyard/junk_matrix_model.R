coucal_vital_rates <- 
  BC_hazard_rate_boot_tidy$vital_rate_ests_boot %>% 
  dplyr::group_by(sex, stage, rate) %>% 
  dplyr::summarise(value = mean(value, na.rm = TRUE)) %>% 
  data.frame()

# Create a list of demographic rates from the survival analyses above
demographic_rates <- list(Egg_survival = coucal_vital_rates[11, 4],
                          F_Nestling_survival = coucal_vital_rates[5, 4],
                          F_Groundling_survival = coucal_vital_rates[4, 4],
                          F_Fledgling_survival = coucal_vital_rates[3, 4],
                          F_Adult_survival = coucal_vital_rates[2, 4],
                          F_Adult_immigration = coucal_vital_rates[1, 4],
                          M_Nestling_survival = coucal_vital_rates[10, 4],
                          M_Groundling_survival = coucal_vital_rates[9, 4],
                          M_Fledgling_survival = coucal_vital_rates[8, 4],
                          M_Adult_survival = coucal_vital_rates[7, 4],
                          M_Adult_immigration = coucal_vital_rates[6, 4],
                          
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
M = demographic_matrix
h = demographic_rates$h
HSR = demographic_rates$HSR
k = demographic_rates$k
iterations = 100
num_boot = 1
species = "BC"
immigrant_pop_size = 100
ISR = BC_ISR
iter_add = 1
n = rep(10, nrow(M))

# Number of stages in matrix
x <- length(n)

# Number of time steps to simulate
t <- iterations

# an empty t by x matrix to store the stage distributions
stage <- matrix(data = numeric(x * t), nrow = x, ncol = t)

# an empty t by x matrix to store the stage-specific frequencies
pop_dist <- matrix(data = numeric(x * t), nrow = x, ncol = t)
rownames(pop_dist) <- rownames(M)
colnames(pop_dist) <- 0:(t - 1)

# an empty t vector to store the population sizes
pop <- numeric(t)

I <- matrix(data = c(0, immigrant_pop_size * (1 - ISR), 
                     0, immigrant_pop_size * ISR), 
            nrow = 4, ncol = 1)

# for loop that goes through each of t time steps
for (i in 1:t) {
  
  # stage distribution at time t
  stage[,i] <- n
  
  # population size at time t
  pop[i] <- sum(n)
  
  # number of male adults at time t = (number alive at t-1) * 
  # (survival rate) + (number of immigrants entering at t)
  M2 <- stage[4, i] * M["M_Adult", "M_Adult"] + I[4, 1]
  
  # number of female adults at time t
  F2 <- stage[2, i] * M["F_Adult", "F_Adult"] + I[2, 1]
  
  # Female freq-dep fecundity of Female chicks
  M["F_1st_year", "F_Adult"] <- ((k * M2) / (M2 + (F2 / h)) * (1 - HSR))
  # Female freq-dep fecundity of Female chicks
  M["F_1st_year", "F_1st_year"] <- ((k * M2) / (M2 + (F2 / h)) * (1 - HSR))
  
  # Female freq-dep fecundity of Male chicks
  M["M_1st_year", "F_Adult"] <- ((k * M2) / (M2 + (F2 / h)) * HSR)
  # Female freq-dep fecundity of Male chicks
  M["M_1st_year", "M_1st_year"] <- ((k * M2) / (M2 + (F2 / h)) * HSR)
  
  # Male freq-dep fecundity of Female chicks
  M["F_1st_year", "M_Adult"] <- ((k * F2) / (M2 + (F2 / h)) * (1 - HSR))
  # Male freq-dep fecundity of Female chicks
  M["F_1st_year", "F_1st_year"] <- ((k * F2) / (M2 + (F2 / h)) * (1 - HSR))
  
  # Male freq-dep fecundity of Male chicks
  M["M_1st_year", "M_Adult"] <- ((k * F2) / (M2 + (F2 / h)) * HSR)
  # Male freq-dep fecundity of Male chicks
  M["M_1st_year", "M_1st_year"] <- ((k * F2) / (M2 + (F2 / h)) * HSR)
  
  # define the new n (i.e., new stage distribution at time t)
  n <- M %*% (n + I)
  
  # number of female 1st years produced
  N_female_1st_years <-
    (n[2] * M["F_1st_year", "F_Adult"] * M["F_Adult", "F_1st_year"]) +
    (n[4] * M["F_1st_year", "M_Adult"] * M["F_Adult", "F_1st_year"])

  # number of male 1st years produced
  N_male_1st_years <-
    (n[2] * M["M_1st_year", "F_Adult"] * M["M_Adult", "M_1st_year"]) +
    (n[4] * M["M_1st_year", "M_Adult"] * M["M_Adult", "M_1st_year"])
  
  sex_ratio_1st_years <- N_male_1st_years / (N_male_1st_years + N_female_1st_years)
  
  n[1]/(n[1]+n[3])
  
  # define rownames of stage matrix
  rownames(stage) <- rownames(M)
  
  # define colnames of stage matrix
  colnames(stage) <- 0:(t - 1)
  
  # calculate the proportional stable stage distribution
  stage <- apply(stage, 2, function(x) x/sum(x))
  
  # define stable stage as the last stage
  stable.stage <- stage[, t]
  
  pop_dist[, i] <- n
  
}

# calc ASR as the proportion of the adult stable stage class that is male
ASR <- stable.stage[4]/(stable.stage[2] + stable.stage[4])
SSR <- stable.stage[3]/(stable.stage[1] + stable.stage[3])

# make a list of results
pop.proj <- list(ASR = ASR,
                 SSR = SSR,
                 lambda = pop[t]/pop[t - 1],
                 pop_size = pop,
                 stage_size = pop_dist,
                 stable.stage = stable.stage,
                 stage.vectors = stage,
                 SSD_M2 = stable.stage[4],
                 SSD_F2 = stable.stage[2],
                 iter = num_boot + ((iter_add - 1) * niter), 
                 species = species)

# print the list as output to the function
pop.proj