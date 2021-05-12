### --- two-sex toy matrix model ----

k = 4
HSR = 0.4955
h = 1/2.9
Egg_survival = 0.32
ISR = 0.738
immigrant_pop_size = 100
F_Nestling_survival = 0.741
M_Nestling_survival = 0.738
F_Groundling_survival = 0.610
M_Groundling_survival = 0.809
F_Fledgling_survival = 0.990
M_Fledgling_survival = 0.972
F_Adult_survival = 0.055
M_Adult_survival = 0.232
iterations = 10
stages <- c("F_1st_year",  "F_Adult",  "M_1st_year",  "M_Adult")

M <- 
  matrix(c(
    
    # top row of matrix
    0, NA, 0, NA, 
    
    # second row of matrix
    (Egg_survival * 
       F_Nestling_survival * 
       F_Groundling_survival *
       F_Fledgling_survival),
    F_Adult_survival,
    0, 0,
    
    # third row of matrix
    0, NA, 0, NA, 
    
    # fourth row of matrix
    0, 0, 
    (Egg_survival * 
       M_Nestling_survival * 
       M_Groundling_survival *
       M_Fledgling_survival),
    M_Adult_survival),
    nrow = length(stages), byrow = TRUE,
    dimnames = list(stages, stages))

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
# for (i in 1:t) {
  i = 2
  # stage distribution at time t
  stage[,i] <- n
  stage[, c(1:10)]
  # population size at time t
  pop[i] <- sum(n)
  pop[c(1:10)]
  
  # number of male adults at time t = (number alive at t-1) * 
  # (survival rate) + (number of immigrants entering at t)
  M2 <- (stage[4, i] * M["M_Adult", "M_Adult"]) + (I[4, 1] * M_Fledgling_survival)
  
  # number of female adults at time t
  F2 <- stage[2, i] * M["F_Adult", "F_Adult"] + I[2, 1]
  
  # Female freq-dep fecundity of Female chicks
  M["F_1st_year", "F_Adult"] <- ((k * M2) / (M2 + (F2 / h)) * (1 - HSR) )
  
  # Female freq-dep fecundity of Male chicks
  M["M_1st_year", "F_Adult"] <- ((k * M2) / (M2 + (F2 / h)) * HSR)
  
  # Male freq-dep fecundity of Female chicks
  M["F_1st_year", "M_Adult"] <- ((k * F2) / (M2 + (F2 / h)) * (1 - HSR))
  
  # Male freq-dep fecundity of Male chicks
  M["M_1st_year", "M_Adult"] <- ((k * F2) / (M2 + (F2 / h)) * HSR)
  
  # define the new n (i.e., new stage distribution at time t)
  n <- M %*% (n + I)
  
  # # number of female 1st years produced
  # N_female_1st_years <-
  #   (n[2] * M["F_1st_year", "F_Adult"] * M["F_Adult", "F_1st_year"]) +
  #   (n[4] * M["F_1st_year", "M_Adult"] * M["F_Adult", "F_1st_year"])
  # 
  # 
  # # number of male 1st years produced
  # N_male_1st_years <-
  #   (n[2] * M["M_1st_year", "F_Adult"] * M["M_Adult", "M_1st_year"]) +
  #   (n[4] * M["M_1st_year", "M_Adult"] * M["M_Adult", "M_1st_year"])
  # # 
  # sex_ratio_1st_years <- N_male_1st_years / (N_male_1st_years + N_female_1st_years)
  
  pop_dist[, i] <- n
  
# }
# define rownames of stage matrix
rownames(stage) <- rownames(M)

# define colnames of stage matrix
colnames(stage) <- 0:(t - 1)

# calculate the proportional stable stage distribution
stage <- apply(stage, 2, function(x) x/sum(x))

# define stable stage as the last stage
stable.stage <- stage[, t]
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
                 SSD_F2 = stable.stage[2]
                 # iter = num_boot + ((iter_add - 1) * niter), 
                 # species = species
                 )

# plot distrubution to assure that it is not chaotic

palette_custom <- c("#FC8D62", "#FC8D62", "#66C2A5", "#66C2A5")
linetype_custom <- c(2, 1, 2, 1)

SSD_prop_plot <-
  pop.proj$stage.vectors %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(iteration = as.numeric(rownames(.))) %>% 
  pivot_longer(-iteration, names_to = "stage", values_to = "proportion") %>% 
  mutate(stage = as.factor(as.character(factor(stage)))) %>% 
  ggplot(data = ., aes(x = iteration, y = proportion, group = stage, color = stage, linetype = stage)) +
  geom_line() +
  scale_color_manual(values = palette_custom) +
  scale_linetype_manual(values = linetype_custom) +
  xlab("Time-step") +
  ylab("Proporational stage distribution") +
  annotate(geom = "text", y = 0.2, x = 0.8*ncol(pop.proj$stage_size),
           label = paste("First-year sex ratio = ", round(pop.proj$ASR, digits = 3)),
           color = "black", size = 3) +
  annotate(geom = "text", y = 0.1, x = 0.8*ncol(pop.proj$stage_size),
           label = paste("Hatching sex ratio = ", round(pop.proj$SSR, digits = 3)),
           color = "black", size = 3)

SSD_pop_plot <-
  pop.proj$stage_size %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(iteration = as.numeric(rownames(.))) %>% 
  pivot_longer(-iteration, names_to = "stage", values_to = "population") %>% 
  mutate(stage = as.factor(as.character(factor(stage)))) %>% 
  ggplot(data = ., aes(x = iteration, y = population, group = stage, color = stage, linetype = stage)) +
  geom_line() +
  scale_color_manual(values = palette_custom) +
  scale_linetype_manual(values = linetype_custom) +
  xlab("Time-step") +
  ylab("Absolute stage distribution") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top") +
  annotate(geom = "text", y = 30, x = 0.8*ncol(pop.proj$stage_size),
           label = paste("Lamda = ", round(pop.proj$lambda, digits = 3)),
           color = "black", size = 3)

plot_result <- 
  ggarrange(SSD_pop_plot, 
            SSD_prop_plot,
            nrow = 2, align = "v",
            heights = c(0.5, 0.5))
plot_result

# print the list as output to the function
pop.proj
# 
# adult_F_immigrants <- immigrant_pop_size * (1 - ISR)
# adult_M_immigrants <- immigrant_pop_size * ISR
# 
# coucal_adult_immigration <- 
#   data.frame(sex = c("Female", "Male"),
#              value = c(adult_F_immigrants, adult_M_immigrants),
#              stage = c("adult"),
#              rate = c("immigration"))
# 
# coucal_egg_survival <- 
#   data.frame(sex = NA,
#              value = egg_survival,
#              stage = c("egg"),
#              rate = c("survival"))
# 
# # Bind the juvenile and adult dataframe with the nestlings
# coucal_vital_rates <- 
#   bind_rows(coucal_egg_survival,
#             coucal_nestling_survival,
#             coucal_groundling_survival,
#             coucal_fledgling_survival,
#             coucal_adult_survival,
#             coucal_adult_immigration) %>% 
#   as.data.frame() %>% 
#   mutate(iter = num_boot + ((iter_add - 1) * niter),
#          species = species)
# 
# # Create a list of demographic rates from the survival analyses above
# demographic_rates <- list(Egg_survival = coucal_vital_rates[1, 2],
#                           F_Nestling_survival = coucal_vital_rates[2, 2],
#                           F_Groundling_survival = coucal_vital_rates[4, 2],
#                           F_Fledgling_survival = coucal_vital_rates[6, 2],
#                           F_Adult_survival = coucal_vital_rates[8, 2],
#                           F_Adult_immigration = coucal_vital_rates[10, 2],
#                           M_Nestling_survival = coucal_vital_rates[3, 2],
#                           M_Groundling_survival = coucal_vital_rates[5, 2],
#                           M_Fledgling_survival = coucal_vital_rates[7, 2],
#                           M_Adult_survival = coucal_vital_rates[9, 2],
#                           M_Adult_immigration = coucal_vital_rates[11, 2],
#                           
#                           # Define hatching sex ratio
#                           HSR = HSR,
#                           
#                           # Define the mating system (h), and clutch size (k)
#                           h = h,
#                           k = k)
# 
# # Build matrix based on rates specified in the list above
# demographic_matrix <- coucal_matrix(demographic_rates)
# 
# # populate sex-specific adult immigration matrix
# immigration_vector <- matrix(data = c(0, demographic_rates$F_Adult_immigration, 
#                                       0, demographic_rates$M_Adult_immigration), 
#                              nrow = 4, ncol = 1)
# 
# # Determine the ASR at the stable stage distribution
# ASR_SSD <- matrix_ASR(M = demographic_matrix,
#                       h = demographic_rates$h,
#                       HSR = demographic_rates$HSR, iterations = 1000,
#                       num_boot = num_boot,
#                       species = species,
#                       immigrant_pop_size = immigrant_pop_size,
#                       ISR = ISR, iter_add = 1)


