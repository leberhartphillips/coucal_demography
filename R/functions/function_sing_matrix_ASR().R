
source("R/project/project_libraries.R")

# coucal_matrix() builds the two-sex Lefkovitch matrix using the vital rates 
# specified in the *demographic_rates* object.

F_ad_surv = 0.799025
M_ad_surv = 0.799025
hazard_rate_boot_tidy = BC_hazard_rate_boot_tidy
k = 4
HSR = 0.4955
h = 1/2.9
species = "BC"
ISR = 0.738
iterations = 10
n = rep(10, 2)
i = 1
# meta_matrix_ASR <-
#   function(hazard_rate_boot_tidy,
#            n = rep(100, 4), 
#            h, k, HSR,
#            iterations = 10,
#            species,
#            popA_F_disp_rate = 1,
#            popA_M_disp_rate = 1,
#            popB_F_disp_rate = 1,
#            popB_M_disp_rate = 1,
#            ISR){

stages <- c("F",  "M")

n <- n * c((1 - ISR), ISR)

# Number of stages in matrix
x <- length(n)

# Number of time steps to simulate
t <- iterations

# an empty t by x matrix to store the stage distributions
stage <- matrix(data = numeric(x * t), nrow = x, ncol = t)

# an empty t by x matrix to store the stage-specific frequencies
pop_dist <- matrix(data = numeric(x * t), nrow = x, ncol = t)
rownames(pop_dist) <- stages
colnames(pop_dist) <- 0:(t - 1)

# an empty t vector to store the population sizes
pop <- numeric(t)

rand_iter <- sample(1:1000, 1)

pop_rates <- 
  hazard_rate_boot_tidy$vital_rate_ests_boot %>% 
  dplyr::filter(iter == rand_iter[1])

# for loop that goes through each of t time steps
for (i in 1:t) {
  
  # Build the 2x2 matrix
  M <- 
    matrix(c(
      
      #### Row 1 ----
      # cell [1, 1] (Female production of females)
      popA_rates[which(popA_rates$stage == "egg"), "value"] *
        popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "nestling"), "value"] *
        popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "groundling"), "value"] *
        popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "fledgling"), "value"] *
        (1 - HSR),
      
      # cell [1, 2] (Male production of females)
      popA_rates[which(popA_rates$stage == "egg"), "value"] *
        popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "nestling"), "value"] *
        popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "groundling"), "value"] *
        popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "fledgling"), "value"] *
        (1 - HSR),
      
      #### Row 2 ----
      # cell [2, 1] (Female production of males that stay locally in PopA)
      popA_rates[which(popA_rates$stage == "egg"), "value"] *
        popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "nestling"), "value"] *
        popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "groundling"), "value"] *
        popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "fledgling"), "value"] *
        HSR,
      
      # cell [2, 2] (Male production of males that stay locally in PopA)
      popA_rates[which(popA_rates$stage == "egg"), "value"] *
        popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "nestling"), "value"] *
        popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "groundling"), "value"] *
        popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "fledgling"), "value"] *
        HSR
    ),
    
    #### tidy matrix ----
    nrow = length(stages), byrow = TRUE,
    dimnames = list(stages, stages))
  
  # stage distribution at time t
  stage[,i] <- n
  
  # population size at time t
  pop[i] <- sum(n)
  
  # number of female at time t
  nF <- stage[1, i]
  
  # number of males at time t
  nM <- stage[2, i]
  
  # PopA Female freq-dep fecundity
  F_fert <- (k * nM) / (nM + (nF / h))
  M_fert <- (k * nF) / (nM + (nF / h))
  bfun <- (2 * k * nM * nF) / (nM + (nF / h))
  
  M["F", "F"] <- M["F", "F"] * bfun * F_ad_surv
  M["M", "F"] <- M["M", "F"] * bfun * M_ad_surv
  
  # PopA Male freq-dep fecundity
  M["M", "M"] <- M["M", "M"] * bfun * M_ad_surv
  M["F", "M"] <- M["F", "M"] * bfun * F_ad_surv

  
  # M["F_1st_year", "F_Adult"] <- ((k * M2) / (M2 + (F2 / h)) * (1 - HSR))
  # 
  # # Female freq-dep fecundity of Male chicks
  # M["M_1st_year", "F_Adult"] <- ((k * M2) / (M2 + (F2 / h)) * HSR)
  # 
  # # Male freq-dep fecundity of Female chicks
  # M["F_1st_year", "M_Adult"] <- ((k * F2) / (M2 + (F2 / h)) * (1 - HSR))
  # 
  # # Male freq-dep fecundity of Male chicks
  # M["M_1st_year", "M_Adult"] <- ((k * F2) / (M2 + (F2 / h)) * HSR)
  
  # define the new n (i.e., new stage distribution at time t)
  n <- M %*% n
  
  # # set carrying capacity to 1000
  # if(n["F_PopA", 1] > 1000){n["F_PopA", 1] <- 1000}
  # if(n["M_PopA", 1] > 1000){n["M_PopA", 1] <- 1000}
  # if(n["F_PopB", 1] > 1000){n["F_PopB", 1] <- 1000}
  # if(n["M_PopB", 1] > 1000){n["M_PopB", 1] <- 1000}
  
  pop_dist[, i] <- n
}
pop_dist
# define rownames of stage matrix
rownames(stage) <- rownames(M)

# define colnames of stage matrix
colnames(stage) <- 0:(t - 1)

# calculate the proportional stable stage distribution
stage <- apply(stage, 2, function(x) x/sum(x))

# define stable stage as the last stage
stable.stage <- stage[, t]

# calc ASR as the proportion of the adult stable stage class that is male
# mASR <- (stable.stage[2] + stable.stage[4]) / (stable.stage[1] + stable.stage[2] + stable.stage[3] + stable.stage[4])
# popA_ASR <- (stable.stage[2]) / (stable.stage[1] + stable.stage[2])
# popB_ASR <- (stable.stage[4]) / (stable.stage[3] + stable.stage[4])

ASR <- (stage[2, t]) / (stage[1, t] + stage[2, t])

# make a list of results
pop.proj <- list(ASR = ASR,
                 lambda = pop[t]/pop[t - 1],
                 pop = pop,
                 stage_size = pop_dist,
                 stable.stage = stable.stage,
                 stage.vectors = stage,
                 # iter = num_boot + ((iter_add - 1) * niter), 
                 species = species, 
                 rand_iter = rand_iter)

# print the list as output to the function
pop.proj
# }

pop.proj$ASR
pop.proj$lambda