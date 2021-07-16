
source("R/project/project_libraries.R")

# coucal_matrix() builds the two-sex Lefkovitch matrix using the vital rates 
# specified in the *demographic_rates* object.

popA_F_disp_rate = 0.9
popA_M_disp_rate = 0.9
popB_F_disp_rate = 0.9
popB_M_disp_rate = 0.9
popA_F_ad_surv = 0.65
popA_M_ad_surv = 0.65
popB_F_ad_surv = 0.65
popB_M_ad_surv = 0.65
hazard_rate_boot_tidy = BC_hazard_rate_boot_tidy
k = 4
HSR = 0.4955
h = 1/2.9
species = "BC"
ISR = 0.738
iterations = 10
n = rep(10, 4)
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
    
    stages <- c("F_PopA",  "M_PopA",  "F_PopB",  "M_PopB")
    
    n <- n * c((1 - ISR), ISR, (1 - ISR), ISR)

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
    popA <- numeric(t)
    popB <- numeric(t)
    
    rand_iter <- sample(1:1000, 2, replace = TRUE)
    
    popA_rates <- 
      hazard_rate_boot_tidy$vital_rate_ests_boot %>% 
      dplyr::filter(iter == rand_iter[1])
    
    popB_rates <- 
      hazard_rate_boot_tidy$vital_rate_ests_boot %>% 
      dplyr::filter(iter == rand_iter[2])
    
    # # for loop that goes through each of t time steps
    for (i in 1:t) {

      # Build the 4x4 matrix
      M <- 
        matrix(c(
          
          #### Row 1 ----
          # cell [1, 1] (Female production of females that stay locally in PopA)
          popA_rates[which(popA_rates$stage == "egg"), "value"] *
            popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "nestling"), "value"] *
            popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "groundling"), "value"] *
            popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "fledgling"), "value"] *
            (1 - HSR) * (1 - popA_F_disp_rate),
          
          # cell [1, 2] (Male production of females that stay locally in PopA)
          popA_rates[which(popA_rates$stage == "egg"), "value"] *
            popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "nestling"), "value"] *
            popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "groundling"), "value"] *
            popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "fledgling"), "value"] *
            (1 - HSR) * (1 - popA_F_disp_rate),
          
          # cell [1, 3] (Female production of females that disperse from PopB)
          popB_rates[which(popB_rates$stage == "egg"), "value"] *
            popB_rates[which(popB_rates$sex == "Female" & popB_rates$stage == "nestling"), "value"] *
            popB_rates[which(popB_rates$sex == "Female" & popB_rates$stage == "groundling"), "value"] *
            popB_rates[which(popB_rates$sex == "Female" & popB_rates$stage == "fledgling"), "value"] *
            (1 - HSR) * popB_F_disp_rate,
          
          # cell [1, 4] (Male production of females that disperse from PopB)
          popB_rates[which(popB_rates$stage == "egg"), "value"] *
            popB_rates[which(popB_rates$sex == "Female" & popB_rates$stage == "nestling"), "value"] *
            popB_rates[which(popB_rates$sex == "Female" & popB_rates$stage == "groundling"), "value"] *
            popB_rates[which(popB_rates$sex == "Female" & popB_rates$stage == "fledgling"), "value"] *
            (1 - HSR) * popB_F_disp_rate,
          
          #### Row 2 ----
          # cell [2, 1] (Female production of males that stay locally in PopA)
          popA_rates[which(popA_rates$stage == "egg"), "value"] *
            popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "nestling"), "value"] *
            popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "groundling"), "value"] *
            popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "fledgling"), "value"] *
            HSR * (1 - popA_M_disp_rate),
          
          # cell [2, 2] (Male production of males that stay locally in PopA)
          popA_rates[which(popA_rates$stage == "egg"), "value"] *
            popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "nestling"), "value"] *
            popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "groundling"), "value"] *
            popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "fledgling"), "value"] *
            HSR * (1 - popA_M_disp_rate),
          
          # cell [2, 3] (Female production of males that disperse from PopB)
          popB_rates[which(popB_rates$stage == "egg"), "value"] *
            popB_rates[which(popB_rates$sex == "Male" & popB_rates$stage == "nestling"), "value"] *
            popB_rates[which(popB_rates$sex == "Male" & popB_rates$stage == "groundling"), "value"] *
            popB_rates[which(popB_rates$sex == "Male" & popB_rates$stage == "fledgling"), "value"] *
            HSR * popB_M_disp_rate,
          
          # cell [2, 4] (Male production of males that disperse from PopB)
          popB_rates[which(popB_rates$stage == "egg"), "value"] *
            popB_rates[which(popB_rates$sex == "Male" & popB_rates$stage == "nestling"), "value"] *
            popB_rates[which(popB_rates$sex == "Male" & popB_rates$stage == "groundling"), "value"] *
            popB_rates[which(popB_rates$sex == "Male" & popB_rates$stage == "fledgling"), "value"] *
            HSR * popB_M_disp_rate,
          
          #### Row 3 ----
          # cell [3, 1] (Female production of females that disperse from PopA)
          popA_rates[which(popA_rates$stage == "egg"), "value"] *
            popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "nestling"), "value"] *
            popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "groundling"), "value"] *
            popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "fledgling"), "value"] *
            (1 - HSR) * popA_F_disp_rate,
          
          # cell [3, 2] (Male production of females that disperse from PopA)
          popA_rates[which(popA_rates$stage == "egg"), "value"] *
            popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "nestling"), "value"] *
            popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "groundling"), "value"] *
            popA_rates[which(popA_rates$sex == "Female" & popA_rates$stage == "fledgling"), "value"] *
            (1 - HSR) * popA_F_disp_rate,
          
          # cell [3, 3] (Female production of females that stay locally in PopB)
          popB_rates[which(popB_rates$stage == "egg"), "value"] *
            popB_rates[which(popB_rates$sex == "Female" & popB_rates$stage == "nestling"), "value"] *
            popB_rates[which(popB_rates$sex == "Female" & popB_rates$stage == "groundling"), "value"] *
            popB_rates[which(popB_rates$sex == "Female" & popB_rates$stage == "fledgling"), "value"] *
            (1 - HSR) * (1 - popB_F_disp_rate),
          
          # cell [3, 4] (Male production of females that stay locally in PopB)
          popB_rates[which(popB_rates$stage == "egg"), "value"] *
            popB_rates[which(popB_rates$sex == "Female" & popB_rates$stage == "nestling"), "value"] *
            popB_rates[which(popB_rates$sex == "Female" & popB_rates$stage == "groundling"), "value"] *
            popB_rates[which(popB_rates$sex == "Female" & popB_rates$stage == "fledgling"), "value"] *
            (1 - HSR) * (1 - popB_F_disp_rate),
          
          #### Row 4 ----
          
          # cell [4, 1] (Female production of males that disperse from PopA)
          popA_rates[which(popA_rates$stage == "egg"), "value"] *
            popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "nestling"), "value"] *
            popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "groundling"), "value"] *
            popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "fledgling"), "value"] *
            HSR * popA_M_disp_rate,
          
          # cell [4, 2] (Male production of males that disperse from PopA)
          popA_rates[which(popA_rates$stage == "egg"), "value"] *
            popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "nestling"), "value"] *
            popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "groundling"), "value"] *
            popA_rates[which(popA_rates$sex == "Male" & popA_rates$stage == "fledgling"), "value"] *
            HSR * popA_M_disp_rate,
          
          # cell [4, 3] (Female production of males that stay locally in PopB)
          popB_rates[which(popB_rates$stage == "egg"), "value"] *
            popB_rates[which(popB_rates$sex == "Male" & popB_rates$stage == "nestling"), "value"] *
            popB_rates[which(popB_rates$sex == "Male" & popB_rates$stage == "groundling"), "value"] *
            popB_rates[which(popB_rates$sex == "Male" & popB_rates$stage == "fledgling"), "value"] *
            HSR * (1 - popB_M_disp_rate),
          
          # cell [4, 4] (Male production of males that stay locally in PopB)
          popB_rates[which(popB_rates$stage == "egg"), "value"] *
            popB_rates[which(popB_rates$sex == "Male" & popB_rates$stage == "nestling"), "value"] *
            popB_rates[which(popB_rates$sex == "Male" & popB_rates$stage == "groundling"), "value"] *
            popB_rates[which(popB_rates$sex == "Male" & popB_rates$stage == "fledgling"), "value"] *
            HSR * (1 - popB_M_disp_rate)
          ),
          
          #### tidy matrix ----
          nrow = length(stages), byrow = TRUE,
          dimnames = list(stages, stages))
      
      # stage distribution at time t
      stage[,i] <- n
      
      # metapopulation size at time t
      pop[i] <- sum(n)
      
      # popA size at time t
      popA[i] <- stage[1, i] + stage[2, i]
      
      # popB size at time t
      popB[i] <- stage[3, i] + stage[4, i]
      
      # number of female at time t
      popA_nF <- stage[1, i]
      
      # number of males at time t
      popA_nM <- stage[2, i]
      
      # number of female at time t
      popB_nF <- stage[3, i]
      
      # number of males at time t
      popB_nM <- stage[4, i]
      
      # PopA Female freq-dep fecundity
      popA_F_fert <- (k * popA_nM) / (popA_nM + (popA_nF / h))
      popA_bfun <- (2 * k * popA_nM * popA_nF) / (popA_nM + (popA_nF / h))
      
      M["F_PopA", "F_PopA"] <- M["F_PopA", "F_PopA"] * popA_bfun * popA_F_ad_surv
      M["M_PopA", "F_PopA"] <- M["M_PopA", "F_PopA"] * popA_bfun * popA_M_ad_surv
      M["F_PopB", "F_PopA"] <- M["F_PopB", "F_PopA"] * popA_bfun * popA_F_ad_surv
      M["M_PopB", "F_PopA"] <- M["M_PopB", "F_PopA"] * popA_bfun * popA_M_ad_surv
      
      # PopA Male freq-dep fecundity
      popA_M_fert <- (k * popA_nF) / (popA_nM + (popA_nF / h))
      M["M_PopA", "M_PopA"] <- M["M_PopA", "M_PopA"] * popA_bfun * popA_M_ad_surv
      M["F_PopA", "M_PopA"] <- M["F_PopA", "M_PopA"] * popA_bfun * popA_F_ad_surv
      M["M_PopB", "M_PopA"] <- M["M_PopB", "M_PopA"] * popA_bfun * popA_M_ad_surv
      M["F_PopB", "M_PopA"] <- M["F_PopB", "M_PopA"] * popA_bfun * popA_F_ad_surv
      
      # popB Female freq-dep fecundity
      popB_F_fert <- (k * popB_nM) / (popB_nM + (popB_nF / h))
      popB_bfun <- (2 * k * popB_nM * popB_nF) / (popB_nM + (popB_nF / h))
      
      M["F_PopB", "F_PopB"] <- M["F_PopB", "F_PopB"] * popB_bfun * popB_F_ad_surv
      M["M_PopB", "F_PopB"] <- M["M_PopB", "F_PopB"] * popB_bfun * popB_M_ad_surv
      M["F_PopA", "F_PopB"] <- M["F_PopA", "F_PopB"] * popB_bfun * popB_F_ad_surv
      M["M_PopA", "F_PopB"] <- M["M_PopA", "F_PopB"] * popB_bfun * popB_M_ad_surv
      
      # popB Male freq-dep fecundity
      popB_M_fert <- (k * popB_nF) / (popB_nM + (popB_nF / h))
      M["M_PopB", "M_PopB"] <- M["M_PopB", "M_PopB"] * popB_bfun * popB_M_ad_surv
      M["F_PopB", "M_PopB"] <- M["F_PopB", "M_PopB"] * popB_bfun * popB_F_ad_surv
      M["M_PopA", "M_PopB"] <- M["M_PopA", "M_PopB"] * popB_bfun * popB_M_ad_surv
      M["F_PopA", "M_PopB"] <- M["F_PopA", "M_PopB"] * popB_bfun * popB_F_ad_surv
      
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
    
    mASR <- (stage[2, t] + stage[4, t]) / (stage[1, t] + stage[2, t] + stage[3, t] + stage[4, t])
    popA_ASR <- (stage[2, t]) / (stage[1, t] + stage[2, t])
    popB_ASR <- (stage[4, t]) / (stage[3, t] + stage[4, t])
    
    # make a list of results
    pop.proj <- list(mASR = mASR,
                     popA_ASR = popA_ASR,
                     popB_ASR = popB_ASR,
                     lambda = pop[t]/pop[t - 1],
                     popA_lambda =  popA[t]/popA[t - 1],
                     popB_lambda =  popB[t]/popB[t - 1],
                     metapop_size = pop,
                     popA_size = popA,
                     popB_size = popB,
                     stage_size = pop_dist,
                     stable.stage = stable.stage,
                     stage.vectors = stage,
                     # iter = num_boot + ((iter_add - 1) * niter), 
                     species = species)
    
    # print the list as output to the function
    pop.proj
  # }
    
    mASR
    popA_ASR
    popB_ASR
    
    
    #### JUNK STUFF? ----
    # # two population matrix model
    # 
    # source("R/project/project_libraries.R")
    # 
    # # coucal_matrix() builds the two-sex Lefkovitch matrix using the vital rates 
    # # specified in the *demographic_rates* object.
    # 
    # meta_coucal_matrix <- 
    #   function(demographic_rates, two_sex = TRUE, two_pop = FALSE){
    #     if(two_sex){
    #       
    #       # Define coucal life-stages of the coucal matrix model
    #       stages <- c("F_1st_year",  "F_Adult",  "M_1st_year",  "M_Adult")
    #       
    #       # Build the 4x4 matrix
    #       result <- 
    #         matrix(c(
    #           
    #           # top row of matrix
    #           0, NA, 0, NA, 
    #           
    #           # second row of matrix
    #           (demographic_rates$Egg_survival * 
    #              demographic_rates$F_Nestling_survival * 
    #              demographic_rates$F_Groundling_survival *
    #              demographic_rates$F_Fledgling_survival),
    #           demographic_rates$F_Adult_survival,
    #           0, 0,
    #           
    #           # third row of matrix
    #           0, NA, 0, NA, 
    #           
    #           # fourth row of matrix
    #           0, 0, 
    #           (demographic_rates$Egg_survival * 
    #              demographic_rates$M_Nestling_survival * 
    #              demographic_rates$M_Groundling_survival *
    #              demographic_rates$M_Fledgling_survival),
    #           demographic_rates$M_Adult_survival),
    #           nrow = length(stages), byrow = TRUE,
    #           dimnames = list(stages, stages))
    #     }
    #     if(two_pop){
    #       # Define coucal life-stages of the coucal matrix model
    #       stages <- c("F_1st_year_A",  "F_Adult_A",  "M_1st_year_A",  "M_Adult_A",
    #                   "F_1st_year_B",  "F_Adult_B",  "M_1st_year_B",  "M_Adult_B")
    #       
    #       # Build the 4x4 matrix
    #       result <- 
    #         matrix(c(
    #           
    #           # top row of matrix
    #           0, NA, 0, NA, 
    #           
    #           # second row of matrix
    #           (demographic_rates$Egg_survival * 
    #              demographic_rates$F_Nestling_survival * 
    #              demographic_rates$F_Groundling_survival *
    #              demographic_rates$F_Fledgling_survival),
    #           demographic_rates$F_Adult_survival,
    #           0, 0,
    #           
    #           # third row of matrix
    #           0, NA, 0, NA, 
    #           
    #           # fourth row of matrix
    #           0, 0, 
    #           (demographic_rates$Egg_survival * 
    #              demographic_rates$M_Nestling_survival * 
    #              demographic_rates$M_Groundling_survival *
    #              demographic_rates$M_Fledgling_survival),
    #           demographic_rates$M_Adult_survival),
    #           nrow = length(stages), byrow = TRUE,
    #           dimnames = list(stages, stages))
    #     }
    #     else{
    #       
    #       # Define coucal life-stages of the Ceuta snowy coucal matrix model
    #       stages <- c("1st_year",  "Adult")
    #       
    #       # Build the 4x4 matrix
    #       result <- 
    #         matrix(c(
    #           
    #           # top row of matrix
    #           0, NA,
    #           
    #           # second row of matrix
    #           (demographic_rates$Egg_survival * 
    #              demographic_rates$Nestling_survival * 
    #              demographic_rates$Groundling_survival *
    #              demographic_rates$Fledgling_survival),
    #           demographic_rates$Adult_survival),
    #           nrow = length(stages), byrow = TRUE,
    #           dimnames = list(stages, stages))
    #     }
    #     result
    #   }