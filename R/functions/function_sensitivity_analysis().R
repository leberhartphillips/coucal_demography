# sensitivity_analysis() takes the vital rate summary of the bootstrap procedure
# and conducts a perturbation analysis on each rate to assess how a proportional
# change in a given vital rate changes the ASR

sensitivity_analysis <-
  function(vital_rate_summary, matrix_str, h = 1, k = 4, 
           HSR, ISR, niter = 1000, ASR, lambda, immigrant_pop_size = 100){
    
    # make a list of all parameters
    vr <-
      list(F_Nestling_survival = vital_rate_summary$F_Nestling_survival,
           F_Groundling_survival = vital_rate_summary$F_Groundling_survival,
           F_Fledgling_survival = vital_rate_summary$F_Fledgling_survival,
           F_Adult_survival = vital_rate_summary$F_Adult_survival,
           M_Nestling_survival = vital_rate_summary$M_Nestling_survival,
           M_Groundling_survival = vital_rate_summary$M_Groundling_survival,
           M_Fledgling_survival = vital_rate_summary$M_Fledgling_survival,
           M_Adult_survival = vital_rate_summary$M_Adult_survival)
    
    # number of stages in the matrix
    no_stages <- sqrt(length(matrix_str))
    
    # Define plover life-stages of the Ceuta snowy plover matrix model
    stages <- c("F_1st_year",  "F_Adult",  "M_1st_year",  "M_Adult")
    
    # an empty t by x matrix
    stage <- matrix(numeric(no_stages * niter), nrow = no_stages)
    
    # an empty t vector to store the population sizes
    pop <- numeric(niter)
    
    # dataframe to store the perturbation results
    ASR_pert_results <-
      data.frame(parameter = c("F_Nestling_survival", "F_Groundling_survival", 
                               "F_Fledgling_survival", "F_Adult_survival",
                               "M_Nestling_survival", "F_Groundling_survival", 
                               "M_Fledgling_survival", "M_Adult_survival",
                               "h", "k", "HSR", "ISR"),
                 sensitivities = numeric(length(vr) + 4),
                 elasticities = numeric(length(vr) + 4))
    
    lambda_pert_results <-
      data.frame(parameter = c("F_Nestling_survival", "F_Groundling_survival", 
                               "F_Fledgling_survival", "F_Adult_survival",
                               "M_Nestling_survival", "F_Groundling_survival", 
                               "M_Fledgling_survival", "M_Adult_survival",
                               "h", "k", "HSR", "ISR"),
                 sensitivities = numeric(length(vr) + 4),
                 elasticities = numeric(length(vr) + 4))
    
    # specifiy how many survival rates there are
    n <- length(vr)
    
    # create vectors of perturbations to test on parameters of the matrix model
    vr_nums <- seq(0, 1, 0.01) # proportional changes in survival and HSR (i.e., between 0 an 1)
    h_nums <- seq(0, 2, 0.02) # proportional changes in h index (i.e., between 0 and 2)
    k_nums <- seq(3, 5, 0.02) # proportional changes in k (i.e, between 3 and 5)
    
    # create empty dataframes to store the perturbation results for ASR and lambda
    vr_pert_ASR <- matrix(numeric(n * length(vr_nums)),
                          ncol = n, dimnames = list(vr_nums, names(vr)))
    
    h_pert_ASR <- matrix(numeric(length(h_nums)),
                         ncol = 1, dimnames = list(h_nums, "h"))
    
    k_pert_ASR <- matrix(numeric(length(k_nums)),
                         ncol = 1, dimnames = list(k_nums, "k"))
    
    HSR_pert_ASR <- matrix(numeric(length(vr_nums)),
                           ncol = 1, dimnames = list(vr_nums, "HSR"))
    
    ISR_pert_ASR <- matrix(numeric(length(vr_nums)),
                           ncol = 1, dimnames = list(vr_nums, "ISR"))
    
    vr_pert_lambda <- matrix(numeric(n * length(vr_nums)),
                             ncol = n, dimnames = list(vr_nums, names(vr)))
    
    h_pert_lambda <- matrix(numeric(length(h_nums)),
                            ncol = 1, dimnames = list(h_nums, "h"))
    
    k_pert_lambda <- matrix(numeric(length(k_nums)),
                            ncol = 1, dimnames = list(k_nums, "k"))
    
    HSR_pert_lambda <- matrix(numeric(length(vr_nums)),
                              ncol = 1, dimnames = list(vr_nums, "HSR"))
    
    ISR_pert_lambda <- matrix(numeric(length(vr_nums)),
                              ncol = 1, dimnames = list(vr_nums, "ISR"))
    
    ##### perturbation of survival rates ####
    for (g in 1:n) # pick a column (i.e., a variable)
    {
      vr2 <- vr # reset the vital rates to the original
      
      for (i in 1:length(vr_nums)) # pick a perturbation level
      {
        
        vr2[[g]] <- vr_nums[i] # specify the vital rate with the new perturbation level
        
        A <- matrix(sapply(matrix_str, eval, vr2, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        I <- matrix(data = c(0, vital_rate_summary$F_Adult_immigration, 
                             0, vital_rate_summary$M_Adult_immigration), 
                    nrow = 4, ncol = 1)
        
        # reset the starting stage distribution for simulation (all with 10 individuals)
        m <- rep(10, no_stages) 
        
        for (j in 1:niter) { # project the matrix through t iteration
          # stage distribution at time t
          stage[,j] <- m
          
          # population size at time t
          pop[j] <- sum(m)
          
          # number of male adults at time t
          M2 <- stage[4, i] * A["M_Adult", "M_Adult"] + I[4, 1]
          
          # number of female adults at time t
          F2 <- stage[2, i] * A["F_Adult", "F_Adult"] + I[2, 1]
          
          # Female freq-dep fecundity of Female chicks
          A["F_1st_year", "F_Adult"] <- ((k * M2) / (M2 + (F2 / h)) * (1 - HSR) )
          
          # Female freq-dep fecundity of Male chicks
          A["M_1st_year", "F_Adult"] <- ((k * M2) / (M2 + (F2 / h)) * HSR)
          
          # Male freq-dep fecundity of Female chicks
          A["F_1st_year", "M_Adult"] <- ((k * F2) / (M2 + (F2 / h)) * (1 - HSR))
          
          # Male freq-dep fecundity of Male chicks
          A["M_1st_year", "M_Adult"] <- ((k * F2) / (M2 + (F2 / h)) * HSR)
          
          # define the new n (i.e., new stage distribution at time t)
          m <- A %*% (m + I)
        }
        
        # define rownames of stage matrix
        rownames(stage) <- rownames(A)
        
        # define colnames of stage matrix
        colnames(stage) <- 0:(niter - 1)
        
        # calculate the proportional stable stage distribution
        stage <- apply(stage, 2, function(x) x / sum(x))
        
        # define stable stage as the last stage
        stable.stage <- stage[, niter]
        
        # calc ASR as the proportion of the adult stable stage class that is male
        vr_pert_ASR[i, g] <- stable.stage[no_stages] / (stable.stage[no_stages/2] + 
                                                          stable.stage[no_stages])
        
        # calc lambda as the pop change in the counts of the last two iterations
        vr_pert_lambda[i, g] <- pop[niter]/pop[niter - 1]
      }
      
      # get the spline function of ASR
      spl_ASR <- 
        smooth.spline(vr_pert_ASR[which(!is.na(vr_pert_ASR[, g])), g] ~ 
                        names(vr_pert_ASR[which(!is.na(vr_pert_ASR[, g])), g]))
      
      # estimate the slope of the tangent of the spline at the vital rate
      ASR_pert_results[g, 2] <- predict(spl_ASR, x = vr[[g]], deriv = 1)$y
      
      # re-scale sensitivity into elasticity
      ASR_pert_results[g, 3] <- vr[[g]] / ASR * ASR_pert_results[g, 2]
      
      # do the same steps but for lambda
      spl_lambda <- 
        smooth.spline(vr_pert_lambda[which(!is.na(vr_pert_lambda[, g])), g] ~ 
                        names(vr_pert_lambda[which(!is.na(vr_pert_lambda[, g])), g]))
      
      lambda_pert_results[g, 2] <- predict(spl_lambda, x = vr[[g]], deriv = 1)$y
      
      lambda_pert_results[g, 3] <- vr[[g]] / lambda * lambda_pert_results[g, 2]
    }
    
    ##### perturbation of the h index parameter ####
    for (i in 1:length(h_nums)) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                  nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                  dimnames = list(stages, stages)) 
      
      I <- matrix(data = c(0, vital_rate_summary$F_Adult_immigration, 
                           0, vital_rate_summary$M_Adult_immigration), 
                  nrow = 4, ncol = 1)
      
      # reset the starting stage distribution for simulation (all with 10 individuals)
      m <- rep(10, no_stages) 
      
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m
        
        # population size at time t
        pop[j] <- sum(m)
        
        # number of male adults at time t
        M2 <- stage[4, i] * A["M_Adult", "M_Adult"] + I[4, 1]
        
        # number of female adults at time t
        F2 <- stage[2, i] * A["F_Adult", "F_Adult"] + I[2, 1]
        
        # Female freq-dep fecundity of Female chicks
        A["F_1st_year", "F_Adult"] <- ((k * M2) / (M2 + (F2 / h_nums[i])) * (1 - HSR) )
        
        # Female freq-dep fecundity of Male chicks
        A["M_1st_year", "F_Adult"] <- ((k * M2) / (M2 + (F2 / h_nums[i])) * HSR)
        
        # Male freq-dep fecundity of Female chicks
        A["F_1st_year", "M_Adult"] <- ((k * F2) / (M2 + (F2 / h_nums[i])) * (1 - HSR))
        
        # Male freq-dep fecundity of Male chicks
        A["M_1st_year", "M_Adult"] <- ((k * F2) / (M2 + (F2 / h_nums[i])) * HSR)
        
        # define the new n (i.e., new stage distribution at time t)
        m <- A %*% (m + I)
      }
      # define rownames of stage matrix
      rownames(stage) <- rownames(A)
      
      # define colnames of stage matrix
      colnames(stage) <- 0:(niter - 1)
      
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x / sum(x))
      
      # define stable stage as the last stage
      stable.stage <- stage[, niter]
      
      # calc ASR as the proportion of the adult stable stage class that is male
      h_pert_ASR[i,] <- stable.stage[no_stages] / (stable.stage[no_stages / 2] + 
                                                     stable.stage[no_stages])
      
      # calc lambda as the pop change in the counts of the last two iterations
      h_pert_lambda[i, ] <- pop[niter]/pop[niter - 1]
      
    }
    # get the spline function of ASR
    spl_ASR <- 
      smooth.spline(h_pert_ASR[which(!is.na(h_pert_ASR)), 1] ~ 
                      names(h_pert_ASR[which(!is.na(h_pert_ASR)), ]))
    
    # estimate the slope of the tangent of the spline at the vital rate
    ASR_pert_results[n + 1, 2] <- predict(spl_ASR, x = h, deriv = 1)$y
    
    # re-scale sensitivity into elasticity
    ASR_pert_results[n + 1, 3] <- h / ASR * ASR_pert_results[n + 1, 2]
    
    # do the same steps but for lambda
    spl_lambda <- 
      smooth.spline(h_pert_lambda[which(!is.na(h_pert_lambda)), 1] ~ 
                      names(h_pert_lambda[which(!is.na(h_pert_lambda)), ]))
    lambda_pert_results[n + 1, 2] <- predict(spl_lambda, x = h, deriv = 1)$y
    lambda_pert_results[n + 1, 3] <- h / lambda * lambda_pert_results[n + 1, 2]
    
    ##### perturbation of k parameter ####
    for (i in 1:length(k_nums)) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                  nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                  dimnames = list(stages, stages))
      
      I <- matrix(data = c(0, vital_rate_summary$F_Adult_immigration, 
                           0, vital_rate_summary$M_Adult_immigration), 
                  nrow = 4, ncol = 1)
      
      # reset the starting stage distribution for simulation (all with 10 individuals)
      m <- rep(10, no_stages) 
      
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m
        
        # population size at time t
        pop[j] <- sum(m)
        
        # number of male adults at time t
        M2 <- stage[4, i] * A["M_Adult", "M_Adult"] + I[4, 1]
        
        # number of female adults at time t
        F2 <- stage[2, i] * A["F_Adult", "F_Adult"] + I[2, 1]
        
        # Female freq-dep fecundity of Female chicks
        A["F_1st_year", "F_Adult"] <- ((k_nums[i] * M2) / (M2 + (F2 / h)) * (1 - HSR) )
        
        # Female freq-dep fecundity of Male chicks
        A["M_1st_year", "F_Adult"] <- ((k_nums[i] * M2) / (M2 + (F2 / h)) * HSR)
        
        # Male freq-dep fecundity of Female chicks
        A["F_1st_year", "M_Adult"] <- ((k_nums[i] * F2) / (M2 + (F2 / h)) * (1 - HSR))
        
        # Male freq-dep fecundity of Male chicks
        A["M_1st_year", "M_Adult"] <- ((k_nums[i] * F2) / (M2 + (F2 / h)) * HSR)
        
        # define the new n (i.e., new stage distribution at time t)
        m <- A %*% (m + I)
      }
      # define rownames of stage matrix
      rownames(stage) <- rownames(A)
      
      # define colnames of stage matrix
      colnames(stage) <- 0:(niter - 1)
      
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x / sum(x))
      
      # define stable stage as the last stage
      stable.stage <- stage[, niter]
      
      # calc ASR as the proportion of the adult stable stage class that is male
      k_pert_ASR[i,] <- stable.stage[no_stages] / (stable.stage[no_stages/2] + 
                                                     stable.stage[no_stages])
      
      # calc lambda as the pop change in the counts of the last two iterations
      k_pert_lambda[i, ] <- pop[niter]/pop[niter - 1]
    }
    
    # get the spline function of ASR
    spl_ASR <- 
      smooth.spline(k_pert_ASR[which(!is.na(k_pert_ASR)), 1] ~ 
                      names(k_pert_ASR[which(!is.na(k_pert_ASR)), ]))
    
    # estimate the slope of the tangent of the spline at the vital rate
    
    ASR_pert_results[n + 2, 2] <- predict(spl_ASR, x = k, deriv = 1)$y
    
    # re-scale sensitivity into elasticity
    ASR_pert_results[n + 2, 3] <- k / ASR * ASR_pert_results[n + 2, 2]
    
    # do the same steps but for lambda
    spl_lambda <- 
      smooth.spline(h_pert_lambda[which(!is.na(h_pert_lambda)), 1] ~ 
                      names(h_pert_lambda[which(!is.na(h_pert_lambda)), ]))
    lambda_pert_results[n + 2, 2] <- predict(spl_lambda, x = k, deriv = 1)$y
    lambda_pert_results[n + 2, 3] <- k / lambda * lambda_pert_results[n + 2, 2]
    
    ##### perturbation of HSR ####
    for (i in 1:length(vr_nums)) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                  nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                  dimnames = list(stages, stages))
      
      I <- matrix(data = c(0, vital_rate_summary$F_Adult_immigration, 
                           0, vital_rate_summary$M_Adult_immigration), 
                  nrow = 4, ncol = 1)
      
      # reset the starting stage distribution for simulation (all with 10 individuals)
      m <- rep(10, no_stages) 
      
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m
        
        # population size at time t
        pop[j] <- sum(m)
        
        # number of male adults at time t
        M2 <- stage[4, i] * A["M_Adult", "M_Adult"] + I[4, 1]
        
        # number of female adults at time t
        F2 <- stage[2, i] * A["F_Adult", "F_Adult"] + I[2, 1]
        
        # Female freq-dep fecundity of Female chicks
        A["F_1st_year", "F_Adult"] <- ((k * M2) / (M2 + (F2 / h)) * (1 - vr_nums[i]) )
        
        # Female freq-dep fecundity of Male chicks
        A["M_1st_year", "F_Adult"] <- ((k * M2) / (M2 + (F2 / h)) * vr_nums[i])
        
        # Male freq-dep fecundity of Female chicks
        A["F_1st_year", "M_Adult"] <- ((k * F2) / (M2 + (F2 / h)) * (1 - vr_nums[i]))
        
        # Male freq-dep fecundity of Male chicks
        A["M_1st_year", "M_Adult"] <- ((k * F2) / (M2 + (F2 / h)) * vr_nums[i])
        
        # define the new n (i.e., new stage distribution at time t)
        m <- A %*% (m + I)
        
      }
      # define rownames of stage matrix
      rownames(stage) <- rownames(A)
      
      # define colnames of stage matrix
      colnames(stage) <- 0:(niter - 1)
      
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x / sum(x))
      
      # define stable stage as the last stage
      stable.stage <- stage[, niter]
      
      # calc ASR as the proportion of the adult stable stage class that is male
      HSR_pert_ASR[i,] <- stable.stage[no_stages] / (stable.stage[no_stages/2] + 
                                                       stable.stage[no_stages])
      
      # calc lambda as the pop change in the counts of the last two iterations
      HSR_pert_lambda[i, ] <- pop[niter] / pop[niter - 1]
      
    }
    # get the spline function of ASR
    spl_ASR <- 
      smooth.spline(HSR_pert_ASR[which(!is.na(HSR_pert_ASR)), 1] ~ 
                      names(HSR_pert_ASR[which(!is.na(HSR_pert_ASR)), ]))
    
    # estimate the slope of the tangent of the spline at the vital rate    
    ASR_pert_results[n + 3, 2] <- predict(spl_ASR, x = HSR, deriv = 1)$y
    
    # re-scale sensitivity into elasticity
    ASR_pert_results[n + 3, 3] <- HSR / ASR * ASR_pert_results[n + 3, 2]
    
    # do the same steps but for lambda
    spl_lambda <- 
      smooth.spline(h_pert_lambda[which(!is.na(h_pert_lambda)), 1] ~ 
                      names(h_pert_lambda[which(!is.na(h_pert_lambda)), ]))
    lambda_pert_results[n + 3, 2] <- predict(spl_lambda, x = HSR, deriv = 1)$y
    lambda_pert_results[n + 3, 3] <- HSR/lambda * lambda_pert_results[n + 3, 2]
    
    ##### perturbation of ISR ####
    for (i in 1:length(vr_nums)) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                  nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                  dimnames = list(stages, stages))
      
      I <- matrix(data = c(0, immigrant_pop_size * (1 - vr_nums[i]), 
                           0, immigrant_pop_size * vr_nums[i]), 
                  nrow = 4, ncol = 1)
      
      # reset the starting stage distribution for simulation (all with 10 individuals)
      m <- rep(10, no_stages) 
      
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m
        
        # population size at time t
        pop[j] <- sum(m)
        
        # number of male adults at time t
        M2 <- stage[4, i] * A["M_Adult", "M_Adult"] + I[4, 1]
        
        # number of female adults at time t
        F2 <- stage[2, i] * A["F_Adult", "F_Adult"] + I[2, 1]
        
        # Female freq-dep fecundity of Female chicks
        A["F_1st_year", "F_Adult"] <- ((k * M2) / (M2 + (F2 / h)) * (1 - HSR) )
        
        # Female freq-dep fecundity of Male chicks
        A["M_1st_year", "F_Adult"] <- ((k * M2) / (M2 + (F2 / h)) * HSR)
        
        # Male freq-dep fecundity of Female chicks
        A["F_1st_year", "M_Adult"] <- ((k * F2) / (M2 + (F2 / h)) * (1 - HSR))
        
        # Male freq-dep fecundity of Male chicks
        A["M_1st_year", "M_Adult"] <- ((k * F2) / (M2 + (F2 / h)) * HSR)
        
        # define the new n (i.e., new stage distribution at time t)
        m <- A %*% (m + I)
        
      }
      # define rownames of stage matrix
      rownames(stage) <- rownames(A)
      
      # define colnames of stage matrix
      colnames(stage) <- 0:(niter - 1)
      
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x / sum(x))
      
      # define stable stage as the last stage
      stable.stage <- stage[, niter]
      
      # calc ASR as the proportion of the adult stable stage class that is male
      ISR_pert_ASR[i,] <- stable.stage[no_stages] / (stable.stage[no_stages/2] + 
                                                       stable.stage[no_stages])
      
      # calc lambda as the pop change in the counts of the last two iterations
      ISR_pert_lambda[i, ] <- pop[niter] / pop[niter - 1]
      
    }
    # get the spline function of ASR
    spl_ASR <- 
      smooth.spline(ISR_pert_ASR[which(!is.na(ISR_pert_ASR)), 1] ~ 
                      names(ISR_pert_ASR[which(!is.na(ISR_pert_ASR)), ]))
    
    # estimate the slope of the tangent of the spline at the vital rate    
    ASR_pert_results[n + 4, 2] <- predict(spl_ASR, x = ISR, deriv = 1)$y
    
    # re-scale sensitivity into elasticity
    ASR_pert_results[n + 4, 3] <- ISR / ASR * ASR_pert_results[n + 4, 2]
    
    # do the same steps but for lambda
    spl_lambda <- 
      smooth.spline(ISR_pert_lambda[which(!is.na(ISR_pert_lambda)), 1] ~ 
                      names(ISR_pert_lambda[which(!is.na(ISR_pert_lambda)), ]))
    lambda_pert_results[n + 4, 2] <- predict(spl_lambda, x = ISR, deriv = 1)$y
    lambda_pert_results[n + 4, 3] <- ISR / lambda * lambda_pert_results[n + 4, 2]
    
    
    
    
    #### store all results into a list ----
    result <- list(ASR_pert_results = ASR_pert_results,
                   lambda_pert_results = lambda_pert_results,
                   ISR_pert_ASR = ISR_pert_ASR,
                   ISR_pert_lambda = ISR_pert_lambda,
                   HSR_pert_ASR = HSR_pert_ASR,
                   HSR_pert_lambda = HSR_pert_lambda,
                   k_pert_ASR = k_pert_ASR,
                   k_pert_lambda = k_pert_lambda,
                   h_pert_ASR = h_pert_ASR,
                   h_pert_lambda = h_pert_lambda,
                   vr_pert_ASR = vr_pert_ASR,
                   vr_pert_lambda = vr_pert_lambda)
    
  }