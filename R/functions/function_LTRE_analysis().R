source("R/project/project_libraries.R")

# LTRE_analysis() compares the relative difference in ASR response of an
# unbiased two-sex matrix vs. the observed sex-specific matrix

LTRE_analysis <-
  function(Mprime_sens, matrix_str, vital_rates, species_name, sex){
    
    # make empty dataframes to store LTRE results for ASR and lambda
    LTRE_ASR <-
      data.frame(parameter = c("Nestling survival", "Groundling survival",
                               "Fledgling survival", "Adult survival",
                               "Hatching sex ratio", "Immigrant sex ratio",
                               "Mating system"),
                 contribution = numeric(7))
    LTRE_lambda <-
      data.frame(parameter = c("Nestling survival", "Groundling survival",
                               "Fledgling survival", "Adult survival",
                               "Hatching sex ratio", "Immigrant sex ratio",
                               "Mating system"),
                 contribution = numeric(7))
    
    # run a for loop to extract the parameter contributions
    if (sex == "male"){
      # male rates scenario
      for(i in 1:nrow(LTRE_ASR))
      {
        LTRE_ASR[i, 2] <-
          
          # survival rates
          ifelse(i < 5, (vital_rates[[i + 4]] - vital_rates[[i]]) * 
                   Mprime_sens$ASR_pert_results$sensitivities[i + 4],
                 
                 # HSR
                 ifelse(i == 5, (vital_rates[[12]] - 0.5) * 
                          Mprime_sens$ASR_pert_results$sensitivities[11],
                        
                        # ISR
                        ifelse(i == 6, (vital_rates[[13]] - 0.5) * 
                                 Mprime_sens$ASR_pert_results$sensitivities[12],
                               
                               # mating system
                               (1 - vital_rates[[10]]) * 
                                 Mprime_sens$ASR_pert_results$sensitivities[9])))
      }
      for(i in 1:nrow(LTRE_lambda))
      {
        # survival rates
        ifelse(i < 5, (vital_rates[[i + 4]] - vital_rates[[i]]) * 
                 Mprime_sens$lambda_pert_results$sensitivities[i + 4],
               
               # HSR
               ifelse(i == 5, (vital_rates[[12]] - 0.5) * 
                        Mprime_sens$lambda_pert_results$sensitivities[11],
                      
                      # ISR
                      ifelse(i == 6, (vital_rates[[13]] - 0.5) * 
                               Mprime_sens$lambda_pert_results$sensitivities[12],
                             
                             # mating system
                             (1 - vital_rates[[10]]) * 
                               Mprime_sens$lambda_pert_results$sensitivities[9])))
      }
    }
    else{
      # female rates scenario
      for(i in 1:nrow(LTRE_ASR))
      {
        LTRE_ASR[i, 2] <-
          
          # survival rates
          ifelse(i < 5, (vital_rates[[i]] - vital_rates[[i + 4]]) * 
                   Mprime_sens$ASR_pert_results$sensitivities[i],
                 
                 # HSR
                 ifelse(i == 5, (vital_rates[[12]] - 0.5) * 
                          Mprime_sens$ASR_pert_results$sensitivities[11],
                        
                        # ISR
                        ifelse(i == 6, (vital_rates[[13]] - 0.5) * 
                                 Mprime_sens$ASR_pert_results$sensitivities[12],
                               
                               # mating system
                               (1 - vital_rates[[10]]) * 
                                 Mprime_sens$ASR_pert_results$sensitivities[9])))
      }
      
      for(i in 1:nrow(LTRE_lambda))
      {
        # survival rates
        ifelse(i < 5, (vital_rates[[i]] - vital_rates[[i + 4]]) * 
                 Mprime_sens$lambda_pert_results$sensitivities[i],
               
               # HSR
               ifelse(i == 5, (vital_rates[[12]] - 0.5) * 
                        Mprime_sens$lambda_pert_results$sensitivities[11],
                      
                      # ISR
                      ifelse(i == 6, (vital_rates[[13]] - 0.5) * 
                               Mprime_sens$lambda_pert_results$sensitivities[12],
                             
                             # mating system
                             (1 - vital_rates[[10]]) * 
                               Mprime_sens$lambda_pert_results$sensitivities[9])))
      }
    }
    LTRE_ASR$parameter <- 
      factor(LTRE_ASR$parameter, levels = c("Adult survival",
                                            "Fledgling survival",
                                            "Groundling survival",
                                            "Nestling survival",
                                            "Hatching sex ratio",
                                            "Immigrant sex ratio",
                                            "Mating system"))
    LTRE_lambda$parameter <- 
      factor(LTRE_lambda$parameter, levels = c("Adult survival",
                                               "Fledgling survival",
                                               "Groundling survival",
                                               "Nestling survival",
                                               "Hatching sex ratio",
                                               "Immigrant sex ratio",
                                               "Mating system"))
    
    LTRE_ASR$model <- "ASR"
    LTRE_ASR$species <- species_name
    
    LTRE_lambda$model <- "lambda"
    LTRE_lambda$species <- species_name
    
    LTRE_results <- list(LTRE_ASR = LTRE_ASR,
                         LTRE_lambda = LTRE_lambda)
    
    LTRE_results
  }
