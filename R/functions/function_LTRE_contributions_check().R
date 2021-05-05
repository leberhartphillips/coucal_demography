source("scripts/01_libraries.R")

# LTRE_contributions_check() takes the results of the LTRE analysis and checks
# the total contribution that each vital rate has on ASR. In general, the 
# sum of the LTRE contribution should roughly equal the difference in ASR between
# the observed scenario and the unbiased scenario.

LTRE_contributions_check <- 
  function(M_matrix_vital_rates, Mprime_sensitivities, M_matrix_ASR, scenario){
    if (scenario == "male")
      contribution_sum <- 
        sum(
          (M_matrix_vital_rates[[1]] - M_matrix_vital_rates[[5]]) * 
            Mprime_sensitivities[1],
          
          (M_matrix_vital_rates[[2]] - M_matrix_vital_rates[[6]]) * 
            Mprime_sensitivities[2],
          
          (M_matrix_vital_rates[[3]] - M_matrix_vital_rates[[7]]) * 
            Mprime_sensitivities[3],
          
          (M_matrix_vital_rates[[4]] - M_matrix_vital_rates[[8]]) * 
            Mprime_sensitivities[4],
          
          (M_matrix_vital_rates[[5]] - M_matrix_vital_rates[[5]]) * 
            Mprime_sensitivities[1],
          
          (M_matrix_vital_rates[[6]] - M_matrix_vital_rates[[6]]) * 
            Mprime_sensitivities[2],
          
          (M_matrix_vital_rates[[7]] - M_matrix_vital_rates[[7]]) * 
            Mprime_sensitivities[3],
          
          (M_matrix_vital_rates[[8]] - M_matrix_vital_rates[[8]]) * 
            Mprime_sensitivities[4],
          
          (M_matrix_vital_rates[[10]] - 1) *
            Mprime_sensitivities[9],
          
          (M_matrix_vital_rates[[12]] - 0.5) * Mprime_sensitivities[11],
          
          (M_matrix_vital_rates[[13]] - 0.5) * Mprime_sensitivities[12]
          
        )
    
    else
      contribution_sum <- 
        sum(
          (M_matrix_vital_rates[[5]] - M_matrix_vital_rates[[1]]) * 
            Mprime_sensitivities[5],
          
          (M_matrix_vital_rates[[6]] - M_matrix_vital_rates[[2]]) * 
            Mprime_sensitivities[6],
          
          (M_matrix_vital_rates[[7]] - M_matrix_vital_rates[[3]]) * 
            Mprime_sensitivities[7],
          
          (M_matrix_vital_rates[[8]] - M_matrix_vital_rates[[4]]) * 
            Mprime_sensitivities[8],
          
          (M_matrix_vital_rates[[1]] - M_matrix_vital_rates[[1]]) * 
            Mprime_sensitivities[5],
          
          (M_matrix_vital_rates[[2]] - M_matrix_vital_rates[[2]]) * 
            Mprime_sensitivities[6],
          
          (M_matrix_vital_rates[[3]] - M_matrix_vital_rates[[3]]) * 
            Mprime_sensitivities[7],
          
          (M_matrix_vital_rates[[4]] - M_matrix_vital_rates[[4]]) * 
            Mprime_sensitivities[8],
          
          (M_matrix_vital_rates[[10]] - 1) *
            Mprime_sensitivities[9],
          
          (M_matrix_vital_rates[[12]] - 0.5) * Mprime_sensitivities[11],
          
          (M_matrix_vital_rates[[13]] - 0.5) * Mprime_sensitivities[12]
        )
    ASR_bias <- abs(M_matrix_ASR - 0.5)
    absolute_difference <- abs(ASR_bias) - abs(contribution_sum)
    
    return(list(contribution_sum = as.vector(contribution_sum), 
                ASR_bias = as.vector(ASR_bias), 
                absolute_difference = as.vector(absolute_difference)))
  }