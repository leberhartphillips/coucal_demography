source("R/project/project_libraries.R")

# coucal_matrix() builds the two-sex Lefkovitch matrix using the vital rates 
# specified in the *demographic_rates* object.

coucal_matrix <- 
  function(demographic_rates, two_sex = TRUE){
    if(two_sex){
      
      # Define coucal life-stages of the coucal matrix model
      stages <- c("F_1st_year",  "F_Adult",  "M_1st_year",  "M_Adult")
      
      # Build the 4x4 matrix
      result <- 
        matrix(c(
          
          # top row of matrix
          0, NA, 0, NA, 
          
          # second row of matrix
          (demographic_rates$Egg_survival * 
             demographic_rates$F_Nestling_survival * 
             demographic_rates$F_Groundling_survival *
             demographic_rates$F_Fledgling_survival),
          demographic_rates$F_Adult_survival,
          0, 0,
          
          # third row of matrix
          0, NA, 0, NA, 
          
          # fourth row of matrix
          0, 0, 
          (demographic_rates$Egg_survival * 
             demographic_rates$M_Nestling_survival * 
             demographic_rates$M_Groundling_survival *
             demographic_rates$M_Fledgling_survival),
          demographic_rates$M_Adult_survival),
          nrow = length(stages), byrow = TRUE,
          dimnames = list(stages, stages))
    }
    else{
      
      # Define coucal life-stages of the Ceuta snowy coucal matrix model
      stages <- c("1st_year",  "Adult")
      
      # Build the 4x4 matrix
      result <- 
        matrix(c(
          
          # top row of matrix
          0, NA,
          
          # second row of matrix
          (demographic_rates$Egg_survival * 
             demographic_rates$Nestling_survival * 
             demographic_rates$Groundling_survival *
             demographic_rates$Fledgling_survival),
          demographic_rates$Adult_survival),
          nrow = length(stages), byrow = TRUE,
          dimnames = list(stages, stages))
    }
    result
  }