source("R/project/project_libraries.R")

# LTRE_analysis() compares the relative difference in ASR response of an
# unbiased two-sex matrix vs. the observed sex-specific matrix

LTRE_analysis_JSR <-
  function(Mprime_sens, VR_dataframe, base_sex){
    if(base_sex == "male"){
    # make empty dataframes to store LTRE results for ASR and lambda
    LTRE_JSR <-
      VR_dataframe %>% 
      mutate(parameter = ifelse(str_detect(parameter, "Female"), str_replace(parameter, "Female", "F"),
                                ifelse(str_detect(parameter, "Male"), str_replace(parameter, "Male", "M"), parameter))) %>% 
      left_join(., Mprime_sens$JSR_pert_results, by = "parameter") %>% 
      mutate(sex = ifelse(str_detect(parameter, "F"), "F", 
                          ifelse(str_detect(parameter, "M"), "M", NA))) %>% 
      dplyr::select(!elasticities) %>% 
      mutate(parameter = ifelse(str_detect(parameter, "F_age_"), str_remove(parameter, "F_"), 
                                ifelse(str_detect(parameter, "M_age_"), str_remove(parameter, "M_"), parameter))) %>% 
      pivot_wider(names_from = sex, values_from = c(obs, sensitivities)) %>% 
      mutate(contribution = ifelse(str_detect(parameter, "age_"), (obs_F - obs_M) * sensitivities_F,
                                     ifelse(parameter == "h", (1 - obs_NA) * sensitivities_NA,
                                                   ifelse(parameter == "HSR", (obs_NA - 0.5) * sensitivities_NA,
                                                          ifelse(parameter == "ISR", (obs_NA - 0.5) * sensitivities_NA, NA))))) %>% 
      filter(!is.na(contribution))
    }
    else{
      LTRE_JSR <-
        VR_dataframe %>% 
        mutate(parameter = ifelse(str_detect(parameter, "Female"), str_replace(parameter, "Female", "F"),
                                  ifelse(str_detect(parameter, "Male"), str_replace(parameter, "Male", "M"), parameter))) %>% 
        left_join(., Mprime_sens$JSR_pert_results, by = "parameter") %>% 
        mutate(sex = ifelse(str_detect(parameter, "F"), "F", 
                            ifelse(str_detect(parameter, "M"), "M", NA))) %>% 
        dplyr::select(!elasticities) %>% 
        mutate(parameter = ifelse(str_detect(parameter, "F_age_"), str_remove(parameter, "F_"), 
                                  ifelse(str_detect(parameter, "M_age_"), str_remove(parameter, "M_"), parameter))) %>% 
        pivot_wider(names_from = sex, values_from = c(obs, sensitivities)) %>% 
        mutate(contribution = ifelse(str_detect(parameter, "age_"), (obs_M - obs_F) * sensitivities_M,
                                     ifelse(parameter == "h", (1 - obs_NA) * sensitivities_NA,
                                            ifelse(parameter == "HSR", (obs_NA - 0.5) * sensitivities_NA,
                                                   ifelse(parameter == "ISR", (obs_NA - 0.5) * sensitivities_NA, NA))))) %>% 
        filter(!is.na(contribution))
    }
    
    LTRE_JSR
  }
