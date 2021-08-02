source("R/project/project_libraries.R")

# make_mprime_matrix() makes a prime-matrix (i.e., a matrix halfway between
# the treatment matrix and the unbiased matrix) from a set of vital rate
# summaries

make_mprime_matrix <- 
  function(survival_rates_boot_summary, species, h, k, HSR, ISR, sex){
    if(sex == "male"){
      list(F_Nestling_survival = (filter(survival_rates_boot_summary,
                                         species == species)[7, 4] + 
                                    filter(survival_rates_boot_summary,
                                           species == species)[14, 4]) / 2,
           F_Groundling_survival = (filter(survival_rates_boot_summary,
                                           species == species)[6, 4] + 
                                      filter(survival_rates_boot_summary,
                                             species == species)[13, 4]) / 2,
           F_Fledgling_survival = (filter(survival_rates_boot_summary,
                                          species == species)[4, 4] + 
                                     filter(survival_rates_boot_summary,
                                            species == species)[11, 4]) / 2,
           F_Adult_survival = (filter(survival_rates_boot_summary,
                                      species == species)[2, 4] + 
                                 filter(survival_rates_boot_summary,
                                        species == species)[9, 4]) / 2,
           M_Nestling_survival = (filter(survival_rates_boot_summary,
                                         species == species)[14, 4] + 
                                    filter(survival_rates_boot_summary,
                                           species == species)[14, 4]) / 2,
           M_Groundling_survival = (filter(survival_rates_boot_summary,
                                           species == species)[13, 4] + 
                                      filter(survival_rates_boot_summary,
                                             species == species)[13, 4]) / 2,
           M_Fledgling_survival = (filter(survival_rates_boot_summary,
                                          species == species)[11, 4] + 
                                     filter(survival_rates_boot_summary,
                                            species == species)[11, 4]) / 2,
           M_Adult_survival = (filter(survival_rates_boot_summary,
                                      species == species)[9, 4] + 
                                 filter(survival_rates_boot_summary,
                                        species == species)[9, 4]) / 2,
           Egg_survival = filter(survival_rates_boot_summary,
                                 species == species)[15, 4],
           
           # Define h (harem size, h = 1 is monogamy) and k (clutch size)
           h = (h + 1) / 2,
           k = k,
           
           # Define primary sex ratio
           HSR = (HSR + 0.5) / 2,
           
           # Define immigrant sex ratio
           ISR = (ISR + 0.5) / 2)
    }
    else{
      list(F_Nestling_survival = (filter(survival_rates_boot_summary,
                                         species == species)[7, 4] + 
                                    filter(survival_rates_boot_summary,
                                           species == species)[7, 4]) / 2,
           F_Groundling_survival = (filter(survival_rates_boot_summary,
                                           species == species)[6, 4] + 
                                      filter(survival_rates_boot_summary,
                                             species == species)[6, 4]) / 2,
           F_Fledgling_survival = (filter(survival_rates_boot_summary,
                                          species == species)[4, 4] + 
                                     filter(survival_rates_boot_summary,
                                            species == species)[4, 4]) / 2,
           F_Adult_survival = (filter(survival_rates_boot_summary,
                                      species == species)[2, 4] + 
                                 filter(survival_rates_boot_summary,
                                        species == species)[2, 4]) / 2,
           M_Nestling_survival = (filter(survival_rates_boot_summary,
                                         species == species)[14, 4] + 
                                    filter(survival_rates_boot_summary,
                                           species == species)[7, 4]) / 2,
           M_Groundling_survival = (filter(survival_rates_boot_summary,
                                           species == species)[13, 4] + 
                                      filter(survival_rates_boot_summary,
                                             species == species)[6, 4]) / 2,
           M_Fledgling_survival = (filter(survival_rates_boot_summary,
                                          species == species)[11, 4] + 
                                     filter(survival_rates_boot_summary,
                                            species == species)[4, 4]) / 2,
           M_Adult_survival = (filter(survival_rates_boot_summary,
                                      species == species)[9, 4] + 
                                 filter(survival_rates_boot_summary,
                                        species == species)[2, 4]) / 2,
           Egg_survival = filter(survival_rates_boot_summary,
                                 species == species)[15, 4],
           
           # Define h (harem size, h = 1 is monogamy) and k (clutch size)
           h = (h + 1) / 2,
           k = k,
           
           # Define primary sex ratio
           HSR = (HSR + 0.5) / 2,
           
           # Define immigrant sex ratio
           ISR = (ISR + 0.5) / 2)
    }
  }
