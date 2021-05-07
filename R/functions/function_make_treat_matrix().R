source("R/project/project_libraries.R")

# make_treat_matrix() makes a treatment matrix from a summary or the bootstrapped
# vital rates

make_treat_matrix <- 
  function(survival_rates_boot_summary, species, h, k, HSR, ISR){
    
    list(F_Nestling_survival = filter(survival_rates_boot_summary,
                                      species == species)[5, 4],
         F_Groundling_survival = filter(survival_rates_boot_summary,
                                        species == species)[4, 4],
         F_Fledgling_survival = filter(survival_rates_boot_summary,
                                       species == species)[3, 4],
         F_Adult_survival = filter(survival_rates_boot_summary,
                                   species == species)[2, 4],
         M_Nestling_survival = filter(survival_rates_boot_summary,
                                      species == species)[10, 4],
         M_Groundling_survival = filter(survival_rates_boot_summary,
                                        species == species)[9, 4],
         M_Fledgling_survival = filter(survival_rates_boot_summary,
                                       species == species)[8, 4],
         M_Adult_survival = filter(survival_rates_boot_summary,
                                   species == species)[7, 4],
         Egg_survival = filter(survival_rates_boot_summary,
                               species == species)[11, 4],
         
         # Define h (harem size, h = 1 is monogamy) and k (clutch size)
         h = h,
         k = k,
         
         # Define primary sex ratio
         HSR = HSR,
         
         # Define immigrant sex ratio
         ISR = ISR)
  }