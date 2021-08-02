source("R/project/project_libraries.R")

# make_treat_matrix() makes a treatment matrix from a summary or the bootstrapped
# vital rates

make_treat_matrix <- 
  function(survival_rates_boot_summary, species_name, h, k, HSR, ISR){
    
    list(F_Nestling_survival = filter(survival_rates_boot_summary,
                                      species == species_name)[7, 4],
         F_Groundling_survival = filter(survival_rates_boot_summary,
                                        species == species_name)[6, 4],
         F_Fledgling_survival = filter(survival_rates_boot_summary,
                                       species == species_name)[4, 4],
         F_Adult_survival = filter(survival_rates_boot_summary,
                                   species == species_name)[2, 4],
         M_Nestling_survival = filter(survival_rates_boot_summary,
                                      species == species_name)[14, 4],
         M_Groundling_survival = filter(survival_rates_boot_summary,
                                        species == species_name)[13, 4],
         M_Fledgling_survival = filter(survival_rates_boot_summary,
                                       species == species_name)[11, 4],
         M_Adult_survival = filter(survival_rates_boot_summary,
                                   species == species_name)[9, 4],
         Egg_survival = filter(survival_rates_boot_summary,
                               species == species_name)[15, 4],
         
         # Define h (harem size, h = 1 is monogamy) and k (clutch size)
         h = h,
         k = k,
         
         # Define primary sex ratio
         HSR = HSR,
         
         # Define immigrant sex ratio
         ISR = ISR)
  }
