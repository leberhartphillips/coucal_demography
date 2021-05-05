source("scripts/01_libraries.R")

# AIC_table_wrangle() summarizes all the AIC-based model selection results
# produced by the bootstrapping procedure. It is run to prepare the data
# for visualization

AIC_table_wrangle <- function(boot_out_list){
  nestling_AIC_tables_boot_summary <- 
    boot_out_list$nestling_AIC_tables_boot %>%
    dplyr::group_by(model, model_no_orig) %>%
    dplyr::summarise(avg_Delta = mean(DeltaAICc),
                     IQR_Delta = IQR(DeltaAICc),
                     avg_Weight = mean(weight),
                     IQR_Weight = IQR(weight), .groups = 'drop') %>% 
    dplyr::arrange(avg_Delta) %>% 
    mutate(#model_number = as.numeric(model),
           stage = "nestling") %>% 
    dplyr::left_join(., select(boot_out_list$nestling_AIC_tables_boot, -model), by = "model_no_orig")
  
  fledgling_AIC_tables_boot_summary <- 
    boot_out_list$fledgling_AIC_tables_boot %>%
    dplyr::group_by(model, model_no_orig) %>%
    dplyr::summarise(avg_Delta = mean(DeltaAICc),
                     IQR_Delta = IQR(DeltaAICc),
                     avg_Weight = mean(weight),
                     IQR_Weight = IQR(weight), .groups = 'drop') %>% 
    dplyr::arrange(avg_Delta) %>% 
    mutate(#model_number = as.numeric(model),
           stage = "fledgling") %>% 
    dplyr::left_join(., select(boot_out_list$fledgling_AIC_tables_boot, -model), by = "model_no_orig")
  
  adult_AIC_tables_boot_summary <- 
    boot_out_list$adult_AIC_tables_boot %>%
    dplyr::group_by(model, model_no_orig) %>%
    dplyr::summarise(avg_Delta = mean(DeltaAICc),
                     IQR_Delta = IQR(DeltaAICc),
                     avg_Weight = mean(weight),
                     IQR_Weight = IQR(weight), .groups = 'drop') %>% 
    dplyr::arrange(avg_Delta) %>% 
    mutate(#model_number = as.numeric(model),
           stage = "adult") %>% 
    dplyr::left_join(., select(boot_out_list$adult_AIC_tables_boot, -model), by = "model_no_orig")
  
  # define the levels of each parameter according to their Delta AIC rank
  boot_out_list$nestling_AIC_tables_boot$S <- 
    factor(boot_out_list$nestling_AIC_tables_boot$S, 
           levels = unique(nestling_AIC_tables_boot_summary$S))
  
  boot_out_list$fledgling_AIC_tables_boot$S <- 
    factor(boot_out_list$fledgling_AIC_tables_boot$S, 
           levels = unique(fledgling_AIC_tables_boot_summary$S))
  
  boot_out_list$fledgling_AIC_tables_boot$p <- 
    factor(boot_out_list$fledgling_AIC_tables_boot$p, 
           levels = unique(fledgling_AIC_tables_boot_summary$p))
  
  boot_out_list$fledgling_AIC_tables_boot$r <- 
    factor(boot_out_list$fledgling_AIC_tables_boot$r, 
           levels = unique(fledgling_AIC_tables_boot_summary$r))
  
  boot_out_list$fledgling_AIC_tables_boot$F <- 
    factor(boot_out_list$fledgling_AIC_tables_boot$F, 
           levels = unique(fledgling_AIC_tables_boot_summary$F))
  
  boot_out_list$adult_AIC_tables_boot$p <- 
    factor(boot_out_list$adult_AIC_tables_boot$p, 
           levels = unique(adult_AIC_tables_boot_summary$p))
  
  boot_out_list$nestling_AIC_tables_boot <- 
    boot_out_list$nestling_AIC_tables_boot %>% 
    mutate(stage = "nestling")
  
  boot_out_list$fledgling_AIC_tables_boot <- 
    boot_out_list$fledgling_AIC_tables_boot %>% 
    mutate(stage = "fledgling")
  
  boot_out_list$adult_AIC_tables_boot <- 
    boot_out_list$adult_AIC_tables_boot %>% 
    mutate(stage = "adult")
  
  boot_out_list
  
}
