source("R/project/project_libraries.R")

# boot_out_wrangle() loads all the four output rds files from the parallel 
# bootstraping proceedure and consolidates them into a single object
boot_out_wrangle <- function(species, niter){
  
  # specify directory string
  boot_1_out_path <- paste0("output/bootstraps/", 
                            species,"_survival_ASR_bootstrap_result_one.rds")
  boot_2_out_path <- paste0("output/bootstraps/", 
                            species,"_survival_ASR_bootstrap_result_two.rds")
  boot_3_out_path <- paste0("output/bootstraps/", 
                            species,"_survival_ASR_bootstrap_result_three.rds")
  boot_4_out_path <- paste0("output/bootstraps/", 
                            species,"_survival_ASR_bootstrap_result_four.rds")
  
  # load each of the four bootstrap output files
  load(boot_1_out_path)
  load(boot_2_out_path)
  load(boot_3_out_path)
  load(boot_4_out_path)
  
  # determine which species to wrangle
  if(species == "BC"){
    survival_ASR_bootstrap_result_one <- BC_survival_ASR_bootstrap_result_one
    survival_ASR_bootstrap_result_two <- BC_survival_ASR_bootstrap_result_two
    survival_ASR_bootstrap_result_three <- BC_survival_ASR_bootstrap_result_three
    survival_ASR_bootstrap_result_four <- BC_survival_ASR_bootstrap_result_four
  }
  else{
    survival_ASR_bootstrap_result_one <- WBC_survival_ASR_bootstrap_result_one
    survival_ASR_bootstrap_result_two <- WBC_survival_ASR_bootstrap_result_two
    survival_ASR_bootstrap_result_three <- WBC_survival_ASR_bootstrap_result_three
    survival_ASR_bootstrap_result_four <- WBC_survival_ASR_bootstrap_result_four
  }
  
  # bind all the summary vital rate output together
  survival_rates_boot <-
    rbind(do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_one["coucal_vital_rates", x]$coucal_vital_rates)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_two["coucal_vital_rates", x]$coucal_vital_rates)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_three["coucal_vital_rates", x]$coucal_vital_rates)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_four["coucal_vital_rates", x]$coucal_vital_rates)))
  
  # bind all the ASR output together
  ASR_boot <-
    rbind(do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_one["ASR_SSD", x]$ASR_SSD$ASR)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_two["ASR_SSD", x]$ASR_SSD$ASR)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_three["ASR_SSD", x]$ASR_SSD$ASR)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_four["ASR_SSD", x]$ASR_SSD$ASR))) %>%
    as.data.frame() %>% 
    dplyr::rename(ASR_estimate = M_Adult) %>% 
    mutate(species = species,
           iter = row.names(.))
  
  # bind all the nestling daily survival rates output together
  nestling_reals_survival_rates_boot <-
    rbind(do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_one["nestling_survival_model_output", x]$nestling_survival_model_output$reals)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_two["nestling_survival_model_output", x]$nestling_survival_model_output$reals)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_three["nestling_survival_model_output", x]$nestling_survival_model_output$reals)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_four["nestling_survival_model_output", x]$nestling_survival_model_output$reals))) %>%
    arrange(iter, sex, age)
  
  # bind all the fledgling daily survival rates output together
  fledgling_reals_survival_rates_boot <-
    rbind(do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_one["fledgling_survival_model_output", x]$fledgling_survival_model_output$reals)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_two["fledgling_survival_model_output", x]$fledgling_survival_model_output$reals)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_three["fledgling_survival_model_output", x]$fledgling_survival_model_output$reals)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_four["fledgling_survival_model_output", x]$fledgling_survival_model_output$reals))) %>%
    arrange(iter, sex, age)
  
  # bind all the adult rates output together
  adult_reals_survival_rates_boot <-
    rbind(do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_one["adult_survival_model_output", x]$adult_survival_model_output$reals)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_two["adult_survival_model_output", x]$adult_survival_model_output$reals)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_three["adult_survival_model_output", x]$adult_survival_model_output$reals)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                                function(x) survival_ASR_bootstrap_result_four["adult_survival_model_output", x]$adult_survival_model_output$reals))) %>%
    arrange(iter, sex)
  
  # bind all the nestling beta estimates output together
  nestling_betas_survival_rates_boot <-
    rbind(do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_one["nestling_survival_model_output", x]$nestling_survival_model_output$betas)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_two["nestling_survival_model_output", x]$nestling_survival_model_output$betas)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_three["nestling_survival_model_output", x]$nestling_survival_model_output$betas)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_four["nestling_survival_model_output", x]$nestling_survival_model_output$betas)))
  
  # bind all the fledgling beta estimates output together
  fledgling_betas_survival_rates_boot <-
    rbind(do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_one["fledgling_survival_model_output", x]$fledgling_survival_model_output$betas)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_two["fledgling_survival_model_output", x]$fledgling_survival_model_output$betas)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_three["fledgling_survival_model_output", x]$fledgling_survival_model_output$betas)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_four["fledgling_survival_model_output", x]$fledgling_survival_model_output$betas)))
  
  # bind all the adult beta estimates output together
  adult_betas_survival_rates_boot <-
    rbind(do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_one["adult_survival_model_output", x]$adult_survival_model_output$betas)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_two["adult_survival_model_output", x]$adult_survival_model_output$betas)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_three["adult_survival_model_output", x]$adult_survival_model_output$betas)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_four["adult_survival_model_output", x]$adult_survival_model_output$betas)))
  
  # bind all the nestling AIC table output together
  nestling_AIC_tables_boot <-
    rbind(do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_one["nestling_survival_model_output", x]$nestling_survival_model_output$AIC_table)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_two["nestling_survival_model_output", x]$nestling_survival_model_output$AIC_table)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_three["nestling_survival_model_output", x]$nestling_survival_model_output$AIC_table)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_four["nestling_survival_model_output", x]$nestling_survival_model_output$AIC_table)))
  
  # bind all the fledgling AIC table output together
  fledgling_AIC_tables_boot <-
    rbind(do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_one["fledgling_survival_model_output", x]$fledgling_survival_model_output$AIC_table)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_two["fledgling_survival_model_output", x]$fledgling_survival_model_output$AIC_table)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_three["fledgling_survival_model_output", x]$fledgling_survival_model_output$AIC_table)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_four["fledgling_survival_model_output", x]$fledgling_survival_model_output$AIC_table)))
  
  # bind all the adult AIC table output together
  adult_AIC_tables_boot <-
    rbind(do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_one["adult_survival_model_output", x]$adult_survival_model_output$AIC_table)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_two["adult_survival_model_output", x]$adult_survival_model_output$AIC_table)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_three["adult_survival_model_output", x]$adult_survival_model_output$AIC_table)),
          do.call(rbind, lapply(seq(from = 1, to = niter, by = 1), 
                                function(x) survival_ASR_bootstrap_result_four["adult_survival_model_output", x]$adult_survival_model_output$AIC_table)))
  
  # tidy all the bound data together into a mega list of results
  boot_out_tidy <- 
    list(survival_rates_boot = survival_rates_boot,
         ASR_boot = ASR_boot,
         nestling_reals_survival_rates_boot = nestling_reals_survival_rates_boot,
         fledgling_reals_survival_rates_boot = fledgling_reals_survival_rates_boot,
         adult_reals_survival_rates_boot = adult_reals_survival_rates_boot,
         nestling_betas_survival_rates_boot = nestling_betas_survival_rates_boot,
         fledgling_betas_survival_rates_boot = fledgling_betas_survival_rates_boot,
         adult_betas_survival_rates_boot = adult_betas_survival_rates_boot,
         nestling_AIC_tables_boot = nestling_AIC_tables_boot,
         fledgling_AIC_tables_boot = fledgling_AIC_tables_boot,
         adult_AIC_tables_boot = adult_AIC_tables_boot)
  
  boot_out_tidy
}