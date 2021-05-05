source("scripts/01_libraries.R")

# load functions
function.sources = list.files(path = "scripts", 
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
sapply(function.sources, source, .GlobalEnv)

# load capture histories
data.sources = list.files(path = "cooked_data", 
                          pattern="*ch.rds$", full.names = TRUE, 
                          ignore.case = TRUE)
sapply(data.sources, load, .GlobalEnv)
nestling_ch_input = Black_Coucal_nestling_Known_ch

fledgling_ch_input = Black_Coucal_fledgling_Burnham_ch
num_boot = 2
bootstrap_name = "BC_boot_test"
species = "BC"
prefix_number = "BC_boot_test_"
iter_add = 1

bootstrap_burnham_model <- 
  function(fledgling_ch_input, num_boot,
           bootstrap_name,
           species,
           prefix_number, iter_add) {
    
    # specify the bootstrapped data samples (from the previous function)
    fledgling_ch <- 
      fledgling_ch_input %>% 
      dplyr::group_by(nest_ID) %>%
      dplyr::sample_n(1) %>% 
      dplyr::mutate(sex = as.factor(sex)) %>%
      dplyr::ungroup() %>% 
      dplyr::select(ch, sex)
    
    
    # process the fledgling data as a "Burnham" joint live-dead encounter data
    coucal_fledgling.proc <- 
      RMark::process.data(data = fledgling_ch,
                          model = "Burnham",
                          groups = "sex")
    
    # make design matrix
    coucal_fledgling.ddl <- 
      RMark::make.design.data(coucal_fledgling.proc)
    
    # run the MARK model
    rmark_model <- 
      RMark::mark(model = "Burnham", groups = "sex", 
                  data = coucal_fledgling.proc, 
                  ddl = coucal_fledgling.ddl,
                  model.parameters = list(S = list(formula = ~sex * (I(Time) + I(Time^2))),
                                          p = list(formula = ~Time),
                                          F = list(formula = ~Time),
                                          r = list(formula = ~1)),
                  delete = TRUE, 
                  wrap = FALSE, threads = 1, brief = TRUE,
                  silent = TRUE, output = FALSE, prefix = prefix_number)
    
    extract_rmark_model_output <- 
      function(rmark_model_output, stage_name){
        # extract and wrangle reals from model output 
        if(stage_name == "nestling"){
          reals <- 
            rmark_model_output$results$real %>% 
            bind_cols(data.frame(str_split_fixed(rownames(.), " ", 
                                                 n = 5)), .) %>% 
            mutate(age = as.integer(unlist(str_extract_all(X3,"[0-9]+"))),
                   sex = as.factor(ifelse(unlist(str_extract_all(X2,"[FM]")) == "F", 
                                          "Female", "Male"))) %>% 
            select(sex, age, estimate) %>% 
            mutate(iter = num_boot + ((iter_add - 1) * niter),
                   species = species)
        }
        else if(stage_name == "fledgling"){
          reals <- 
            rmark_model_output$results$real %>% 
            bind_cols(data.frame(str_split_fixed(rownames(.), " ", 
                                                 n = 5)), .) %>% 
            mutate(age = as.integer(unlist(str_extract_all(X4,"[0-9]+"))),
                   sex = as.factor(ifelse(unlist(str_extract_all(X2,"[FM]")) == "F", 
                                          "Female", "Male")),
                   parameter = as.factor(unlist(str_extract_all(X1,"[SpFr]")))) %>% 
            select(parameter, sex, age, estimate) %>% 
            mutate(iter = num_boot + ((iter_add - 1) * niter),
                   species = species)
        }
        else{
          reals <- 
            rmark_model_output[[mod_num]]$results$real %>% 
            bind_cols(data.frame(str_split_fixed(rownames(.), " ", n = 5)), .) %>% 
            filter(X1 == "Phi") %>% 
            mutate(sex = as.factor(ifelse(unlist(str_extract_all(X2,"[FM]")) == "F", "Female", "Male"))) %>% 
            select(sex, estimate) %>% 
            mutate(iter = num_boot + ((iter_add - 1) * niter),
                   species = species)
        }
        
        # extract and wrangle the beta estimates from the model output
        betas <- 
          rmark_model_output$results$beta %>% 
          mutate(statistic = rownames(.),
                 iter = num_boot + ((iter_add - 1) * niter),
                 species = species) %>% 
          mutate(parameter = ifelse(grepl(x = statistic, pattern = "S:"), "S",
                                    ifelse(grepl(x = statistic, pattern = "p:"), "p",
                                           ifelse(grepl(x = statistic, pattern = "r:"), "r",
                                                  ifelse(grepl(x = statistic, pattern = "F:"), "F", 
                                                         ifelse(grepl(x = statistic, pattern = "Phi:"), "Phi", "XXX"))))),
                 variable = str_sub(statistic, start = 3, end = nchar(statistic))) %>% 
          # mutate(variable = ifelse(grepl(x = statistic, pattern = "\\<Intercept\\>"), "Intercept",
          #                          ifelse(grepl(x = statistic, pattern = "\\<sexM:I(Time)\\>"), "sexM:linear-age",
          #                                 ifelse(grepl(x = statistic, pattern = "\\<sexM:I(Time^2)\\>"), "sexM:quadratic-age",
          #                                        ifelse(grepl(x = statistic, pattern = "\\<sexM:Time\\>"), "sexM:Time",
          #                                               ifelse(grepl(x = statistic, pattern = "\\<sexM\\>"), "sexM",
          #                                                      ifelse(grepl(x = statistic, pattern = "\\<Time\\>"), "Time",
          #                                                             ifelse(grepl(x = statistic, pattern = "\\<I(Time)\\>"), "linear-age", 
          #                                                                    ifelse(grepl(x = statistic, pattern = "\\<I(Time^2)\\>"), "quad-age","XXX"))))))))) %>% 
          select(-statistic)
        
        # consolidate the AIC model selection results
        # AIC_table <- 
        #   rmark_model_output$results$ %>% 
        #   mutate(iter = num_boot + ((iter_add - 1) * niter),
        #          species = species,
        #          model_no_orig = as.numeric(rownames(.))) %>% 
        #   mutate(model_no_rank = as.numeric(rownames(.)))
        
        # consolidate the all output into a list
        survival_model_output_list <- 
          list(reals = reals,
               betas = betas,
               lnl = rmark_model_output$results$lnl,
               deviance = rmark_model_output$results$deviance,
               deviance.df = rmark_model_output$results$deviance.df,
               npar = rmark_model_output$results$npar,
               npar.unadjusted = rmark_model_output$results$npar.unadjusted,
               n = rmark_model_output$results$n,
               AICc = rmark_model_output$results$AICc,
               AICc.unadjusted = rmark_model_output$results$AICc.unadjusted,
               singular = rmark_model_output$results$singular)
        
        survival_model_output_list
      }
    
    # extract and format survival rates 
    fledgling_rmark_out <- 
      extract_rmark_model_output(rmark_model_output = rmark_model, 
                                 stage_name = "fledgling")
    
    bootstrap_results_list <- 
      list(fledgling_survival_model_output = fledgling_rmark_out,
           bootstrapped_data = fledgling_ch)
    
    # extract file directory
    file_directory <- paste0(getwd(), "/output/bootstraps/single_models/raw/", bootstrap_name, "_", num_boot, ".Rds")
    
    # save the results of the iteration
    save(bootstrap_results_list, file = file_directory)
    
    bootstrap_results_list
  }

bootstrap_known_model <- 
  function(ch_input, num_boot,
           bootstrap_name,
           species,
           prefix_number, iter_add) {
    
    # specify the bootstrapped data samples (from the previous function)
    boot_ch <- 
      ch_input %>% 
      dplyr::group_by(nest_ID) %>%
      dplyr::sample_n(1) %>% 
      dplyr::mutate(sex = as.factor(sex)) %>%
      dplyr::ungroup() %>% 
      dplyr::select(ch, sex)
    
    
    # process the fledgling data as a "Burnham" joint live-dead encounter data
    boot_known.proc <- 
      RMark::process.data(data = boot_ch,
                          model = "Known",
                          groups = "sex")
    
    # make design matrix
    boot_known.ddl <- 
      RMark::make.design.data(boot_known.proc)
    
    # run the MARK model
    rmark_model <- 
      RMark::mark(model = "Known", groups = "sex", 
                  data = boot_known.proc, 
                  ddl = boot_known.ddl,
                  model.parameters = list(S = list(formula = ~sex * (I(Time) + I(Time^2)))),
                  delete = TRUE, 
                  wrap = FALSE, threads = 1, brief = TRUE,
                  silent = TRUE, output = FALSE, prefix = prefix_number)
    
    extract_rmark_model_output <- 
      function(rmark_model_output, stage_name){
        # extract and wrangle reals from model output 
        if(stage_name == "nestling"){
          reals <- 
            rmark_model_output$results$real %>% 
            bind_cols(data.frame(str_split_fixed(rownames(.), " ", 
                                                 n = 5)), .) %>% 
            mutate(age = as.integer(unlist(str_extract_all(X3,"[0-9]+"))),
                   sex = as.factor(ifelse(unlist(str_extract_all(X2,"[FM]")) == "F", 
                                          "Female", "Male"))) %>% 
            select(sex, age, estimate) %>% 
            mutate(iter = num_boot + ((iter_add - 1) * niter),
                   species = species)
        }
        else if(stage_name == "fledgling"){
          reals <- 
            rmark_model_output$results$real %>% 
            bind_cols(data.frame(str_split_fixed(rownames(.), " ", 
                                                 n = 5)), .) %>% 
            mutate(age = as.integer(unlist(str_extract_all(X4,"[0-9]+"))),
                   sex = as.factor(ifelse(unlist(str_extract_all(X2,"[FM]")) == "F", 
                                          "Female", "Male")),
                   parameter = as.factor(unlist(str_extract_all(X1,"[SpFr]")))) %>% 
            select(parameter, sex, age, estimate) %>% 
            mutate(iter = num_boot + ((iter_add - 1) * niter),
                   species = species)
        }
        else{
          reals <- 
            rmark_model_output[[mod_num]]$results$real %>% 
            bind_cols(data.frame(str_split_fixed(rownames(.), " ", n = 5)), .) %>% 
            filter(X1 == "Phi") %>% 
            mutate(sex = as.factor(ifelse(unlist(str_extract_all(X2,"[FM]")) == "F", "Female", "Male"))) %>% 
            select(sex, estimate) %>% 
            mutate(iter = num_boot + ((iter_add - 1) * niter),
                   species = species)
        }
        
        # extract and wrangle the beta estimates from the model output
        betas <- 
          rmark_model_output$results$beta %>% 
          mutate(statistic = rownames(.),
                 iter = num_boot + ((iter_add - 1) * niter),
                 species = species) %>% 
          mutate(parameter = ifelse(grepl(x = statistic, pattern = "S:"), "S",
                                    ifelse(grepl(x = statistic, pattern = "p:"), "p",
                                           ifelse(grepl(x = statistic, pattern = "r:"), "r",
                                                  ifelse(grepl(x = statistic, pattern = "F:"), "F", 
                                                         ifelse(grepl(x = statistic, pattern = "Phi:"), "Phi", "XXX"))))),
                 variable = str_sub(statistic, start = 3, end = nchar(statistic))) %>% 
          # mutate(variable = ifelse(grepl(x = statistic, pattern = "\\<Intercept\\>"), "Intercept",
          #                          ifelse(grepl(x = statistic, pattern = "\\<sexM:I(Time)\\>"), "sexM:linear-age",
          #                                 ifelse(grepl(x = statistic, pattern = "\\<sexM:I(Time^2)\\>"), "sexM:quadratic-age",
          #                                        ifelse(grepl(x = statistic, pattern = "\\<sexM:Time\\>"), "sexM:Time",
          #                                               ifelse(grepl(x = statistic, pattern = "\\<sexM\\>"), "sexM",
          #                                                      ifelse(grepl(x = statistic, pattern = "\\<Time\\>"), "Time",
          #                                                             ifelse(grepl(x = statistic, pattern = "\\<I(Time)\\>"), "linear-age", 
          #                                                                    ifelse(grepl(x = statistic, pattern = "\\<I(Time^2)\\>"), "quad-age","XXX"))))))))) %>% 
          select(-statistic)
        
        # consolidate the AIC model selection results
        # AIC_table <- 
        #   rmark_model_output$results$ %>% 
        #   mutate(iter = num_boot + ((iter_add - 1) * niter),
        #          species = species,
        #          model_no_orig = as.numeric(rownames(.))) %>% 
        #   mutate(model_no_rank = as.numeric(rownames(.)))
        
        # consolidate the all output into a list
        survival_model_output_list <- 
          list(reals = reals,
               betas = betas,
               lnl = rmark_model_output$results$lnl,
               deviance = rmark_model_output$results$deviance,
               deviance.df = rmark_model_output$results$deviance.df,
               npar = rmark_model_output$results$npar,
               npar.unadjusted = rmark_model_output$results$npar.unadjusted,
               n = rmark_model_output$results$n,
               AICc = rmark_model_output$results$AICc,
               AICc.unadjusted = rmark_model_output$results$AICc.unadjusted,
               singular = rmark_model_output$results$singular)
        
        survival_model_output_list
      }
    
    # extract and format survival rates 
    known_rmark_out <- 
      extract_rmark_model_output(rmark_model_output = rmark_model, 
                                 stage_name = "nestling")
    
    bootstrap_results_list <- 
      list(known_survival_model_output = known_rmark_out,
           bootstrapped_data = boot_ch)
    
    # extract file directory
    file_directory <- paste0(getwd(), "/output/bootstraps/single_models/raw/", bootstrap_name, "_", num_boot, ".Rds")
    
    # save the results of the iteration
    save(bootstrap_results_list, file = file_directory)
    
    bootstrap_results_list
  }

run_bootstrap_burnham_model <- 
  function(num_boot, fledgling_ch_input,
           bootstrap_name,
           species, iter_add, prefix_number){
  # run the survival analysis and ASR deduction on the sampled data
  result <- bootstrap_burnham_model(fledgling_ch_input = fledgling_ch_input, 
                                   num_boot = num_boot, 
                                   bootstrap_name = bootstrap_name,
                                   species = species,
                                   iter_add = iter_add,
                                   prefix_number = prefix_number)
  result
  }

run_bootstrap_known_model <- 
  function(num_boot, ch_input,
           bootstrap_name,
           species, iter_add, prefix_number){
    # run the survival analysis and ASR deduction on the sampled data
    result <- bootstrap_known_model(ch_input = fledgling_ch_input, 
                                      num_boot = num_boot, 
                                      bootstrap_name = bootstrap_name,
                                      species = species,
                                      iter_add = iter_add,
                                      prefix_number = prefix_number)
    result
  }

# boot_out_wrangle() loads all the four output rds files from the parallel 
# bootstraping proceedure and consolidates them into a single object
single_mod_boot_out_wrangle <- 
  function(species, niter, rds_file = "_bootstrap_burnham_model_bootstrap.rds"){
  
  # specify directory string
  boot_path <- paste0("output/bootstraps/single_models/combined/", 
                            species, rds_file)
  
  # load each of the four bootstrap output files
  load(boot_path)

  # determine which species to wrangle
  if(species == "BC"){
    survival_ASR_bootstrap <- BC_bootstrap_burnham_model_bootstrap
  }
  else{
    survival_ASR_bootstrap <- WBC_bootstrap_burnham_model_bootstrap
  }
  
  # bind all the fledgling daily survival rates output together
  fledgling_reals_survival_rates_boot <-
   do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                         function(x) survival_ASR_bootstrap["fledgling_survival_model_output", x]$fledgling_survival_model_output$reals)) %>%
    arrange(iter, parameter, sex, age)
  
  # bind all the fledgling beta estimates output together
  fledgling_betas_survival_rates_boot <-
    do.call(rbind, lapply(seq(from = 1, to = niter, by = 1),
                          function(x) survival_ASR_bootstrap["fledgling_survival_model_output", x]$fledgling_survival_model_output$betas))
  
  # tidy all the bound data together into a mega list of results
  boot_out_tidy <- 
    list(fledgling_reals_survival_rates_boot = fledgling_reals_survival_rates_boot,
         fledgling_betas_survival_rates_boot = fledgling_betas_survival_rates_boot)
  
  boot_out_tidy
  }

niter = 1000
# set.seed(14)

rnorm(n = 2, mean = 0, sd = 1)

nchar(Black_Coucal_fledgling_Burnham_ch$ch)/2

# run bootstrap procedure on Black Coucals
BC_bootstrap_burnham_model_bootstrap <-
  pbsapply(1:niter, run_bootstrap_burnham_model, 
           fledgling_ch_input = Black_Coucal_fledgling_Burnham_ch,
           bootstrap_name = "BC_boot_test",
           species = "BC",
           iter_add = 1,
           prefix_number = "BC_boot_test_")

save(BC_bootstrap_burnham_model_bootstrap,
     file = "output/bootstraps/single_models/combined/BC_bootstrap_burnham_model_bootstrap.rds")

# run bootstrap procedure on Black Coucals
BC_bootstrap_known_model_bootstrap_F <-
  pbsapply(1:niter, run_bootstrap_known_model, 
           ch_input = Black_Coucal_fledgling_Burnham_ch,
           bootstrap_name = "BC_boot_test_known_",
           species = "BC",
           iter_add = 1,
           prefix_number = "BC_boot_test_known_")

save(BC_bootstrap_known_model_bootstrap,
     file = "output/bootstraps/single_models/combined/BC_bootstrap_known_model_bootstrap.rds")

# run bootstrap procedure on Black Coucals
BC_bootstrap_known_model_bootstrap_N <-
  pbsapply(1:niter, run_bootstrap_known_model, 
           ch_input = Black_Coucal_nestling_Known_ch,
           bootstrap_name = "BC_boot_test_known_N_",
           species = "BC",
           iter_add = 1,
           prefix_number = "BC_boot_test_known_N_")

save(BC_bootstrap_known_model_bootstrap_N,
     file = "output/bootstraps/single_models/combined/BC_bootstrap_known_model_bootstrap_N.rds")

# clean up the output from the bootstrap procedure and save as rds
BC_single_mod_boot_out <- 
  single_mod_boot_out_wrangle(species = "BC", niter = 1000)

# clean up the output from the bootstrap procedure and save as rds
BC_fledgling_single_mod_boot_out <- 
  single_mod_boot_out_wrangle(species = "BC", niter = 1000, rds_file = "_bootstrap_known_model_bootstrap.rds")

# clean up the output from the bootstrap procedure and save as rds
BC_nestling_single_mod_boot_out <- 
  single_mod_boot_out_wrangle(species = "BC", niter = 1000, rds_file = "_bootstrap_known_model_bootstrap_N.rds")

test <- 
  BC_single_mod_boot_out$fledgling_reals_survival_rates_boot %>% 
  filter(parameter == "S")

test <- 
  BC_nestling_single_mod_boot_out$fledgling_reals_survival_rates_boot %>% 
  filter(parameter == "S")

brewer.pal(6, "Dark2")
ggplot() +
  luke_theme +
  theme(legend.position = "none") +
  geom_vline(xintercept = 36.85 - 15, linetype = "dashed", alpha = 0.5, color = "orange") +
  annotate("rect", xmin = 34.71 - 15, xmax = 38.99 - 15, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "orange") +
  annotate("rect", xmin = 32.73 - 15, xmax = 36.94 - 15, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "green") +
  
  geom_vline(xintercept = 34.85 - 15, linetype = "dashed", alpha = 0.5, color = "green") +
  geom_line(data = test, 
            aes(x = age, y = estimate, 
                group = interaction(iter, sex), 
                color = sex),
            alpha = 0.05) +
  scale_colour_brewer(palette = "Dark2", direction = -1) + 
  ylab("Estimated daily survival rate") +
  xlab("Age (Days since hatching)") + 
  scale_y_continuous(limits = c(0.9, 1),
                     breaks = seq(from = 0.9, to = 1, by = 0.025)) +
  scale_x_continuous(limits = c(0, 55), 
                     breaks = seq(0, 55, by = 5), 
                     expand = c(0.01, 0.01),
                     labels = as.character(seq(15, 70, by = 5)))


  #+
  geom_line(data = mod_dist_fits, aes(x = age, y = fit, color = sex),
            lwd = 0.5) +
  geom_ribbon(data = mod_dist_fits, aes(x = age, ymax = upper, ymin = lower, fill = sex),
              lwd = 1, alpha = 0.25) +
  facet_grid(species ~ ., 
             labeller = labeller(species = species.labs)) +
  ylab(expression(paste("Distance moved between observations (log10)" %+-%  "95% CI", sep = ""))) +
  xlab("Age since leaving nest") +
  scale_color_manual(values = plot_palette_sex,
                     name = "Sex",
                     breaks = c("Female", "Male"),
                     labels = c("Female", "Male")) +
  scale_fill_manual(values = plot_palette_sex,
                    name = "Sex",
                    breaks = c("Female", "Male"),
                    labels = c("Female", "Male"))
  