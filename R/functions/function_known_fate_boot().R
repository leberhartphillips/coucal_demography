known_fate_boot <- 
  function(ch_input, num_boot,
           bootstrap_name,
           species, stage_name,
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
        # if(stage_name == "nestling"){
          reals <- 
            rmark_model_output$results$real %>% 
            bind_cols(data.frame(str_split_fixed(rownames(.), " ", 
                                                 n = 5)), .) %>% 
            mutate(age = as.integer(unlist(str_extract_all(X3,"[0-9]+"))),
                   sex = as.factor(ifelse(unlist(str_extract_all(X2,"[FM]")) == "F", 
                                          "Female", "Male")),
                   parameter = as.factor(unlist(str_extract_all(X1,"[SpFr]"))),
                   stage = stage_name) %>% 
            select(parameter, sex, stage, age, estimate) %>% 
            mutate(iter = num_boot + ((iter_add - 1) * niter),
                   species = species)
        # }
        # else if(stage_name == "fledgling"){
        #   reals <- 
        #     rmark_model_output$results$real %>% 
        #     bind_cols(data.frame(str_split_fixed(rownames(.), " ", 
        #                                          n = 5)), .) %>% 
        #     mutate(age = as.integer(unlist(str_extract_all(X4,"[0-9]+"))),
        #            sex = as.factor(ifelse(unlist(str_extract_all(X2,"[FM]")) == "F", 
        #                                   "Female", "Male")),
        #            parameter = as.factor(unlist(str_extract_all(X1,"[SpFr]")))) %>% 
        #     select(parameter, sex, age, estimate) %>% 
        #     mutate(iter = num_boot + ((iter_add - 1) * niter),
        #            species = species)
        # }
        # else{
        #   reals <- 
        #     rmark_model_output[[mod_num]]$results$real %>% 
        #     bind_cols(data.frame(str_split_fixed(rownames(.), " ", n = 5)), .) %>% 
        #     filter(X1 == "Phi") %>% 
        #     mutate(sex = as.factor(ifelse(unlist(str_extract_all(X2,"[FM]")) == "F", "Female", "Male"))) %>% 
        #     select(sex, estimate) %>% 
        #     mutate(iter = num_boot + ((iter_add - 1) * niter),
        #            species = species)
        # }
        
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
                                 stage_name = stage_name)
    
    bootstrap_results_list <- 
      list(known_fate_model_output = known_rmark_out,
           bootstrapped_data = boot_ch)
    
    # extract file directory
    file_directory <- paste0(getwd(), "/output/bootstraps/single_models/raw/", bootstrap_name, "_", num_boot, ".Rds")
    
    # save the results of the iteration
    save(bootstrap_results_list, file = file_directory)
    
    bootstrap_results_list
  }