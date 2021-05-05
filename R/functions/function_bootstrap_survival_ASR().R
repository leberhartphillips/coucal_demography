source("scripts/01_libraries.R")
source("scripts/02_func_matrix_ASR().R")
source("scripts/02_func_coucal_matrix().R")

# bootstrap_survival_ASR() runs the survival analyses and estimates the ASR of 
# the bootstrapped sample created from bootstrap_data(). 
# coucal_boot_list :  the output list from bootstrap_data()
# ISR : the immigration sex ratio (i.e., the observed adult sex ratio during the
# breeding season)
# immigrant_pop_size : the number of immigrants entering the population at the 
# start of each season
# default is a vector with 10 individuals in each stage)
# h : harem size index
# k : clutch size
# HSR : hatching sex ratio
# num_boot : bootstrap number in the loop (do not specify)
# species : coucal species name
# iter_add : addition to bootstrap number in the loop (do not specify)
# bootstrap_name : ID of bootstrap iteration (do not specify)
# flight_age : age at first flight in days
# first_year : the first year of the study period (e.g., 2001)
# egg_survival : hatching probability

bootstrap_survival_ASR <- 
  function(coucal_boot_list, num_boot, iter_add,
           egg_survival,
           ISR,
           immigrant_pop_size = 100,
           HSR = 0.5,
           h,
           k,
           flight_age,
           first_year,
           bootstrap_name,
           species,
           prefix_number) {
    
    # specify the bootstrapped data samples (from the previous function)
    nestling_ch <- coucal_boot_list[["nestling_boot"]]
    fledgling_ch <- coucal_boot_list[["fledgling_boot"]]
    adult_ch <- coucal_boot_list[["adult_boot"]]
    
    # clean up capture histories
    fledgling_ch <- 
      fledgling_ch %>% 
      ungroup() %>% 
      as.data.frame()
    
    nestling_ch <-    
      nestling_ch %>% 
      ungroup() %>% 
      as.data.frame()
    
    adult_ch <-    
      adult_ch %>% 
      ungroup() %>% 
      as.data.frame() %>% 
      mutate(sex = as.factor(sex))
    
    # process the nestling data as a "Known" fate analysis
    coucal_nestling.proc <- process.data(data = nestling_ch,
                                         model = "Known",
                                         groups = c("sex"))
    
    # make design matrix
    coucal_nestling.ddl <- 
      make.design.data(coucal_nestling.proc)
    
    # create polynomial variations of time
    time <- c(1:(coucal_nestling.proc$nocc[1]+1))
    Quadratic <- time^2
    Cubic <- time^3
    quad_time <- data.frame(time, Quadratic, Cubic)
    
    # merge polynomials to the design matrix
    coucal_nestling.ddl$S <- 
      merge_design.covariates(coucal_nestling.ddl$S, quad_time, 
                              bygroup = FALSE, bytime = TRUE)
    
    # process the fledgling data as a "Burnham"  analysis
    coucal_fledgling.proc <- process.data(data = fledgling_ch,
                                          model = "Burnham",
                                          groups = c("sex"))
    
    # make design matrix
    coucal_fledgling.ddl <- 
      make.design.data(coucal_fledgling.proc)
    
    # create polynomial variations of time
    time <- c(1:(coucal_fledgling.proc$nocc[1]+1))
    Quadratic <- time^2
    Cubic <- time^3
    quad_time <- data.frame(time, Quadratic, Cubic)
    
    # merge polynomials to the design matrix
    coucal_fledgling.ddl$S <- 
      merge_design.covariates(coucal_fledgling.ddl$S, quad_time, 
                              bygroup = FALSE, bytime = TRUE)
    coucal_fledgling.ddl$F <- 
      merge_design.covariates(coucal_fledgling.ddl$F, quad_time, 
                              bygroup = FALSE, bytime = TRUE)
    coucal_fledgling.ddl$r <- 
      merge_design.covariates(coucal_fledgling.ddl$r, quad_time, 
                              bygroup = FALSE, bytime = TRUE)
    
    # process the adult data as a "CJS"  analysis
    coucal_adult.proc <- RMark::process.data(adult_ch, 
                                             model = "CJS",
                                             groups = c("sex"),
                                             begin.time = first_year)
    
    # Create the design matrix from the processed mark-recapture datasets
    coucal_adult.ddl <- RMark::make.design.data(coucal_adult.proc)
    
    nestling_survival_analysis_run <- 
      function(proc_data, design_data){
        # apriori model components for
        # S (survival probability):
        
        # sex- and age-specific model
        S.sexXTime = list(formula = ~sex * Time)
        
        # sex- and age-specific model (quadratic)
        S.sexXQuad = list(formula = ~sex * Quadratic)
        
        # sex- and age-specific model (cubic)
        S.sexXCube = list(formula = ~sex * Cubic)
        
        # create model list for all a priori models above
        cml <- create.model.list("Known")
        
        # run the models in program MARK
        model.list <-  RMark::mark.wrapper(cml, data = proc_data, 
                                           ddl = design_data, delete = TRUE, 
                                           wrap = FALSE, threads = 1, brief = TRUE,
                                           silent = TRUE, output = FALSE, prefix = prefix_number)
        
        # output the model list and store the results
        return(model.list)
      }
    
    fledgling_survival_analysis_run <- 
      function(proc_data, design_data){
        # apriori model components for
        # S (survival probability):
        
        # sex- and age-specific model
        S.sexXTime = list(formula = ~sex * Time)
        
        # sex- and age-specific model (quadratic)
        S.sexXQuad = list(formula = ~sex * Quadratic)
        
        # sex- and age-specific model (cubic)
        S.sexXCube = list(formula = ~sex * Cubic)
        
        # p (encounter probability):
        # null model
        p.dot = list(formula = ~1)
        
        # sex-specific model
        p.sex = list(formula = ~sex)
        
        # age-specific model
        p.Time = list(formula = ~Time)
        
        # F (site fidelity probability):
        # null model
        F.dot = list(formula = ~1)
        
        # sex-specific model
        F.sex = list(formula = ~sex)
        
        # age-specific model
        F.Time = list(formula = ~Time)
        
        # r (recovery probability)
        # null model
        r.dot = list(formula = ~1)
        
        # sex-specific model
        r.sex = list(formula = ~sex)
        
        # age-specific model
        r.time = list(formula = ~Time)
        
        # create model list for all a priori models above
        cml <- create.model.list("Burnham")
        
        # run the models in program MARK
        model.list <-  RMark::mark.wrapper(cml, data = proc_data, 
                                           ddl = design_data, delete = TRUE, 
                                           wrap = FALSE, threads = 1, brief = TRUE,
                                           silent = TRUE, output = FALSE, prefix = prefix_number)
        
        # output the model list and sotre the results
        return(model.list)
      }
    
    adult_survival_analysis_run = 
      function(proc_data, design_data){
        # sex- and stage-specific survival:
        Phi.sex = list(formula = ~ sex) 
        
        # Models exploring variation in encounter probability
        # constant:
        p.dot = list(formula =  ~ 1)
        
        # sex-dependent:
        p.sex = list(formula =  ~ sex)
        
        # factorial variation across year:
        p.Time = list(formula =  ~ time)
        
        # additive effects of sex and factorial year:
        p.sex_Time = list(formula =  ~ sex + time)
        
        # create a list of candidate models for all the a models above that begin with 
        # either "Phi." or "p."
        cml <-  RMark::create.model.list("CJS")
        
        # specify the data, design matrix, delete unneeded output files, and 
        # run the models in Program MARK
        model.list <-  RMark::mark.wrapper(cml, data = proc_data, 
                                           ddl = design_data, delete = TRUE, 
                                           wrap = FALSE, threads = 1, brief = TRUE,
                                           silent = TRUE, output = FALSE, prefix = prefix_number)
        
        # output the model list and sotre the results
        return(model.list)
      }
    
    # Run the models on the bootstrapped data
    nestling_survival_analysis_out <-
      nestling_survival_analysis_run(proc_data = coucal_nestling.proc,
                                     design_data = coucal_nestling.ddl)
    
    fledgling_survival_analysis_out <-
      fledgling_survival_analysis_run(proc_data = coucal_fledgling.proc,
                                      design_data = coucal_fledgling.ddl)
    
    adult_survival_analysis_out <-
      adult_survival_analysis_run(proc_data = coucal_adult.proc,
                                  design_data = coucal_adult.ddl)
    
    extract_top_model_output <- 
      function(rmark_output, stage_name, top_model = TRUE, mod_num){
        # Find the model number for the first ranked model of the AIC table
        if(top_model == TRUE){
          mod_num <- 
            as.numeric(rownames(rmark_output$model.table[1,]))
        }
        
        else{
          mod_num <- mod_num
        }
        
        # extract and wrangle reals from model output 
        if(stage_name == "nestling"){
          reals <- 
            rmark_output[[mod_num]]$results$real %>% 
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
            rmark_output[[mod_num]]$results$real %>% 
            bind_cols(data.frame(str_split_fixed(rownames(.), " ", 
                                                 n = 5)), .) %>% 
            mutate(age = as.integer(unlist(str_extract_all(X4,"[0-9]+"))),
                   sex = as.factor(ifelse(unlist(str_extract_all(X2,"[FM]")) == "F", 
                                          "Female", "Male"))) %>% 
            filter(X1 == "S") %>% 
            select(sex, age, estimate) %>% 
            mutate(iter = num_boot + ((iter_add - 1) * niter),
                   species = species)
        }
        else{
          reals <- 
            rmark_output[[mod_num]]$results$real %>% 
            bind_cols(data.frame(str_split_fixed(rownames(.), " ", n = 5)), .) %>% 
            filter(X1 == "Phi") %>% 
            mutate(sex = as.factor(ifelse(unlist(str_extract_all(X2,"[FM]")) == "F", "Female", "Male"))) %>% 
            select(sex, estimate) %>% 
            mutate(iter = num_boot + ((iter_add - 1) * niter),
                   species = species)
        }
        
        # extract and wrangle the beta estimates from the model output
        betas <- 
          rmark_output[[mod_num]]$results$beta %>% 
          mutate(statistic = rownames(.),
                 iter = num_boot + ((iter_add - 1) * niter),
                 species = species) %>% 
          mutate(parameter = ifelse(grepl(x = statistic, pattern = "S:"), "S",
                                    ifelse(grepl(x = statistic, pattern = "p:"), "p",
                                           ifelse(grepl(x = statistic, pattern = "r:"), "r",
                                                  ifelse(grepl(x = statistic, pattern = "F:"), "F", 
                                                         ifelse(grepl(x = statistic, pattern = "Phi:"), "Phi", "XXX"))))),
                 variable = ifelse(grepl(x = statistic, pattern = "Intercept"), "Intercept",
                                   ifelse(grepl(x = statistic, pattern = "sexM:Cubic"), "sexM:Cubic",
                                          ifelse(grepl(x = statistic, pattern = "sexM:Quadratic"), "sexM:Quadratic",
                                                 ifelse(grepl(x = statistic, pattern = "sexM:Time"), "sexM:Time",
                                                        ifelse(grepl(x = statistic, pattern = "sexM"), "sexM",
                                                               ifelse(grepl(x = statistic, pattern = "Time"), "Time",
                                                                      ifelse(grepl(x = statistic, pattern = "Cubic"), "Cubic", 
                                                                             ifelse(grepl(x = statistic, pattern = "Quadratic"), "Quadratic","XXX"))))))))) %>% 
          select(-statistic)
        
        # consolidate the AIC model selection results
        AIC_table <- 
          rmark_output$model.table %>% 
          mutate(iter = num_boot + ((iter_add - 1) * niter),
                 species = species,
                 model_no_orig = as.numeric(rownames(.))) %>% 
          mutate(model_no_rank = as.numeric(rownames(.)))
        
        # consolidate the all output into a list
        survival_model_output_list <- 
          list(reals = reals,
               betas = betas,
               AIC_table = AIC_table)
        
        survival_model_output_list
      }
    
    # extract and format survival rates from juvenile and adult model output
    nestling_survival_model_output_list <- 
      extract_top_model_output(rmark_output = nestling_survival_analysis_out, 
                               stage_name = "nestling")
    
    fledgling_survival_model_output_list <- 
      extract_top_model_output(rmark_output = fledgling_survival_analysis_out, 
                               stage_name = "fledgling")
    
    adult_survival_model_output_list <- 
      extract_top_model_output(rmark_output = adult_survival_analysis_out, 
                               stage_name = "adult")
    
    # transform the daily nestling survival (DCS) to apparent fledgling success
    # by calculating the product of all DCS estimates:
    coucal_nestling_survival <-
      nestling_survival_model_output_list$reals %>% 
      group_by(sex) %>% 
      dplyr::summarise(value = prod(estimate), .groups = 'drop') %>%
      mutate(stage = "nestling",
             rate = "survival")
    
    coucal_groundling_survival <-
      fledgling_survival_model_output_list$reals %>% 
      filter(age < flight_age) %>% 
      group_by(sex) %>% 
      dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
      mutate(stage = "groundling",
             rate = "survival")
    
    coucal_fledgling_survival <-
      fledgling_survival_model_output_list$reals %>% 
      filter(age >= flight_age) %>% 
      group_by(sex) %>% 
      dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
      mutate(stage = "fledgling",
             rate = "survival")
    
    coucal_adult_survival <-
      adult_survival_model_output_list$reals %>%
      select(-species, -iter) %>% 
      mutate(stage = "adult",
             rate = "survival") %>% 
      dplyr::rename(value = estimate)
    
    adult_F_immigrants <- immigrant_pop_size * (1 - ISR)
    adult_M_immigrants <- immigrant_pop_size * ISR
    
    coucal_adult_immigration <- 
      data.frame(sex = c("Female", "Male"),
                 value = c(adult_F_immigrants, adult_M_immigrants),
                 stage = c("adult"),
                 rate = c("immigration"))
    
    coucal_egg_survival <- 
      data.frame(sex = NA,
                 value = egg_survival,
                 stage = c("egg"),
                 rate = c("survival"))
    
    # Bind the juvenile and adult dataframe with the nestlings
    coucal_vital_rates <- 
      bind_rows(coucal_egg_survival,
                coucal_nestling_survival,
                coucal_groundling_survival,
                coucal_fledgling_survival,
                coucal_adult_survival,
                coucal_adult_immigration) %>% 
      as.data.frame() %>% 
      mutate(iter = num_boot + ((iter_add - 1) * niter),
             species = species)
    
    # Create a list of demographic rates from the survival analyses above
    demographic_rates <- list(Egg_survival = coucal_vital_rates[1, 2],
                              F_Nestling_survival = coucal_vital_rates[2, 2],
                              F_Groundling_survival = coucal_vital_rates[4, 2],
                              F_Fledgling_survival = coucal_vital_rates[6, 2],
                              F_Adult_survival = coucal_vital_rates[8, 2],
                              F_Adult_immigration = coucal_vital_rates[10, 2],
                              M_Nestling_survival = coucal_vital_rates[3, 2],
                              M_Groundling_survival = coucal_vital_rates[5, 2],
                              M_Fledgling_survival = coucal_vital_rates[7, 2],
                              M_Adult_survival = coucal_vital_rates[9, 2],
                              M_Adult_immigration = coucal_vital_rates[11, 2],
                              
                              # Define hatching sex ratio
                              HSR = HSR,
                              
                              # Define the mating system (h), and clutch size (k)
                              h = h,
                              k = k)
    
    # Build matrix based on rates specified in the list above
    demographic_matrix <- coucal_matrix(demographic_rates)
    
    # populate sex-specific adult immigration matrix
    immigration_vector <- matrix(data = c(0, demographic_rates$F_Adult_immigration, 
                                          0, demographic_rates$M_Adult_immigration), 
                                 nrow = 4, ncol = 1)
    
    # Determine the ASR at the stable stage distribution
    ASR_SSD <- matrix_ASR(M = demographic_matrix,
                          h = demographic_rates$h,
                          HSR = demographic_rates$HSR, iterations = 1000,
                          num_boot = num_boot,
                          species = species,
                          immigrant_pop_size = immigrant_pop_size,
                          ISR = ISR, iter_add = 1)
    
    # Extract ASR
    ASR_estimate <- ASR_SSD$ASR
    
    # make a list of all the results from this iteration
    bootstrap_results_list <- 
      list(nestling_survival_model_output = nestling_survival_model_output_list, 
           fledgling_survival_model_output = fledgling_survival_model_output_list,
           adult_survival_model_output = adult_survival_model_output_list,
           coucal_vital_rates = coucal_vital_rates,
           ASR_SSD = ASR_SSD,
           bootstrapped_data = coucal_boot_list)
    
    # extract file directory
    file_directory <- paste0(getwd(), "/output/bootstraps/raw/", bootstrap_name, "_", num_boot, ".Rds")
    
    # save the results of the iteration
    save(bootstrap_results_list, file = file_directory)
    
    bootstrap_results_list
  }