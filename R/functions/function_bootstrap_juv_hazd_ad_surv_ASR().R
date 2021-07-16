source("R/project/project_libraries.R")
source("R/functions/function_matrix_ASR().R")
source("R/functions/function_coucal_matrix().R")

# bootstrap_hazard_ASR() extracts the stage and sex-specific hazard rates from the 
# bootstrapped sample and feeds them into the matrix model to estimates the ASR.
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

bootstrap_juv_hazd_ad_surv_ASR <- 
  function(coucal_boot_list, # result from bootstrap_hazard_data()
           num_boot, # bootstrap number in the loop (do not specify)
           iter_add, # addition to bootstrap number in the loop (do not specify)
           egg_surv_dist = c(1, 0), # mean and sd of egg survival
           ISR_dist = c(0.5, 0), # mean and sd of immigrant ASR
           immigrant_pop_size = 100, # starting population of simulation
           HSR_dist = c(0.5, 0), # mean and sd of hatching sex ratio
           h_dist = c(1, 0), # mean and sd of mating system
           k_dist = c(4, 0), # mean and sd of clutch size
           fledge_age_distF, # mean and sd of female fledging age
           flight_age_distF, # mean and sd of female first flight age
           fledge_age_distM, # mean and sd of male fledging age
           flight_age_distM, # mean and sd of male first flight age
           # first_year, 
           bootstrap_name, # ID of bootstrap iteration (do not specify)
           species, # coucal species name
           prefix_number, # prefix to give outfile (to prevent overwriting)
           alpha_value = 1.4, # alpha for hazard curve spline estimation
           adult_surival_boot_out) { # vector of 1000 bootstraps of adult survival rates
    
    # specify the bootstrapped data samples (from the previous function)
    offspring_data <- coucal_boot_list[["offspring_boot"]]
    
    # clean up capture histories
    offspring_data <-    
      offspring_data %>% 
      ungroup() %>% 
      as.data.frame()
    
    # fit smoothed spline of hazard function for either sex
    M_haz_ss <- sshzd(Surv(exit, event, entry) ~ exit, 
                      data = filter(offspring_data, sex == "M"), 
                      alpha = alpha_value)

    F_haz_ss <- sshzd(Surv(exit, event, entry) ~ exit, 
                      data = filter(offspring_data, sex == "F"), 
                      alpha = alpha_value)

    haz_ss_lambda <- list(Male_haz_ss_lambda = M_haz_ss$lambda,
                          Female_haz_ss_lambda = F_haz_ss$lambda)
    
    # extract fitted estimates from the spline function
    M_haz_ss_curve <- 
      hzdcurve.sshzd(object = M_haz_ss, time = coucal_boot_list[["time_vector"]], se = TRUE)
    
    F_haz_ss_curve <- 
      hzdcurve.sshzd(object = F_haz_ss, time = coucal_boot_list[["time_vector"]], se = TRUE)
    
    haz_ss_curve <- 
      expand.grid(species = species, 
                  age = coucal_boot_list[["time_vector"]],
                  sex = c("Male", "Female")) %>% 
      mutate(fit = c(M_haz_ss_curve$fit, F_haz_ss_curve$fit),
             se = c(M_haz_ss_curve$se, F_haz_ss_curve$se)) %>% 
      mutate(estimate = 1 - fit,
             upper = 1 - fit * exp(1.96 * se),
             lower = 1 - fit / exp(1.96 * se),
             iter = num_boot)
    
    # # process the adult data as a "CJS"  analysis
    # coucal_adult.proc <- RMark::process.data(adult_ch, 
    #                                          model = "CJS",
    #                                          groups = c("sex"),
    #                                          begin.time = first_year)
    # 
    # # Create the design matrix from the processed mark-recapture datasets
    # coucal_adult.ddl <- RMark::make.design.data(coucal_adult.proc)
    # 
    # adult_survival_analysis_run = 
    #   function(proc_data, design_data){
    #     # sex- and stage-specific survival:
    #     Phi.sex = list(formula = ~ sex) 
    #     
    #     # Models exploring variation in encounter probability
    #     # constant:
    #     p.dot = list(formula =  ~ 1)
    #     
    #     # sex-dependent:
    #     p.sex = list(formula =  ~ sex)
    #     
    #     # factorial variation across year:
    #     p.Time = list(formula =  ~ time)
    #     
    #     # additive effects of sex and factorial year:
    #     p.sex_Time = list(formula =  ~ sex + time)
    #     
    #     # create a list of candidate models for all the a models above that begin with 
    #     # either "Phi." or "p."
    #     cml <-  RMark::create.model.list("CJS")
    #     
    #     # specify the data, design matrix, delete unneeded output files, and 
    #     # run the models in Program MARK
    #     model.list <-  RMark::mark.wrapper(cml, data = proc_data, 
    #                                        ddl = design_data, delete = TRUE, 
    #                                        wrap = FALSE, threads = 1, brief = TRUE,
    #                                        silent = TRUE, output = FALSE, prefix = prefix_number)
    #     
    #     # output the model list and sotre the results
    #     return(model.list)
    #   }
    # 
    # # Run the models on the bootstrapped data
    # adult_survival_analysis_out <-
    #   adult_survival_analysis_run(proc_data = coucal_adult.proc,
    #                               design_data = coucal_adult.ddl)
    # 
    # extract_top_model_output <- 
    #   function(rmark_output, stage_name, top_model = TRUE, mod_num){
    #     # Find the model number for the first ranked model of the AIC table
    #     if(top_model == TRUE){
    #       mod_num <- 
    #         as.numeric(rownames(rmark_output$model.table[1,]))
    #     }
    #     
    #     else{
    #       mod_num <- mod_num
    #     }
    #     
    #     # extract and wrangle reals from model output 
    #     reals <- 
    #         rmark_output[[mod_num]]$results$real %>% 
    #         bind_cols(data.frame(str_split_fixed(rownames(.), " ", n = 5)), .) %>% 
    #         filter(X1 == "Phi") %>% 
    #         mutate(sex = as.factor(ifelse(unlist(str_extract_all(X2,"[FM]")) == "F", "Female", "Male"))) %>% 
    #         dplyr::select(sex, estimate) %>% 
    #         mutate(iter = num_boot + ((iter_add - 1) * niter),
    #                species = species)
    #     
    #     # extract and wrangle the beta estimates from the model output
    #     betas <- 
    #       rmark_output[[mod_num]]$results$beta %>% 
    #       mutate(statistic = rownames(.),
    #              iter = num_boot + ((iter_add - 1) * niter),
    #              species = species) %>% 
    #       mutate(parameter = ifelse(grepl(x = statistic, pattern = "p:"), "p",
    #                                 ifelse(grepl(x = statistic, pattern = "Phi:"), "Phi", "XXX")),
    #              variable = ifelse(grepl(x = statistic, pattern = "Intercept"), "Intercept", 
    #                                ifelse(grepl(x = statistic, pattern = "sexM"), "sexM", 
    #                                       ifelse(grepl(x = statistic, pattern = "Time"), "Time","XXX")))) %>% 
    #       dplyr::select(-statistic)
    #     
    #     # consolidate the AIC model selection results
    #     AIC_table <- 
    #       rmark_output$model.table %>% 
    #       mutate(iter = num_boot + ((iter_add - 1) * niter),
    #              species = species,
    #              model_no_orig = as.numeric(rownames(.))) %>% 
    #       mutate(model_no_rank = as.numeric(rownames(.)))
    #     
    #     # consolidate the all output into a list
    #     survival_model_output_list <- 
    #       list(reals = reals,
    #            betas = betas,
    #            AIC_table = AIC_table)
    #     
    #     survival_model_output_list
    #   }
    # 
    # # extract and format survival rates from juvenile and adult model output
    # adult_survival_model_output_list <- 
    #   extract_top_model_output(rmark_output = adult_survival_analysis_out)
    
    # transform the daily nestling survival (DCS) to apparent fledgling success
    # by calculating the product of all DCS estimates:
    
    fledge_ageM <- rnorm(1, fledge_age_distM[1], fledge_age_distM[2])
    fledge_ageF <- rnorm(1, fledge_age_distF[1], fledge_age_distF[2])
    
    flight_ageM <- rnorm(1, flight_age_distM[1], flight_age_distM[2])
    flight_ageF <- rnorm(1, flight_age_distF[1], flight_age_distF[2])
    
    coucal_nestling_survivalF <-
      haz_ss_curve %>% 
      filter(age <= fledge_ageF & sex == "Female") %>% 
      group_by(sex) %>%
      dplyr::summarise(value = prod(estimate), .groups = 'drop') %>%
      mutate(stage = "nestling",
             rate = "survival")
    
    coucal_nestling_survivalM <-
      haz_ss_curve %>% 
      filter(age <= fledge_ageM & sex == "Male") %>% 
      group_by(sex) %>%
      dplyr::summarise(value = prod(estimate), .groups = 'drop') %>%
      mutate(stage = "nestling",
             rate = "survival")
    
    coucal_nestling_survival <- 
      bind_rows(coucal_nestling_survivalF, coucal_nestling_survivalM)
    
    coucal_groundling_survivalF <-
      haz_ss_curve %>% 
      filter(age < (flight_ageF + fledge_ageF) & age > fledge_ageF & sex == "Female") %>% 
      group_by(sex) %>% 
      dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
      mutate(stage = "groundling",
             rate = "survival")
    
    coucal_groundling_survivalM <-
      haz_ss_curve %>% 
      filter(age < (flight_ageM + fledge_ageM) & age > fledge_ageM & sex == "Male") %>% 
      group_by(sex) %>% 
      dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
      mutate(stage = "groundling",
             rate = "survival")
    
    coucal_groundling_survival <- 
      bind_rows(coucal_groundling_survivalF, coucal_groundling_survivalM)
    
    coucal_fledgling_survivalF <-
      haz_ss_curve %>% 
      filter(age >= (flight_ageF + fledge_ageF) & sex == "Female") %>% 
      group_by(sex) %>% 
      dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
      mutate(stage = "fledgling",
             rate = "survival")
    
    coucal_fledgling_survivalM <-
      haz_ss_curve %>% 
      filter(age >= (flight_ageM + fledge_ageM) & sex == "Male") %>% 
      group_by(sex) %>% 
      dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
      mutate(stage = "fledgling",
             rate = "survival")
    
    coucal_fledgling_survival <- 
      bind_rows(coucal_fledgling_survivalF, coucal_fledgling_survivalM)
    
    coucal_adult_survival <-
      filter(adult_surival_boot_out, iter == sample(1:1000, 1)) %>%
      dplyr::select(-species, -iter) %>%
      mutate(stage = "adult",
             rate = "survival")
    
    # coucal_adult_survival <-
    #   adult_survival_model_output_list$reals %>%
    #   dplyr::select(-species, -iter) %>% 
    #   mutate(stage = "adult",
    #          rate = "survival") %>% 
    #   dplyr::rename(value = estimate)
    
    ISR <- rnorm(1, ISR_dist[1], ISR_dist[2])
    
    adult_F_immigrants <- immigrant_pop_size * (1 - ISR)
    adult_M_immigrants <- immigrant_pop_size * ISR
    
    coucal_adult_immigration <- 
      data.frame(sex = c("Female", "Male"),
                 value = c(adult_F_immigrants, adult_M_immigrants),
                 stage = c("adult"),
                 rate = c("immigration"))
    
    egg_survival <- rnorm(1, egg_surv_dist[1], egg_surv_dist[2])
    
    coucal_egg_survival <- 
      data.frame(sex = NA,
                 value = egg_survival,
                 stage = c("egg"),
                 rate = c("survival"))
    
    coucal_mating_system <- 
      data.frame(sex = NA,
                 value = 1/rnorm(1, h_dist[1], h_dist[2]),
                 stage = c("h"),
                 rate = c("fecundity"))
    
    coucal_clutch_size <- 
      data.frame(sex = NA,
                 value = rnorm(1, k_dist[1], k_dist[2]),
                 stage = c("k"),
                 rate = c("fecundity"))
    
    coucal_HSR <- 
      data.frame(sex = NA,
                 value = rnorm(1, HSR_dist[1], HSR_dist[2]),
                 stage = c("HSR"),
                 rate = c("fecundity"))
    
    coucal_flight_age <- 
      data.frame(sex = c("Female", "Male"),
                 value = c(flight_ageF, flight_ageM),
                 stage = c("flight_age"),
                 rate = c("development"))
    
    coucal_fledge_age <- 
      data.frame(sex = c("Female", "Male"),
                 value = c(fledge_ageF, fledge_ageM),
                 stage = c("fledge_age"),
                 rate = c("development"))
    
    coucal_ISR <- 
      data.frame(sex = NA,
                 value = c(ISR),
                 stage = c("ISR"),
                 rate = c("immigration"))
    
    # Bind the juvenile and adult dataframe with the nestlings
    coucal_vital_rates <- 
      bind_rows(coucal_egg_survival,
                coucal_nestling_survival,
                coucal_groundling_survival,
                coucal_fledgling_survival,
                coucal_adult_survival,
                coucal_adult_immigration,
                coucal_mating_system,
                coucal_clutch_size,
                coucal_HSR,
                coucal_flight_age,
                coucal_fledge_age,
                coucal_ISR) %>% 
      as.data.frame() %>% 
      mutate(iter = num_boot, #+ ((iter_add - 1) * niter),
             species = species) %>% 
      arrange(sex, rate, stage)
    
    # Create a list of demographic rates from the survival analyses above
    demographic_rates <- list(Egg_survival = pull(filter(coucal_vital_rates, stage == "egg"), value),
                              F_Nestling_survival = pull(filter(coucal_vital_rates, 
                                                                stage == "nestling" & rate == "survival" & sex == "Female"), 
                                                         value),
                              F_Groundling_survival = pull(filter(coucal_vital_rates, 
                                                                  stage == "groundling" & rate == "survival" & sex == "Female"), 
                                                           value),
                              F_Fledgling_survival = pull(filter(coucal_vital_rates, 
                                                                 stage == "fledgling" & rate == "survival" & sex == "Female"), 
                                                          value),
                              F_Adult_survival = pull(filter(coucal_vital_rates, 
                                                             stage == "adult" & rate == "survival" & sex == "Female"), 
                                                      value),
                              F_Adult_immigration = pull(filter(coucal_vital_rates, 
                                                                stage == "adult" & rate == "immigration" & sex == "Female"), 
                                                         value),
                              M_Nestling_survival = pull(filter(coucal_vital_rates, 
                                                                stage == "nestling" & rate == "survival" & sex == "Male"), 
                                                         value),
                              M_Groundling_survival = pull(filter(coucal_vital_rates, 
                                                                  stage == "groundling" & rate == "survival" & sex == "Male"), 
                                                           value),
                              M_Fledgling_survival = pull(filter(coucal_vital_rates, 
                                                                 stage == "fledgling" & rate == "survival" & sex == "Male"), 
                                                          value),
                              M_Adult_survival = pull(filter(coucal_vital_rates, 
                                                             stage == "adult" & rate == "survival" & sex == "Male"), 
                                                      value),
                              M_Adult_immigration = pull(filter(coucal_vital_rates, 
                                                                stage == "adult" & rate == "immigration" & sex == "Male"), 
                                                         value),
                              
                              # Define hatching sex ratio
                              HSR = pull(filter(coucal_vital_rates, stage == "HSR"), value),
                              
                              # Define the mating system (h), and clutch size (k)
                              h = pull(filter(coucal_vital_rates, stage == "h"), value),
                              k = pull(filter(coucal_vital_rates, stage == "k"), value))
    
    # Build matrix based on rates specified in the list above
    demographic_matrix <- coucal_matrix(demographic_rates)
    
    # populate sex-specific adult immigration matrix
    immigration_vector <- matrix(data = c(0, demographic_rates$F_Adult_immigration, 
                                          0, demographic_rates$M_Adult_immigration), 
                                 nrow = 4, ncol = 1)
    
    # Determine the ASR at the stable stage distribution
    ASR_SSD <- matrix_ASR(M = demographic_matrix,
                          h = demographic_rates$h,
                          HSR = demographic_rates$HSR, 
                          iterations = 100,
                          num_boot = num_boot,
                          species = species,
                          immigrant_pop_size = immigrant_pop_size,
                          ISR = ISR, 
                          iter_add = 1)
    
    # Extract ASR
    ASR_estimate <- ASR_SSD$ASR
    
    # make a list of all the results from this iteration
    bootstrap_results_list <- 
      list(offspring_hazard_lambda = haz_ss_lambda, 
           offspring_hazard_rates = haz_ss_curve,
           coucal_vital_rates = coucal_vital_rates,
           ASR_SSD = ASR_SSD,
           bootstrapped_data = coucal_boot_list)
    
    # extract file directory
    file_directory <- paste0(getwd(), "/output/bootstraps/hazard/raw/", bootstrap_name, "_", num_boot, ".Rds")
    
    # save the results of the iteration
    save(bootstrap_results_list, file = file_directory)
    
    bootstrap_results_list
  }
