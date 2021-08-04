source("R/project/project_libraries.R")
source("R/project/project_plotting.R")

# source("scripts/02_func_matrix_ASR().R")
# source("scripts/02_func_coucal_matrix().R")

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

# load output
load("../coucal_ASR/coucal_ASR/output/wrangled/White-browed_Coucal_bootstrap_result_clean.rds")

CI <- 0.95

sex_specific_adult_survival <- 
  WBC_boot_out$survival_rates_boot %>% 
  filter(stage == "adult" & rate == "survival") %>% 
  group_by(sex) %>% 
  dplyr::summarise(mean_annual_survival = mean(value, na.rm = TRUE),
                   var_annual_survival = var(value, na.rm = TRUE),
                   median_annual_surival = median(value, na.rm = TRUE),
                   sd_annual_survival = sd(value, na.rm = TRUE),
                   lcl_annual_survival = stats::quantile(value, (1 - CI)/2, na.rm = TRUE),
                   ucl_annual_survival = stats::quantile(value, 1 - (1 - CI)/2, na.rm = TRUE),
                   max_annual_survival = max(value),
                   min_annual_survival = min(value)) %>% 
  mutate(sex_plot = ifelse(sex == "Male", 1.8, 1.2),
         sex_lab = ifelse(sex == "Male", 1, 2))

WBC_boot_out$survival_rates_boot %>% 
  filter(stage == "adult" & rate == "survival") %>% 
  ggplot() +
  geom_jitter(aes(y = value, x = sex, 
                  color = sex),
              width = 0.1, height = 0.1,
              alpha = 0.75, 
              fill = "white", shape = 16) +
  geom_pointrange(data = sex_specific_adult_survival, 
                  aes(y = median_annual_surival, x = sex_plot, 
                      ymin = (lcl_annual_survival), 
                      ymax = (ucl_annual_survival),
                      fill = sex),
                  size = 0.8, shape = 21, color = "black") +
  luke_theme +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  ylab("Annual adult survival rate (± 95% CI)") +
  scale_color_manual(values = plot_palette_sex) +
  scale_fill_manual(values = plot_palette_sex)

bootstrap_adult_survival_CJS <- 
  function(coucal_boot_list, num_boot, iter_add,
           # egg_survival,
           # ISR,
           # immigrant_pop_size = 100,
           # HSR = 0.5,
           # h,
           # k,
           # flight_age,
           first_year,
           bootstrap_name,
           species,
           prefix_number) {
    
    # specify the bootstrapped data samples (from the previous function)
    adult_ch <- coucal_boot_list[["adult_boot"]]
    
    # clean up capture histories
    adult_ch <-    
      adult_ch %>% 
      ungroup() %>% 
      as.data.frame() %>% 
      mutate(sex = as.factor(sex))
    
    # process the adult data as a "CJS"  analysis
    coucal_adult.proc <- RMark::process.data(adult_ch, 
                                             model = "CJS",
                                             groups = c("sex"),
                                             begin.time = first_year)
    
    # Create the design matrix from the processed mark-recapture datasets
    coucal_adult.ddl <- RMark::make.design.data(coucal_adult.proc)
    
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
    adult_survival_analysis_out <-
      adult_survival_analysis_run(proc_data = coucal_adult.proc,
                                  design_data = coucal_adult.ddl)
    
    extract_top_model_output <- 
      function(rmark_output, top_model = TRUE, mod_num){
        # Find the model number for the first ranked model of the AIC table
        if(top_model == TRUE){
          mod_num <- 
            as.numeric(rownames(rmark_output$model.table[1,]))
        }
        
        else{
          mod_num <- mod_num
        }
        
        # extract and wrangle reals from model output 
        reals <- 
            rmark_output[[mod_num]]$results$real %>% 
            bind_cols(data.frame(str_split_fixed(rownames(.), " ", n = 5)), .) %>% 
            filter(X1 == "Phi") %>% 
            mutate(sex = as.factor(ifelse(unlist(str_extract_all(X2,"[FM]")) == "F", "Female", "Male"))) %>% 
            select(sex, estimate) %>% 
            mutate(iter = num_boot + ((iter_add - 1) * niter),
                   species = species)
        
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
    
    # extract and format survival rates from model output
    adult_survival_model_output_list <- 
      extract_top_model_output(rmark_output = adult_survival_analysis_out, 
                               stage_name = "adult")
    
    # make a list of all the results from this iteration
    bootstrap_results_list <- 
      list(adult_survival_model_output = adult_survival_model_output_list,
           bootstrapped_data = coucal_boot_list)
    
    # extract file directory
    file_directory <- paste0(getwd(), "/output/bootstraps/raw/", bootstrap_name, "_", num_boot, ".Rds")
    
    # save the results of the iteration
    save(bootstrap_results_list, file = file_directory)
    
    bootstrap_results_list
  }
