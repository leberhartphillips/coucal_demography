# load packages
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

set.seed(14051987)
BC_fledgling_ch <- 
  Black_Coucal_fledgling_Burnham_ch %>% 
  group_by(nest_ID) %>%
  sample_n(1) %>% 
  mutate(sex = as.factor(sex)) %>%
  ungroup() %>% 
  dplyr::select(ch, sex)

# process the fledgling data as a "Burnham" analysis
BC_coucal_fledgling.proc <- RMark::process.data(data = BC_fledgling_ch,
                                                model = "Burnham",
                                                groups = c("sex"))

# make design matrix
BC_coucal_fledgling.ddl <- 
  RMark::make.design.data(BC_coucal_fledgling.proc)

burnham_p_analysis_run <- 
  function(proc_data, design_data){
    # apriori model components for
    # S (survival probability):
    # null model
    S.dot = list(formula = ~1)
    
    # p (encounter probability):
    # null model
    p.dot = list(formula = ~1)
    
    # sex-specific model
    p.sex = list(formula = ~sex)
    
    # age-specific model
    p.Time = list(formula = ~Time)
    
    # # sex- and age-specific model
    # p.sexXTime = list(formula = ~sex * Time)
    # 
    # # sex- and age-specific model (quadratic)
    # p.sexXQuad = list(formula = ~sex * (I(Time) + I(Time^2)))
    # 
    # # sex- and age-specific model
    # p.sex_Time = list(formula = ~sex + Time)
    # 
    # # sex- and age-specific model (quadratic)
    # p.sex_Quad = list(formula = ~sex + (I(Time) + I(Time^2)))
    
    # F (site fidelity probability):
    # null model
    F.dot = list(formula = ~1)
    
    # r (recovery probability)
    # null model
    r.dot = list(formula = ~1)
    
    # create model list for all a priori models above
    cml <- RMark::create.model.list("Burnham")
    
    # run the models in program MARK
    model.list <-  RMark::mark.wrapper(cml, data = proc_data, 
                                       ddl = design_data, delete = TRUE, 
                                       wrap = FALSE, threads = 1, brief = TRUE,
                                       silent = TRUE, output = FALSE, prefix = prefix_number)
    
    # output the model list and sotre the results
    return(model.list)
  }

# determine structure p: additive quadratic trend of age and sex
fledgling_p_analysis_out_BC <-
  burnham_p_analysis_run(proc_data = BC_coucal_fledgling.proc,
                         design_data = BC_coucal_fledgling.ddl)

burnham_F_analysis_run <- 
  function(proc_data, design_data){
    # apriori model components for
    # S (survival probability):
    # null model
    S.dot = list(formula = ~1)
    
    # p (encounter probability):
    # sex- and age-specific model (quadratic)
    # p.sex_Quad = list(formula = ~sex + (I(Time) + I(Time^2)))
    p.Time = list(formula = ~Time)
    
    
    # F (site fidelity probability):
    # null model
    F.dot = list(formula = ~1)
    
    # sex-specific model
    F.sex = list(formula = ~sex)
    
    # age-specific model
    F.Time = list(formula = ~Time)
    
    # sex- and age-specific model
    F.sexXTime = list(formula = ~sex * Time)
    
    # sex- and age-specific model (quadratic)
    F.sexXQuad = list(formula = ~sex * (I(Time) + I(Time^2)))
    
    # sex- and age-specific model
    F.sex_Time = list(formula = ~sex + Time)
    
    # sex- and age-specific model (quadratic)
    F.sex_Quad = list(formula = ~sex + (I(Time) + I(Time^2)))
    
    # r (recovery probability)
    # null model
    r.dot = list(formula = ~1)
    
    # create model list for all a priori models above
    cml <- RMark::create.model.list("Burnham")
    
    # run the models in program MARK
    model.list <-  RMark::mark.wrapper(cml, data = proc_data, 
                                       ddl = design_data, delete = TRUE, 
                                       wrap = FALSE, threads = 1, brief = TRUE,
                                       silent = TRUE, output = FALSE, prefix = prefix_number)
    
    # output the model list and sotre the results
    return(model.list)
  }

# determine structure F: linear age trend
fledgling_F_analysis_out_BC <-
  burnham_F_analysis_run(proc_data = BC_coucal_fledgling.proc,
                         design_data = BC_coucal_fledgling.ddl)

burnham_r_analysis_run <- 
  function(proc_data, design_data){
    # apriori model components for
    # S (survival probability):
    # null model
    S.dot = list(formula = ~1)
    
    # p (encounter probability):
    # sex- and age-specific model (quadratic)
    # p.sex_Quad = list(formula = ~sex + (I(Time) + I(Time^2)))
    p.Time = list(formula = ~Time)
    
    # F (site fidelity probability):
    # age-specific model
    F.Time = list(formula = ~Time)
    
    # r (recovery probability)
    # null model
    r.dot = list(formula = ~1)
    
    # sex-specific model
    r.sex = list(formula = ~sex)
    
    # age-specific model
    r.Time = list(formula = ~Time)
    
    # sex- and age-specific model
    r.sexXTime = list(formula = ~sex * Time)
    
    # sex- and age-specific model (quadratic)
    r.sexXQuad = list(formula = ~sex * (I(Time) + I(Time^2)))
    
    # sex- and age-specific model
    r.sex_Time = list(formula = ~sex + Time)
    
    # sex- and age-specific model (quadratic)
    r.sex_Quad = list(formula = ~sex + (I(Time) + I(Time^2)))
    
    # create model list for all a priori models above
    cml <- RMark::create.model.list("Burnham")
    
    # run the models in program MARK
    model.list <-  RMark::mark.wrapper(cml, data = proc_data, 
                                       ddl = design_data, delete = TRUE, 
                                       wrap = FALSE, threads = 1, brief = TRUE,
                                       silent = TRUE, output = FALSE, prefix = prefix_number)
    
    # output the model list and sotre the results
    return(model.list)
  }

# determine structure r: constant
fledgling_r_analysis_out_BC <-
  burnham_r_analysis_run(proc_data = BC_coucal_fledgling.proc,
                         design_data = BC_coucal_fledgling.ddl)


burnham_S_analysis_run <- 
  function(proc_data, design_data){
    # apriori model components for
    # S (survival probability):
    # null model
    S.dot = list(formula = ~1)
    
    # sex- and age-specific model
    S.sexXTime = list(formula = ~sex * Time)
    
    # sex- and age-specific model (quadratic)
    S.sexXQuad = list(formula = ~sex * (I(Time) + I(Time^2)))
    
    # sex- and age-specific model
    S.sex_Time = list(formula = ~sex + Time)
    
    # sex- and age-specific model (quadratic)
    S.sex_Quad = list(formula = ~sex + (I(Time) + I(Time^2)))
    
    # p (encounter probability):
    # sex- and age-specific model (quadratic)
    # p.sex_Quad = list(formula = ~sex + (I(Time) + I(Time^2)))
    p.Time = list(formula = ~Time)
    
    # F (site fidelity probability):
    # age-specific model
    F.Time = list(formula = ~Time)
    
    # r (recovery probability)
    # null model
    r.dot = list(formula = ~1)
    
    # create model list for all a priori models above
    cml <- RMark::create.model.list("Burnham")
    
    # run the models in program MARK
    model.list <-  RMark::mark.wrapper(cml, data = proc_data, 
                                       ddl = design_data, delete = TRUE, 
                                       wrap = FALSE, threads = 1, brief = TRUE,
                                       silent = TRUE, output = FALSE, prefix = prefix_number)
    
    # output the model list and sotre the results
    return(model.list)
  }

# determine structure S: linear age trend
fledgling_S_analysis_out_BC <-
  burnham_S_analysis_run(proc_data = BC_coucal_fledgling.proc,
                         design_data = BC_coucal_fledgling.ddl)

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
        select(sex, age, estimate, se, lcl, ucl) %>% 
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
        select(sex, age, estimate, se, lcl, ucl) %>% 
        mutate(iter = num_boot + ((iter_add - 1) * niter),
               species = species)
    }
    else if(stage_name == "fledgling_CJS"){
      reals <- 
        rmark_output[[mod_num]]$results$real %>% 
        bind_cols(data.frame(str_split_fixed(rownames(.), " ", 
                                             n = 5)), .) %>% 
        mutate(age = as.integer(unlist(str_extract_all(X4,"[0-9]+"))),
               sex = as.factor(ifelse(unlist(str_extract_all(X2,"[FM]")) == "F", 
                                      "Female", "Male"))) %>% 
        filter(X1 == "Phi") %>% 
        select(sex, age, estimate, se, lcl, ucl) %>% 
        mutate(iter = num_boot + ((iter_add - 1) * niter),
               species = species)
    }
    else{
      reals <- 
        rmark_output[[mod_num]]$results$real %>% 
        bind_cols(data.frame(str_split_fixed(rownames(.), " ", n = 5)), .) %>% 
        filter(X1 == "Phi") %>% 
        mutate(sex = as.factor(ifelse(unlist(str_extract_all(X2,"[FM]")) == "F", "Female", "Male"))) %>% 
        select(sex, estimate, se, lcl, ucl) %>% 
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

num_boot = 1
iter_add = 1
niter = 1
species = "BC"
flight_age = 36

fledgling_S_model_output_list_BC <-
  extract_top_model_output(rmark_output = fledgling_S_analysis_out_BC, 
                           stage_name = "fledgling", top_model = TRUE)
fledgling_S_model_output_list_BC$betas

ggplot(data = fledgling_S_model_output_list_BC$reals) +
  geom_line(aes(x = age, y = estimate, group = sex, color = sex)) +
  geom_ribbon(aes(x = age, ymin = lcl, ymax = ucl, group = sex, fill = sex), alpha = 0.2) +
  theme(legend.position = "top")
