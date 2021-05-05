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

BC_nestling_ch <- Black_Coucal_nestling_Known_ch
WBC_nestling_ch <- White_browed_Coucal_nestling_Known_ch

BC_fledgling_ch <- Black_Coucal_fledgling_Burnham_ch
WBC_fledgling_ch <- White_browed_Coucal_fledgling_Burnham_ch

Black_Coucal_fledgling_CJS_ch <- 
  Black_Coucal_fledgling_Burnham_ch %>% 
  mutate(ch = gsub("(.).", "\\1", ch))

prefix_number <- "test"

Black_Coucal_fledgling_CJS_ch

# process the nestling data as a "Known" fate analysis
BC_coucal_nestling.proc <- RMark::process.data(data = BC_nestling_ch,
                                     model = "Known",
                                     groups = c("sex"))

# make design matrix
BC_coucal_nestling.ddl <- 
  RMark::make.design.data(BC_coucal_nestling.proc)

# process the nestling data as a "Known" fate analysis
WBC_coucal_nestling.proc <- RMark::process.data(data = WBC_nestling_ch,
                                        model = "Known",
                                        groups = c("sex"))

# make design matrix
WBC_coucal_nestling.ddl <- 
  RMark::make.design.data(WBC_coucal_nestling.proc)

# process the fledgling data as a "Burnham"  analysis
BC_coucal_fledgling.proc <- RMark::process.data(data = BC_fledgling_ch,
                                      model = "Burnham",
                                      groups = c("sex"))

release.gof(BC_coucal_fledgling.proc)

# make design matrix
BC_coucal_fledgling.ddl <- 
  RMark::make.design.data(BC_coucal_fledgling.proc)

BC_fledgling_ch %>% 
  mutate(recovery = ifelse(str_detect(ch, "11"), 1, 0)) %>% 
  group_by(recovery) %>% 
  dplyr::summarise(n())

# process the fledgling data as a "Burnham"  analysis
WBC_coucal_fledgling.proc <- RMark::process.data(data = WBC_fledgling_ch,
                                         model = "Burnham",
                                         groups = c("sex"))

# make design matrix
WBC_coucal_fledgling.ddl <- 
  RMark::make.design.data(WBC_coucal_fledgling.proc)

WBC_fledgling_ch %>% 
  mutate(recovery = ifelse(str_detect(ch, "11"), 1, 0)) %>% 
  group_by(recovery) %>% 
  dplyr::summarise(n())

# process the fledgling data as a "Burnham"  analysis
BC_coucal_fledgling_CJS.proc <- 
  RMark::process.data(data = Black_Coucal_fledgling_CJS_ch,
                      model = "CJS",
                      groups = c("sex"))

# make design matrix
BC_coucal_fledgling_CJS.ddl <- 
  RMark::make.design.data(BC_coucal_fledgling_CJS.proc)


nestling_survival_analysis_run <- 
  function(proc_data, design_data){
    # apriori model components for
    # S (survival probability):
    
    # sex- and age-specific model
    S.sexXTime = list(formula = ~sex * Time)
    
    # sex- and age-specific model (quadratic)
    S.sexXQuad = list(formula = ~sex * (I(Time) + I(Time^2)))
    
    # sex- and age-specific model (cubic)
    S.sexXCube = list(formula = ~sex * (I(Time) + I(Time^2) + I(Time^3)))
    
    # create model list for all a priori models above
    cml <- RMark::create.model.list("Known")
    
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
    # null model
    # S.dot = list(formula = ~1)

    # sex- and age-specific model
    # S.sexXTime = list(formula = ~sex * Time)

    # sex- and age-specific model (quadratic)
    S.sexXQuad = list(formula = ~sex * (I(Time) + I(Time^2)))

    # sex- and age-specific model (cubic)
    # S.sexXCube = list(formula = ~sex * (I(Time) + I(Time^2) + I(Time^3)))
    # 
    # p (encounter probability):
    # # null model
    # p.dot = list(formula = ~1)
    # 
    # # sex-specific model
    # p.sex = list(formula = ~sex)

    # age-specific model
    p.Time = list(formula = ~Time)

    # F (site fidelity probability):
    # # null model
    # F.dot = list(formula = ~1)
    # 
    # # sex-specific model
    # F.sex = list(formula = ~sex)

    # age-specific model
    F.Time = list(formula = ~Time)
    
    # r (recovery probability)
    # null model
    r.dot = list(formula = ~1)

    # # sex-specific model
    # r.sex = list(formula = ~sex)
    # 
    # # age-specific model
    # r.time = list(formula = ~Time)
    
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

fledgling_survival_analysis_run_CJS <-  
  function(proc_data, design_data){
    # null model
    # Phi.dot = list(formula = ~1)
    
    # # sex- and age-specific model
    # Phi.sexXTime = list(formula = ~sex * Time)
    # 
    # # sex- and age-specific model (quadratic)
    Phi.sexXQuad = list(formula = ~sex * (I(Time) + I(Time^2)))
    # 
    # # sex- and age-specific model (cubic)
    # Phi.sexXCube = list(formula = ~sex * (I(Time) + I(Time^2) + I(Time^3)))
    # # sex- and stage-specific survival:
    # Phi.sex = list(formula = ~ sex)
    
    # Phi.sex_Time = list(formula = ~sex + Time)
    # 
    # # # sex- and age-specific model (quadratic)
    # Phi.sex_Quad = list(formula = ~sex + (I(Time) + I(Time^2)))
    # # 
    # # # sex- and age-specific model (cubic)
    # Phi.sex_Cube = list(formula = ~sex + (I(Time) + I(Time^2) + I(Time^3)))
    
    # Models exploring variation in encounter probability
    # constant:
    p.dot = list(formula =  ~ 1)
    # 
    # # sex-dependent:
    p.sex = list(formula =  ~ sex)
    # 
    # # factorial variation across year:
    p.Time = list(formula =  ~ Time)
    # 
    # # # sex- and age-specific model (quadratic)
    p.sexXQuad = list(formula = ~sex * (I(Time) + I(Time^2)))
    # # 
    # # # sex- and age-specific model (cubic)
    # p.sexXCube = list(formula = ~sex * (I(Time) + I(Time^2) + I(Time^3)))
    # 
    # # # sex- and age-specific model (quadratic)
    # p.sex_Quad = list(formula = ~sex + (I(Time) + I(Time^2)))
    # 
    # # sex- and age-specific model (cubic)
    # p.sex_Cube = list(formula = ~sex + (I(Time) + I(Time^2) + I(Time^3)))
    # 
    # # additive effects of sex and factorial year:
    p.sex_Time = list(formula =  ~ sex + Time)
    # 
    p.sexxTime = list(formula =  ~ sex * Time)
    
    
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

nestling_survival_analysis_out.BC <-
  nestling_survival_analysis_run(proc_data = BC_coucal_nestling.proc,
                                 design_data = BC_coucal_nestling.ddl)

# first determine structure for p with everything else constant: Time
fledgling_survival_analysis_out_BC.p <-
  fledgling_survival_analysis_run(proc_data = BC_coucal_fledgling.proc,
                                  design_data = BC_coucal_fledgling.ddl)

# first determine structure for p with everything else constant: Time
fledgling_survival_analysis_out_BC.CJS <-
  fledgling_survival_analysis_run_CJS(proc_data = BC_coucal_fledgling_CJS.proc,
                                  design_data = BC_coucal_fledgling_CJS.ddl)

# second determine structure for F with S and p constant and p as Time: Time
fledgling_survival_analysis_out_BC.F <-
  fledgling_survival_analysis_run(proc_data = BC_coucal_fledgling.proc,
                                  design_data = BC_coucal_fledgling.ddl)

# third determine structure for r with S as constant and p as Time and F as Time: constant
fledgling_survival_analysis_out_BC.r <-
  fledgling_survival_analysis_run(proc_data = BC_coucal_fledgling.proc,
                                  design_data = BC_coucal_fledgling.ddl)

# fourth determine structure for S with r as constant and p as Time and F as Time: constant
fledgling_survival_analysis_out_BC.S <-
  fledgling_survival_analysis_run(proc_data = BC_coucal_fledgling.proc,
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

head(fledgling_survival_analysis_out_BC.S$model.table)

fledgling_survival_model_output_list_BC <-
  extract_top_model_output(rmark_output = fledgling_survival_analysis_out_BC.S, 
                           stage_name = "fledgling", top_model = FALSE, mod_num = 79)
fledgling_survival_model_output_list_BC$betas

fledgling_survival_model_output_list_BC_CJS <- 
  extract_top_model_output(rmark_output = fledgling_survival_analysis_out_BC.CJS, 
                           stage_name = "fledgling_CJS", top_model = FALSE, mod_num = 2)

nestling_survival_model_output_list_BC <- 
  extract_top_model_output(rmark_output = nestling_survival_analysis_out.BC, 
                           stage_name = "nestling", top_model = TRUE)

fledgling_survival_analysis_out_BC.CJS[[10]]$results$real

ggplot(data = fledgling_survival_model_output_list_BC_CJS$reals) +
  geom_line(aes(x = age, y = estimate, group = sex, color = sex)) +
  geom_ribbon(aes(x = age, ymin = lcl, ymax = ucl, group = sex, fill = sex), alpha = 0.2) +
  theme(legend.position = "top")

ggplot(data = nestling_survival_model_output_list_BC$reals) +
  geom_line(aes(x = age, y = estimate, group = sex, color = sex)) +
  geom_ribbon(aes(x = age, ymin = lcl, ymax = ucl, group = sex, fill = sex), alpha = 0.2) +
  theme(legend.position = "top")

BC_nestling_survival <-
  nestling_survival_model_output_list_BC$reals %>% 
  # filter(age < flight_age) %>% 
  group_by(sex) %>% 
  dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
  mutate(stage = "nestling",
         rate = "survival")

BC_groundling_survival <-
  fledgling_survival_model_output_list_BC$reals %>% 
  filter(age < flight_age) %>% 
  group_by(sex) %>% 
  dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
  mutate(stage = "groundling",
         rate = "survival")

BC_fledgling_survival <-
  fledgling_survival_model_output_list_BC$reals %>% 
  filter(age >= flight_age) %>% 
  group_by(sex) %>% 
  dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
  mutate(stage = "fledgling",
         rate = "survival")


# first determine structure for p with everything else constant: Time
fledgling_survival_analysis_out_WBC.p <-
  fledgling_survival_analysis_run(proc_data = WBC_coucal_fledgling.proc,
                                  design_data = WBC_coucal_fledgling.ddl)

# second determine structure for F with S and p constant and p as Time: Time
fledgling_survival_analysis_out_WBC.F <-
  fledgling_survival_analysis_run(proc_data = WBC_coucal_fledgling.proc,
                                  design_data = WBC_coucal_fledgling.ddl)

# third determine structure for r with S as constant and p as Time and F as Time: constant
fledgling_survival_analysis_out_WBC.r <-
  fledgling_survival_analysis_run(proc_data = WBC_coucal_fledgling.proc,
                                  design_data = WBC_coucal_fledgling.ddl)

# fourth determine structure for S with r as constant and p as Time and F as Time: constant
fledgling_survival_analysis_out_WBC.S <-
  fledgling_survival_analysis_run(proc_data = WBC_coucal_fledgling.proc,
                                  design_data = WBC_coucal_fledgling.ddl)
num_boot = 1
iter_add = 1
niter = 1
species = "WBC"
flight_age = 32
fledgling_survival_model_output_list_WBC <- 
  extract_top_model_output(rmark_output = fledgling_survival_analysis_out_WBC.S, 
                           stage_name = "fledgling", top_model = FALSE, mod_num = 4)

ggplot(data = fledgling_survival_model_output_list_WBC$reals) +
  geom_line(aes(x = age, y = estimate, group = sex, color = sex)) +
  geom_ribbon(aes(x = age, ymin = lcl, ymax = ucl, group = sex, fill = sex), alpha = 0.2) +
  theme(legend.position = "top")

WBC_groundling_survival <-
  fledgling_survival_model_output_list_WBC$reals %>% 
  filter(age < flight_age) %>% 
  group_by(sex) %>% 
  dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
  mutate(stage = "groundling",
         rate = "survival")

WBC_fledgling_survival <-
  fledgling_survival_model_output_list_WBC$reals %>% 
  filter(age >= flight_age) %>% 
  group_by(sex) %>% 
  dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
  mutate(stage = "fledgling",
         rate = "survival")