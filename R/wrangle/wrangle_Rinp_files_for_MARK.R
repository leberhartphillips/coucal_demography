## ----RMark export to MARK -----------------------------------------------
################################################################################
### Create the Rinp files for MARK

# function to create the .Rinp file needed to calculate median c-hat in MARK
Rinp_prep <- function(proc_file, design_data, directory_name,
                      Phi_structure, p_structure, project_name){
  
  # make a general model (i.e. most parameterized) and run that model 
  # in MARK through RMark
  general.model <- mark(data = proc_file, ddl = design_data,
                        model.parameters = list(Phi = Phi_structure,
                                                p = p_structure))
  
  # specify the directory where the results should go
  setwd(directory_name)
  
  # export the results as a .Rinp file
  export.MARK(x = proc_file, project.name = project_name, 
              general.model, replace = TRUE)
}

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

# CJS model ####
BC_Burnham_global <- 
  RMark::mark(data = BC_coucal_fledgling.proc,
              ddl = BC_coucal_fledgling.ddl, 
              begin.time = 2001, model = "Burnham", groups = "sex",
              model.parameters = list(Phi = list(formula = ~ sex), 
                                      p = list(formula =  ~ 1)),
              wrap = FALSE, threads = 1, brief = TRUE,
              silent = TRUE, output = FALSE, delete = TRUE)

BC_CJS_global <- 
  RMark::mark(model = "CJS", groups = "sex", 
              # model.name = "BC_Burnham_global", filename = "BC_Burnham_global",
              data = BC_coucal_fledgling_CJS.proc, 
              ddl = BC_coucal_fledgling_CJS.ddl,
              model.parameters = list(Phi = list(formula = ~sex * Time),
                                      p = list(formula = ~Time)))
              # model.parameters = list(S = list(formula = ~sex * Time),
              #                         p = list(formula = ~Time),
              #                         r = list(formula = ~Time),
              #                         F = list(formula = ~Time)))

BC_Burnham_global <- 
  RMark::mark(model = "Burnham", groups = "sex", 
              # model.name = "BC_Burnham_global", filename = "BC_Burnham_global",
              data = BC_coucal_fledgling.proc, 
              ddl = BC_coucal_fledgling.ddl,
              model.parameters = list(S = list(formula = ~sex * Time),
                                      p = list(formula = ~Time),
                                      r = list(formula = ~1),
                                      F = list(formula = ~Time)))

BC_Burnham_top <- 
  RMark::mark(model = "Burnham", groups = "sex", 
              # model.name = "BC_Burnham_global", filename = "BC_Burnham_global",
              data = BC_coucal_fledgling.proc, 
              ddl = BC_coucal_fledgling.ddl,
              model.parameters = list(S = list(formula = ~sex * Time),
                                      p = list(formula = ~Time),
                                      r = list(formula = ~1),
                                      F = list(formula = ~Time)))

# specify the directory where the results should go
setwd(directory_name)

# export the results as a .Rinp file
RMark::export.MARK(x = BC_coucal_fledgling_CJS.proc, project.name = "BC_CJS_global_export", 
                   model = BC_CJS_global, replace = TRUE)
RMark::export.MARK(x = BC_coucal_fledgling.proc, project.name = "BC_Burnham_global_export", 
                   model = BC_Burnham_global, replace = TRUE)
##### Dipper test #####
data(dipper, package = "marked")

dipper.proc <- 
  RMark::process.data(data = dipper,
                      model = "CJS",
                      groups = "sex")
dipper.ddl <- 
  RMark::make.design.data(dipper.proc)

dipper_CJS_global <- 
  RMark::mark(model = "CJS", groups = "sex", 
              # model.name = "BC_Burnham_global", filename = "BC_Burnham_global",
              data = dipper.proc, 
              ddl = dipper.ddl,
              model.parameters = list(Phi = list(formula = ~sex * Time),
                                      p = list(formula = ~Time)))
RMark::export.MARK(x = dipper.proc, project.name = "dipper_test_export", 
                   model = dipper_CJS_global, replace = TRUE)