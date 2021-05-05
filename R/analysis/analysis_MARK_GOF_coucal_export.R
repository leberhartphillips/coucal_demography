# load packages
source("scripts/01_libraries.R")

# load capture histories
data.sources = list.files(path = "cooked_data", 
                          pattern="*ch.rds$", full.names = TRUE, 
                          ignore.case = TRUE)
sapply(data.sources, load, .GlobalEnv)

BC_fledgling_ch <- Black_Coucal_fledgling_Burnham_ch

fledgling_boot <- 
  fledgling %>% 
  group_by(nest_ID) %>%
  sample_n(1)

write_csv(BC_fledgling_ch, file = "cooked_data/BC_fledgling_ch.csv")

# process the fledgling data as a "Burnham"  analysis
BC_coucal_fledgling.proc <- RMark::process.data(data = BC_fledgling_ch,
                                                model = "Burnham",
                                                groups = c("sex"))

# make design matrix
BC_coucal_fledgling.ddl <- 
  RMark::make.design.data(BC_coucal_fledgling.proc)

# run global model
BC_Burnham_global <- 
  RMark::mark(model = "Burnham", groups = "sex", 
              model.name = "BC_Burnham_global", filename = "BC_Burnham_global",
              data = BC_coucal_fledgling.proc, 
              ddl = BC_coucal_fledgling.ddl,
              model.parameters = list(S = list(formula = ~sex * Time),
                                      p = list(formula = ~Time),
                                      r = list(formula = ~1),
                                      F = list(formula = ~Time)))

# Export as a file to be read by Program MARK
RMark::export.MARK(x = BC_coucal_fledgling.proc, project.name = "BC_Burnham_global", 
                   model = BC_Burnham_global, replace = TRUE)