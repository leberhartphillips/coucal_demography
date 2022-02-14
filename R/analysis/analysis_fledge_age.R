#### Script Preparation ----
# load packages
source("R/project/project_libraries.R")

# load functions
function.sources = list.files(path = "R/functions",
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
try (sapply(function.sources, source), silent = TRUE)

#### Data Wrangle ----
# load raw data
fledge_dat <- 
  read.csv("data/raw/Coucal_chick_survival_2001-2019_20200129.csv", 
           header = TRUE, stringsAsFactors = FALSE, na.strings = c("", " ", "NA")) %>%
  dplyr::select(species, sex, Fledge_age, nest_ID, year, lab_no, Ring_ID, site) %>% 
  filter(site == "Kapunga" & species != "CTC") %>% 
  mutate(Nst_No = str_replace_all(string = nest_ID, fixed(" "), ""),
         Ring_ID = str_replace_all(string = Ring_ID, fixed(" "), ""),
         lab_no = str_replace_all(string = lab_no, fixed(" "), "")) %>% 
  filter(!is.na(Fledge_age)) %>% 
  mutate(Ind_ID = paste(Nst_No, Ring_ID, lab_no, sep = "_")) %>%
  mutate(sex_plot = ifelse(sex == "M", 2.2, 0.8))

# # assess normality of "days_since_fledging" variable
# ggplot(data = fledge_dat) +
#   geom_histogram(aes(Fledge_age), binwidth = 1) +
#   facet_grid(sex ~ species)

# check for repeated measures within nest
# fledge_dat %>% 
#   group_by(nest_ID) %>% 
#   dplyr::summarise(n_ = n())
 
#### Modeling (BC) ----
# "Fledge_age" as dependent variable, interaction with sex and species
mod_fledge_age_BC <- 
  lmer(Fledge_age ~ sex + 
         (1 | nest_ID) + (1 | year), 
       data = filter(fledge_dat, species == "BC"))

new_data_BC <- 
  expand.grid(sex = c("F","M"))

predict_BC <- 
  predict(mod_fledge_age_BC, 
          newdata = new_data_BC, 
          re.form = NA, 
          se.fit = TRUE, 
          nsim = 1000)

#### Modeling (WBC) ----
# "Fledge_age" as dependent variable, interaction with sex and species
mod_fledge_age_WBC <- 
  lmer(Fledge_age ~ sex + 
         (1 | nest_ID) + (1 | year), 
       data = filter(fledge_dat, species == "WBC"))

new_data_WBC <- 
  expand.grid(sex = c("F","M"))

predict_WBC <- 
  predict(mod_fledge_age_WBC, 
          newdata = new_data_WBC, 
          re.form = NA, 
          se.fit = TRUE, 
          nsim = 1000)

# plot effect sizes of model
# plot(allEffects(fledge_age_mod))

# tidy up estimates
coucal_fledge_age <- 
  data.frame(trait = c("fledge_age"),
             species = c("BC", "BC", "WBC", "WBC"),
             sex = c("F", "M", "F", "M"),
             mean = c(predict_BC$fit[1],
                      predict_BC$fit[2],
                      predict_WBC$fit[1],
                      predict_WBC$fit[2]),
             CI_low = c(predict_BC$ci.fit[1, 1],
                        predict_BC$ci.fit[1, 2],
                        predict_WBC$ci.fit[1, 1],
                        predict_WBC$ci.fit[1, 2]),
             CI_high = c(predict_BC$ci.fit[2, 1],
                         predict_BC$ci.fit[2, 2],
                         predict_WBC$ci.fit[2, 1],
                         predict_WBC$ci.fit[2, 2]),
             n_inds = c(filter(fledge_dat, species == "BC" & sex == "F") %>% 
                          summarise(n_ = n_distinct(Ind_ID)) %>% 
                          pull(n_),
                        filter(fledge_dat, species == "BC" & sex == "M") %>% 
                          summarise(n_ = n_distinct(Ind_ID)) %>% 
                          pull(n_),
                        filter(fledge_dat, species == "WBC" & sex == "F") %>% 
                          summarise(n_ = n_distinct(Ind_ID)) %>% 
                          pull(n_),
                        filter(fledge_dat, species == "BC" & sex == "M") %>% 
                          summarise(n_ = n_distinct(Ind_ID)) %>% 
                          pull(n_)),
             n_nests = c(filter(fledge_dat, species == "BC" & sex == "F") %>% 
                           summarise(n_ = n_distinct(Nst_No)) %>% 
                           pull(n_),
                         filter(fledge_dat, species == "BC" & sex == "M") %>% 
                           summarise(n_ = n_distinct(Nst_No)) %>% 
                           pull(n_),
                         filter(fledge_dat, species == "WBC" & sex == "F") %>% 
                           summarise(n_ = n_distinct(Nst_No)) %>% 
                           pull(n_),
                         filter(fledge_dat, species == "WBC" & sex == "M") %>% 
                           summarise(n_ = n_distinct(Nst_No)) %>% 
                           pull(n_)),
             n_years = c(filter(fledge_dat, species == "BC" & sex == "F") %>% 
                           summarise(n_ = n_distinct(year)) %>% 
                           pull(n_),
                         filter(fledge_dat, species == "BC" & sex == "M") %>% 
                           summarise(n_ = n_distinct(year)) %>% 
                           pull(n_),
                         filter(fledge_dat, species == "WBC" & sex == "F") %>% 
                           summarise(n_ = n_distinct(year)) %>% 
                           pull(n_),
                         filter(fledge_dat, species == "WBC" & sex == "M") %>% 
                           summarise(n_ = n_distinct(year)) %>% 
                           pull(n_))) %>% 
  mutate(sd = ifelse(!is.na(CI_low), 
                         approx_sd(x1 = CI_low, x2 = CI_high),
                         CI_low))

rm(fledge_dat, new_data_WBC, new_data_BC, predict_WBC, predict_BC, 
   mod_fledge_age_WBC, mod_fledge_age_BC)