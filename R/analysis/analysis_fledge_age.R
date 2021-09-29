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

# extract model coefficients
mod_fledge_age_BC_coefs <- 
  model_parameters(mod_fledge_age_BC) %>%
  as.data.frame(.)

# plot effect sizes of model
# plot(allEffects(fledge_age_mod))

# tidy up estimates
coucal_fledge_age <- 
  data.frame(trait = c("fledge_age"),
             species = c("BC", "BC", "WBC", "WBC"),
             sex = c("F", "M", "F", "M"),
             mean = c(fledge_age_mod_coefs[1, c(2)],
                                   (fledge_age_mod_coefs[1, c(2)] + fledge_age_mod_coefs[2, c(2)]),
                                   (fledge_age_mod_coefs[1, c(2)] + fledge_age_mod_coefs[3, c(2)]),
                                   (fledge_age_mod_coefs[1, c(2)] + fledge_age_mod_coefs[3, c(2)] + fledge_age_mod_coefs[4, c(2)] + fledge_age_mod_coefs[2, c(2)])),
             CI_low = c(fledge_age_mod_coefs[1, c(5)],
                        (fledge_age_mod_coefs[1, c(5)] + fledge_age_mod_coefs[2, c(5)]),
                        (fledge_age_mod_coefs[1, c(5)] + fledge_age_mod_coefs[3, c(5)]),
                        (fledge_age_mod_coefs[1, c(5)] + fledge_age_mod_coefs[3, c(5)] + fledge_age_mod_coefs[4, c(5)] + fledge_age_mod_coefs[2, c(5)])),
             CI_high = c(fledge_age_mod_coefs[1, c(6)],
                         (fledge_age_mod_coefs[1, c(6)] + fledge_age_mod_coefs[2, c(6)]),
                         (fledge_age_mod_coefs[1, c(6)] + fledge_age_mod_coefs[3, c(6)]),
                         (fledge_age_mod_coefs[1, c(6)] + fledge_age_mod_coefs[3, c(6)] + fledge_age_mod_coefs[4, c(6)] + fledge_age_mod_coefs[2, c(6)])),
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

rm(fledge_dat, fledge_age_mod, fledge_age_mod_coefs)