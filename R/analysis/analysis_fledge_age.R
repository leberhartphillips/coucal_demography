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
  dplyr::select(species, sex, Fledge_age, Fledge_mass, Fledge_tarsus, nest_ID, 
                year, lab_no, Ring_ID, site, brood_size, age_diff, lay_date) %>% 
  filter(site == "Kapunga" & species != "CTC") %>% 
  mutate(Nst_No = str_replace_all(string = nest_ID, fixed(" "), ""),
         Ring_ID = str_replace_all(string = Ring_ID, fixed(" "), ""),
         lab_no = str_replace_all(string = lab_no, fixed(" "), "")) %>% 
  filter(!is.na(Fledge_age)) %>% 
  mutate(Ind_ID = paste(Nst_No, Ring_ID, lab_no, sep = "_")) %>%
  mutate(sex_plot = ifelse(sex == "M", 2.2, 0.8)) %>% 
  mutate(
    # scale the Julian lay date by year
    lay_date_std = scale_by(lay_date ~ year, ., scale = 0)) %>% 
  
  # make the scaled date variable numeric class
  mutate(lay_date_std_num = as.numeric(lay_date_std))

# assess normality of "days_since_fledging" variable
ggplot(data = fledge_dat) +
  geom_histogram(aes(Fledge_age), binwidth = 1) +
  facet_grid(sex ~ species)

# check for repeated measures within nest
fledge_dat %>%
  group_by(nest_ID) %>%
  dplyr::summarise(n_ = n())
 
# Modelling proceedure:
# "Fledge_age" as dependent variable nest_ID and year as random effects.
# Each sex was modeled seperartely to estimate the sex-specific fledge age as an 
# input for the demographic model.
# To control for the effect of body size on age we evaluated the effects of 
# fledge tarsus and fledge mass on fledge age. Fledge mass and fledge tarsus were 
# highly colinear, but fledge tarsus had a better correlation with fledge age and 
# so fledge mass was dropped.
# To control for seasonal effects on fledge age, we evaluated the effect of
# standardized lay day date, however this had virtualy no effect and to make the 
# model parsimonious, we dropped this too.
# Final model included "Fledge_tarsus" as the lone independent variable

#### Modeling (BC) ----
# subset to BC males
male_BC <- filter(fledge_dat, species == "BC" & sex == "M")

# assess collinearity
cor.test(male_BC$Fledge_tarsus, male_BC$Fledge_mass) # highly colinear: 0.63
cor.test(male_BC$Fledge_age, male_BC$Fledge_mass) # 0.24
cor.test(male_BC$Fledge_age, male_BC$Fledge_tarsus) # 0.45 ... better than mass

# run model
mod_fledge_age_BC_male <- 
  lmer(Fledge_age ~ Fledge_tarsus +
         (1 | nest_ID) + (1 | year), 
       data = male_BC)

# get the average size of a male tarsus for the model prediction
new_data_BC_male <- 
  expand.grid(Fledge_tarsus = 
                mean(male_BC$Fledge_tarsus, 
                     na.rm = TRUE))

# simulate prediction
predict_BC_male <- 
  predict(mod_fledge_age_BC_male, 
          newdata = new_data_BC_male, 
          re.form = NA, 
          se.fit = TRUE, 
          nsim = 1000)

# subset to BC females
female_BC <- filter(fledge_dat, species == "BC" & sex == "F")

# assess collinearity
cor.test(female_BC$Fledge_tarsus, female_BC$Fledge_mass) # highly colinear: 0.70
cor.test(female_BC$Fledge_age, female_BC$Fledge_mass) # 0.47
cor.test(female_BC$Fledge_age, female_BC$Fledge_tarsus) # 0.60 ... better than mass

# run model
mod_fledge_age_BC_female <- 
  lmer(Fledge_age ~ Fledge_tarsus +
         (1 | nest_ID) + (1 | year), 
       data = female_BC)

# get the average size of a female tarsus for the model prediction
new_data_BC_female <- 
  expand.grid(Fledge_tarsus = 
                mean(female_BC$Fledge_tarsus, 
                     na.rm = TRUE))

# simulate prediction
predict_BC_female <- 
  predict(mod_fledge_age_BC_female, 
          newdata = new_data_BC_female, 
          re.form = NA, 
          se.fit = TRUE, 
          nsim = 1000)

predict_BC_male$fit[1]
predict_BC_female$fit[1]

#### Modeling (WBC) ----
# subset to WBC males
male_WBC <- filter(fledge_dat, species == "WBC" & sex == "M")

# assess collinearity
cor.test(male_WBC$Fledge_tarsus, male_WBC$Fledge_mass) # highly colinear: 0.63
cor.test(male_WBC$Fledge_age, male_WBC$Fledge_mass) # 0.32
cor.test(male_WBC$Fledge_age, male_WBC$Fledge_tarsus) # 0.53 ... better than mass

# run model
mod_fledge_age_WBC_male <- 
  lmer(Fledge_age ~ Fledge_tarsus +
         (1 | nest_ID) + (1 | year), 
       data = male_WBC)

# get the average size of a male tarsus for the model prediction
new_data_WBC_male <- 
  expand.grid(Fledge_tarsus = 
                mean(male_WBC$Fledge_tarsus, 
                     na.rm = TRUE))

# simulate prediction
predict_WBC_male <- 
  predict(mod_fledge_age_WBC_male, 
          newdata = new_data_WBC_male, 
          re.form = NA, 
          se.fit = TRUE, 
          nsim = 1000)

# subset to BC females
female_WBC <- filter(fledge_dat, species == "WBC" & sex == "F")

# assess collinearity
cor.test(female_WBC$Fledge_tarsus, female_WBC$Fledge_mass) # highly colinear: 0.69
cor.test(female_WBC$Fledge_age, female_WBC$Fledge_mass) # 0.50
cor.test(female_WBC$Fledge_age, female_WBC$Fledge_tarsus) # 0.52 ... better than mass

# run model
mod_fledge_age_WBC_female <- 
  lmer(Fledge_age ~ Fledge_tarsus +
         (1 | nest_ID) + (1 | year), 
       data = female_WBC)

# get the average size of a female tarsus for the model prediction
new_data_WBC_female <- 
  expand.grid(Fledge_tarsus = 
                mean(female_WBC$Fledge_tarsus, 
                     na.rm = TRUE))

# simulate prediction
predict_WBC_female <- 
  predict(mod_fledge_age_WBC_female, 
          newdata = new_data_WBC_female, 
          re.form = NA, 
          se.fit = TRUE, 
          nsim = 1000)

predict_WBC_male$fit[1]
predict_WBC_female$fit[1]

# tidy up estimates
coucal_fledge_age <- 
  data.frame(trait = c("fledge_age"),
             species = c("BC", "BC", "WBC", "WBC"),
             sex = c("F", "M", "F", "M"),
             mean = c(predict_BC_female$fit[1],
                      predict_BC_male$fit[1],
                      predict_WBC_female$fit[1],
                      predict_WBC_male$fit[1]),
             CI_low = c(predict_BC_female$ci.fit[1, 1],
                        predict_BC_male$ci.fit[1, 1],
                        predict_WBC_female$ci.fit[1, 1],
                        predict_WBC_male$ci.fit[1, 1]),
             CI_high = c(predict_BC_female$ci.fit[2, 1],
                         predict_BC_male$ci.fit[2, 1],
                         predict_WBC_female$ci.fit[2, 1],
                         predict_WBC_male$ci.fit[2, 1]),
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