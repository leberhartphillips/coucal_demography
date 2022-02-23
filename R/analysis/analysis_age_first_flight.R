# load packages
source("R/project/project_libraries.R")

# load functions
function.sources = list.files(path = "R/functions",
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
try (sapply(function.sources, source), silent = TRUE)

# load raw data
flight_dat <- 
  read.csv("data/raw/coucals_first_flights_data_per_individual2020.csv", 
           header = TRUE, stringsAsFactors = FALSE, na.strings = c("", " ", "NA")) %>% 
  mutate(Nst_No = str_replace_all(string = Nst_No, fixed(" "), ""),
         Ind_ID = str_replace_all(string = Ind_ID, fixed(" "), "")) %>%
  dplyr::rename(sex = Sex,
                nest_ID = Nst_No,
                year = Year,
                species = Spp,
                flight_age = days_since_fledging,
                Ring_ID = Ind_ID) %>% 
  mutate(sex = ifelse(sex == "Male", "M", ifelse(sex == "Female", "F", sex))) %>% 
  mutate(sex_plot = ifelse(sex == "M", 2.2, 0.8))

# assess normality of "days_since_fledging" variable
ggplot(data = flight_dat) +
  geom_histogram(aes(flight_age), binwidth = 3) +
  facet_grid(sex ~ species)

# check for repeated measures within nest
flight_dat %>%
  group_by(nest_ID) %>%
  dplyr::summarise(n_ = n())

# see here for a blog about how to bootstrap the model to estimate 
# the 95% confidence interval around the predictions:
# http://www.remkoduursma.com/post/2017-06-15-bootpredictlme4/

# install and load the "bootpredictlme4" package
# devtools::install_github("remkoduursma/bootpredictlme4")
# library(bootpredictlme4)

#### Modeling (WBC) ----
# subset to BC males
male_BC <- filter(flight_dat, species == "BC" & sex == "M")

# assess collinearity
cor.test(male_BC$Fledge_tarsus, male_BC$Fledge_mass) # highly colinear: 0.48
cor.test(male_BC$flight_age, male_BC$Fledge_mass) # 0.24
cor.test(male_BC$flight_age, male_BC$Fledge_tarsus) # 0.45 ... better than mass

mod_flight_age_BC_male <- 
  lmer(flight_age ~ 1 +
         (1 | nest_ID) + (1 | year), 
       data = male_BC)

# simulate prediction
predict_BC_male <- 
  predict(mod_flight_age_BC_male, 
          #newdata = new_data_BC_male, 
          re.form = NA, 
          se.fit = TRUE, 
          nsim = 1000)

summary(mod_flight_age_BC_male)

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
          #newdata = new_data_BC_male, 
          re.form = NA, 
          se.fit = TRUE, 
          nsim = 1000)

# "flight_age" as dependent variable, 
# "sex", "fledge_age", and "fledge_mass" as independent
mod_flight_age_WBC <-
  lmer(flight_age ~ sex + Fledge_age + Fledge_mass +
         (1 | nest_ID) + (1 | year), 
       data = filter(flight_dat, species == "WBC"))

model_parameters(mod_flight_age_WBC)

new_data_WBC <- 
  expand.grid(sex = c("Female","Male"),
              Fledge_age = mean(flight_dat[flight_dat$species == "WBC",]$Fledge_age),
              Fledge_mass = mean(flight_dat[flight_dat$species == "WBC",]$Fledge_mass))

predict_WBC <- 
  predict(mod_flight_age_WBC, 
          newdata = new_data_WBC, 
          re.form = NA, 
          se.fit = TRUE, 
          nsim = 1000)

plot(allEffects(mod_flight_age_WBC))

mod_flight_age_WBC2 <-
  lmer(flight_age ~ sex + 
         (1 | nest_ID) + (1 | year), 
       data = filter(flight_dat, species == "WBC"))

new_data_WBC2 <- 
  expand.grid(sex = c("Female","Male"))
              #Fledge_age = mean(flight_dat[flight_dat$species == "WBC",]$Fledge_age),
              #Fledge_mass = mean(flight_dat[flight_dat$species == "WBC",]$Fledge_mass))

predict_WBC2 <- 
  predict(mod_flight_age_WBC2, 
          newdata = new_data_WBC2, 
          re.form = NA, 
          se.fit = TRUE, 
          nsim = 1000)

# boot_WBC_flight_age <- 
#   bootMer(mod_flight_age_WBC, nsim = 1000, 
#           FUN = function(x) predict(x, newdata = new_data_BC, re.form = NA))

#### Modeling (BC) ----
# "flight_age" as dependent variable, 
# "sex", "fledge_age", and "fledge_mass" as independent
mod_flight_age_BC <- 
  lmer(flight_age ~ sex + Fledge_age + Fledge_mass +
         (1 | nest_ID) + (1 | year), 
       data = filter(flight_dat, species == "BC"))

new_data_BC <- 
  expand.grid(sex = c("Female","Male"),
              Fledge_age = mean(flight_dat[flight_dat$species == "BC",]$Fledge_age),
              Fledge_mass = mean(flight_dat[flight_dat$species == "BC",]$Fledge_mass))

predict_BC <- 
  predict(mod_flight_age_BC, 
          newdata = new_data_BC, 
          re.form = NA, 
          se.fit = TRUE, 
          nsim = 1000)

plot(allEffects(mod_flight_age_BC))

# boot_BC_flight_age <- 
#   bootMer(mod_flight_age_BC, nsim = 1000, 
#           FUN = function(x) predict(x, newdata = new_data_BC, re.form = NA))
# 
# hist(boot_BC_flight_age$t, 
#      breaks=seq(16,26,by=0.5),
#      # ylim=c(0,25),
#      main="", xlab="Reaction time at 5 Days (ms)",
#      col="cornflowerblue")
# box()

# plot effect sizes of model
# plot(allEffects(flight_mod))

# source("R/analysis/analysis_flight_age.R")
# 
# coucal_flight_age

coucal_flight_age <- 
  data.frame(trait = c("flight_age"),
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
             n_inds = c(filter(flight_dat, species == "BC" & sex == "Female") %>% 
                          summarise(n_ = n_distinct(Ring_ID)) %>% 
                          pull(n_),
                        filter(flight_dat, species == "BC" & sex == "Male") %>% 
                          summarise(n_ = n_distinct(Ring_ID)) %>% 
                          pull(n_),
                        filter(flight_dat, species == "WBC" & sex == "Female") %>% 
                          summarise(n_ = n_distinct(Ring_ID)) %>% 
                          pull(n_),
                        filter(flight_dat, species == "BC" & sex == "Male") %>% 
                          summarise(n_ = n_distinct(Ring_ID)) %>% 
                          pull(n_)),
             n_nests = c(filter(flight_dat, species == "BC" & sex == "Female") %>% 
                           summarise(n_ = n_distinct(nest_ID)) %>% 
                           pull(n_),
                         filter(flight_dat, species == "BC" & sex == "Male") %>% 
                           summarise(n_ = n_distinct(nest_ID)) %>% 
                           pull(n_),
                         filter(flight_dat, species == "WBC" & sex == "Female") %>% 
                           summarise(n_ = n_distinct(nest_ID)) %>% 
                           pull(n_),
                         filter(flight_dat, species == "WBC" & sex == "Male") %>% 
                           summarise(n_ = n_distinct(nest_ID)) %>% 
                           pull(n_)),
             n_years = c(filter(flight_dat, species == "BC" & sex == "Female") %>% 
                           summarise(n_ = n_distinct(year)) %>% 
                           pull(n_),
                         filter(flight_dat, species == "BC" & sex == "Male") %>% 
                           summarise(n_ = n_distinct(year)) %>% 
                           pull(n_),
                         filter(flight_dat, species == "WBC" & sex == "Female") %>% 
                           summarise(n_ = n_distinct(year)) %>% 
                           pull(n_),
                         filter(flight_dat, species == "WBC" & sex == "Male") %>% 
                           summarise(n_ = n_distinct(year)) %>% 
                           pull(n_))) %>% 
  mutate(sd = ifelse(!is.na(CI_low), 
                     approx_sd(x1 = CI_low, x2 = CI_high),
                     CI_low))

# coucal_flight_age_plus_fledge <-
#   data.frame(trait = c("flight_age"),
#              species = c("BC", "BC", "WBC", "WBC"),
#              sex = c("F", "M", "F", "M"),
#              mean = c(mod_flight_age_coefs[1, c(2)] + pull(filter(coucal_flight_age, species == "BC" & sex == "F"), mean),
#                       (mod_flight_age_coefs[1, c(2)] + mod_flight_age_coefs[2, c(2)] + pull(filter(coucal_flight_age, species == "BC" & sex == "M"), mean)),
#                       (mod_flight_age_coefs[1, c(2)] + mod_flight_age_coefs[3, c(2)] + pull(filter(coucal_flight_age, species == "WBC" & sex == "F"), mean)),
#                       (mod_flight_age_coefs[1, c(2)] + mod_flight_age_coefs[3, c(2)] + mod_flight_age_coefs[4, c(2)] + mod_flight_age_coefs[2, c(2)]) + pull(filter(coucal_flight_age, species == "WBC" & sex == "M"), mean)),
#              CI_low = c(mod_flight_age_coefs[1, c(5)] + pull(filter(coucal_flight_age, species == "BC" & sex == "F"), mean),
#                         (mod_flight_age_coefs[1, c(5)] + mod_flight_age_coefs[2, c(5)]) + pull(filter(coucal_flight_age, species == "BC" & sex == "M"), mean),
#                         (mod_flight_age_coefs[1, c(5)] + mod_flight_age_coefs[3, c(5)]) + pull(filter(coucal_flight_age, species == "WBC" & sex == "F"), mean),
#                         (mod_flight_age_coefs[1, c(5)] + mod_flight_age_coefs[3, c(5)] + mod_flight_age_coefs[4, c(5)] + mod_flight_age_coefs[2, c(5)]) + pull(filter(coucal_flight_age, species == "WBC" & sex == "M"), mean)),
#              CI_high = c(mod_flight_age_coefs[1, c(6)] + pull(filter(coucal_flight_age, species == "BC" & sex == "F"), mean),
#                          (mod_flight_age_coefs[1, c(6)] + mod_flight_age_coefs[2, c(6)]) + pull(filter(coucal_flight_age, species == "BC" & sex == "M"), mean),
#                          (mod_flight_age_coefs[1, c(6)] + mod_flight_age_coefs[3, c(6)]) + pull(filter(coucal_flight_age, species == "WBC" & sex == "F"), mean),
#                          (mod_flight_age_coefs[1, c(6)] + mod_flight_age_coefs[3, c(6)] + mod_flight_age_coefs[4, c(6)] + mod_flight_age_coefs[2, c(6)]) + pull(filter(coucal_flight_age, species == "WBC" & sex == "M"), mean)),
#              n_inds = c(filter(flight_dat, species == "BC" & Sex == "Female") %>%
#                           summarise(n_ = n_distinct(Ind_ID)) %>%
#                           pull(n_),
#                         filter(flight_dat, Spp == "BC" & Sex == "Male") %>%
#                           summarise(n_ = n_distinct(Ind_ID)) %>%
#                           pull(n_),
#                         filter(flight_dat, Spp == "WBC" & Sex == "Female") %>%
#                           summarise(n_ = n_distinct(Ind_ID)) %>%
#                           pull(n_),
#                         filter(flight_dat, Spp == "BC" & Sex == "Male") %>%
#                           summarise(n_ = n_distinct(Ind_ID)) %>%
#                           pull(n_)),
#              n_nests = c(filter(flight_dat, Spp == "BC" & Sex == "Female") %>%
#                            summarise(n_ = n_distinct(Nst_No)) %>%
#                            pull(n_),
#                          filter(flight_dat, Spp == "BC" & Sex == "Male") %>%
#                            summarise(n_ = n_distinct(Nst_No)) %>%
#                            pull(n_),
#                          filter(flight_dat, Spp == "WBC" & Sex == "Female") %>%
#                            summarise(n_ = n_distinct(Nst_No)) %>%
#                            pull(n_),
#                          filter(flight_dat, Spp == "WBC" & Sex == "Male") %>%
#                            summarise(n_ = n_distinct(Nst_No)) %>%
#                            pull(n_)),
#              n_years = c(filter(flight_dat, Spp == "BC" & Sex == "Female") %>%
#                            summarise(n_ = n_distinct(Year)) %>%
#                            pull(n_),
#                          filter(flight_dat, Spp == "BC" & Sex == "Male") %>%
#                            summarise(n_ = n_distinct(Year)) %>%
#                            pull(n_),
#                          filter(flight_dat, Spp == "WBC" & Sex == "Female") %>%
#                            summarise(n_ = n_distinct(Year)) %>%
#                            pull(n_),
#                          filter(flight_dat, Spp == "WBC" & Sex == "Male") %>%
#                            summarise(n_ = n_distinct(Year)) %>%
#                            pull(n_))) %>%
#   mutate(sd = ifelse(!is.na(CI_low),
#                      approx_sd(x1 = CI_low, x2 = CI_high),
#                      CI_low))

#### "Age_1st_Flight" as dependent variable ----
# perhaps not as informative as "days_since_fledging" since it doesn't control
# for the age when the chick left the nest
# flight_mod2 <- lmer(Age_1st_Flight ~ Sex * Spp +
#                       (1 | Nst_No) + (1 | Year),
#                     data = flight_dat)
# 
# flight_mod2_coefs <-
#   model_parameters(flight_mod2) %>%
#   as.data.frame(.)
# 
# plot(allEffects(flight_mod2))

rm(flight_dat, mod_flight_age_WBC, mod_flight_age_BC, new_data_BC, new_data_WBC,
   predict_BC, predict_WBC)