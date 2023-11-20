# Survival analysis of fledglings based on body mass
## load (and install if necessary) packages need for project

# a vector of all the packages needed in the project
packages_required_in_project <- c("RMark",
                                  "tidyverse",
                                  "readxl",
                                  "BaSTA",
                                  "pbapply",
                                  "RColorBrewer",
                                  "grid",
                                  "Rmisc",
                                  "gss",
                                  "arm",
                                  "partR2",
                                  "parameters",
                                  "standardize",
                                  "colorBlindness",
                                  "ggthemes",
                                  "patchwork",
                                  "gt",
                                  "rptR",
                                  "tidybayes",
                                  "broom.mixed",
                                  "effects",
                                  "patchwork",
                                  "devtools",
                                  "unmarked",
                                  "R2ucare",
                                  "marked",
                                  "merTools",
                                  "bootpredictlme4",
                                  "extrafont",
                                  "survminer",
                                  "bdsmatrix",
                                  "coxme",
                                  "epitools",
                                  "survival",
                                  "magrittr",
                                  "ggpubr",
                                  "reshape2",
                                  "sp",
                                  "adehabitatLT",
                                  "reshape",
                                  "coefplot",
                                  "mapview",
                                  "lubridate",
                                  "effects",
                                  "rms")

# of the required packages, check if some need to be installed
new.packages <- 
  packages_required_in_project[!(packages_required_in_project %in% 
                                   installed.packages()[,"Package"])]

# install all packages that are not locally available
if(length(new.packages)) install.packages(new.packages)

# load all the packages into the current R session
lapply(packages_required_in_project, require, character.only = TRUE)

# import data
status_dat_fledglings <- 
  # read raw data
  read.csv("data/raw/Coucal_chick_survival_2001-2019_20200129.csv", 
           header = TRUE, stringsAsFactors = FALSE, na.strings = c("", " ", "NA")) %>% 
  
  # rename ring_ID column
  dplyr::rename(ring_ID = Ring_ID) %>% 
  
  # make all entries lower case for consistency
  mutate(Fledged_status = tolower(Fledged.),
         site = tolower(site)) %>% 
  
  # select variables of interest
  dplyr::select(species, ring_ID, lab_no, sex, year, site, nest_ID, pref_age, 
                Fledged_status, entry, exit, event, ageC, lay_date, hatch_order,
                Relative_hatch_order, age_diff, clutch_size, brood_size, No_fledgelings, 
                Fledge_age, Fledge_mass, Fledge_tarsus, Fledge_FS, Fledge_FSO) %>% 
  
  # remove all white space from data
  mutate(across(everything(), ~str_trim(.x))) %>% 
  mutate(across(everything(), ~str_replace_all(.x, fixed(" "), ""))) %>% 
  
  # specify empty data as NA
  mutate(across(everything(), ~gsub("^$|^ $", NA, .x))) %>% 
  
  # exclude all individuals that died in the nest
  filter(Fledged_status == "yes") %>% 
  
  # classify columns
  mutate(sex = as.factor(sex),
         ageC = as.numeric(ageC),
         postf_age = as.numeric(postf_age),
         postf_status = as.numeric(postf_status),
         hatch_order = as.numeric(hatch_order),
         pref_age = as.numeric(pref_age),
         Relative_hatch_order = as.numeric(Relative_hatch_order),
         age_diff = as.numeric(age_diff),
         clutch_size = as.numeric(clutch_size),
         brood_size = as.numeric(brood_size),
         No_fledgelings = as.numeric(No_fledgelings),
         Fledge_age = as.numeric(Fledge_age),
         Fledge_mass = as.numeric(Fledge_mass),
         Fledge_tarsus = as.numeric(Fledge_tarsus),
         Fledge_FS = as.numeric(Fledge_FS),
         Fledge_FSO = as.numeric(Fledge_FSO)) %>% 
  
  # remove rows with missing sex and age, and post-fledged status info
  filter(!is.na(sex) & !is.na(ageC) & !is.na(postf_status)) %>% 
  
  # make a unique id for each individual
  mutate(ind_ID = paste(nest_ID, lab_no, ring_ID, sep = "_"),
         
         # create the age of entry into the data (all at age 15)
         entry = 15,
         
         # specify that the minimum age is 15 for death and censoring
         exit = ifelse(postf_age <= 15, 16, postf_age)) %>% 
  
  # rename the postf_age column as "exit" to specify the the individual died 
  # or was last censored, rename the postf_status as "event" to specify if 
  # the individual died (1) or was censored (0)
  dplyr::rename(event = postf_status) #%>%
  
  # consolidate to variables of interest
  # dplyr::select(species, ind_ID, nest_ID, year, sex, entry, exit, event, 
  #               hatch_order, Relative_hatch_order, age_diff, clutch_size, 
  #               brood_size, No_fledgelings, Fledge_age, Fledge_mass, 
  #               Fledge_tarsus, Fledge_FS, Fledge_FSO)

# Use colnames and tolower to convert all column names to lowercase
colnames(status_dat_fledglings) <- tolower(colnames(status_dat_fledglings))

status_dat_fledglings %>% 
  group_by(sex, event) %>% 
  dplyr::summarise(n_distinct(ind_id))

# are females more likely to hatch earlier?
m1 <- glmer(sex ~ relative_hatch_order + (1 | nest_id), 
               data = status_dat_fledglings %>% filter(species == "BC"), 
               family = binomial)

# is tarsus length a good predictor of body mass for males and females?
m2 <- lmer(fledge_mass ~ sex*fledge_tarsus + (1 | nest_id), 
               data = status_dat_fledglings %>% filter(species == "BC"))

coefplot(m2)

# make a dataframe of the residuals and the unique capture events
mod_BC <- lmer(fledge_mass ~ fledge_tarsus + (1 | nest_id), 
               data = status_dat_fledglings %>% filter(species == "BC")) # 125
mod_WBC <- lmer(fledge_mass ~ fledge_tarsus + (1 | nest_id), 
                data = status_dat_fledglings %>% filter(species == "WBC")) # 106

status_dat_fledglings %>% filter(species == "BC") %>% nrow()
status_dat_fledglings %>% filter(species == "WBC") %>% nrow()

# extract fitted values
mod_BC_fits <- 
  as.data.frame(effect(term = "fledge_tarsus", mod = mod_BC, 
                       xlevels = list(fledge_tarsus = seq(status_dat_fledglings %>% filter(species == "BC") %>% pull(fledge_tarsus) %>% min(na.rm = TRUE), 
                                                          status_dat_fledglings %>% filter(species == "BC") %>% pull(fledge_tarsus) %>% max(na.rm = TRUE), 1)))) %>% 
  mutate(species = "BC")

mod_WBC_fits <- 
  as.data.frame(effect(term = "fledge_tarsus", mod = mod_WBC, 
                       xlevels = list(fledge_tarsus = seq(status_dat_fledglings %>% filter(species == "WBC") %>% pull(fledge_tarsus) %>% min(na.rm = TRUE), 
                                                          status_dat_fledglings %>% filter(species == "WBC") %>% pull(fledge_tarsus) %>% max(na.rm = TRUE), 1)))) %>% 
  mutate(species = "WBC")

mod_fits <- 
  bind_rows(mod_BC_fits, mod_WBC_fits)

ggplot() +
  geom_point(data = status_dat_fledglings,
             aes(y = fledge_mass, x = fledge_tarsus)) +
  geom_line(data = mod_fits, aes(x = fledge_tarsus, y = fit),
            lwd = 0.5) +
  facet_grid(species ~ .)

status_dat_fledglings_res_BC <- 
  bind_cols(status_dat_fledglings %>% filter(species == "BC"),
            data.frame(mod_res = residuals(mod_BC)))

status_dat_fledglings_res_WBC <- 
  bind_cols(status_dat_fledglings %>% filter(species == "WBC"),
            data.frame(mod_res = residuals(mod_WBC)))

status_dat_fledglings_res <- 
  bind_rows(status_dat_fledglings_res_BC,
            status_dat_fledglings_res_WBC) %>% 
  distinct()

cox_BC_res <- 
  coxme(Surv(entry, exit, event) ~ mod_res + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings_res %>% filter(species == "BC"))  # postfledging survival
summary(cox_BC_res)
plot(predict(cox_BC_res, type = "lp"), status_dat_fledglings_res %>% filter(species == "BC") %>% pull(mod_res))
plot(predict(cox_BC_res, type = "risk"), status_dat_fledglings_res %>% filter(species == "BC") %>% pull(mod_res))
ggplot(Predict(cox_BC_res, mod_res))


cox_BC_res <- 
  cph(Surv(entry, exit, event) ~ mod_res, data =  status_dat_fledglings_res %>% filter(species == "BC")) 
survplot(cox_BC_res, levels.only = TRUE)
ggplot(Predict(cox_BC_res, mod_res))


plot(survfit(cox_BC_res))

cox_WBC_res <- 
  coxme(Surv(entry, exit, event) ~ mod_res*sex + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings_res %>% filter(species == "WBC"))  # postfledging survival
summary(cox_WBC_res)
predict(cox_BC_res)

status_dat_fledglings %>% 
  group_by(lab_no) %>% 
  dplyr::summarise(n = n()) %>% 
  arrange(desc(n))

cap_05_09_std_pca_ten_res

mod_res <-
  lmer(mod_res ~ date_deviance + first_date + last_date +
         (1 | capture_id) + (1 | year_),
       data = cap_05_09_std_pca_ten_res)

tidy_mod_weight <-
  tidy(mod_res, conf.int = TRUE, conf.method = "boot", nsim = 1000)

ggplot() +
  geom_density(data = status_dat_fledglings, aes(x = mod_res, fill = sex), alpha = 0.4)

ggplot() +
  geom_density(data = status_dat_fledglings, aes(x = fledge_mass, fill = sex), alpha = 0.4)

ggplot() +
  geom_density(data = status_dat_fledglings, aes(x = fledge_tarsus, fill = sex), alpha = 0.4)

ggplot() +
  geom_density(data = status_dat_fledglings, aes(x = hatch_order, fill = sex), alpha = 0.4)

cox_BC_sex <- 
  coxme(Surv(entry, exit, event) ~ sex + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings %>% filter(species == "BC"))  # postfledging survival
summary(cox_BC_sex)

stan.fit <- sshzd(Surv(futime,status)~futime+age,data=stan)
## Evaluate fitted hazard
hzdrate.sshzd(stan.fit,data.frame(futime=c(10,20),age=c(20,30)))
## Plot lambda(t,age=20)
tt <- seq(0,60,leng=101)
hh <- hzdcurve.sshzd(stan.fit,tt,data.frame(age=50))
plot(tt,hh,type="l")

cox_BC_hatch_order <- 
  coxme(Surv(entry, exit, event) ~ relative_hatch_order + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings_res %>% filter(species == "BC"))  # postfledging survival
summary(cox_BC_hatch_order)

cox_BC_hatch_order <- 
  coxme(Surv(entry, exit, event) ~ relative_hatch_order + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings %>% filter(species == "BC"))  # postfledging survival
summary(cox_BC_hatch_order)

cox_BC_fledge_mass <- 
  coxme(Surv(entry, exit, event) ~ fledge_mass + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings %>% filter(species == "BC"))  # postfledging survival
summary(cox_BC_fledge_mass)

cox_BC_fledge_tarsus <- 
  coxme(Surv(entry, exit, event) ~ fledge_tarsus + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings %>% filter(species == "BC"))  # postfledging survival
summary(cox_BC_fledge_tarsus)

cox_BC_sex_hatch_order <- 
  coxme(Surv(entry, exit, event) ~ sex + relative_hatch_order + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings %>% filter(species == "BC"))  # postfledging survival
summary(cox_BC_sex_hatch_order)

cox_BC_sexXhatch_order <- 
  coxme(Surv(entry, exit, event) ~ sex * relative_hatch_order + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings %>% filter(species == "BC"))  # postfledging survival
summary(cox_BC_sexXhatch_order)

cox_BC_sex_fledge_mass <- 
  coxme(Surv(entry, exit, event) ~ sex + fledge_mass + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings %>% filter(species == "BC"))  # postfledging survival
summary(cox_BC_sex_fledge_mass)

cox_BC_sex_fledge_tarsus <- 
  coxme(Surv(entry, exit, event) ~ sex + fledge_tarsus + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings %>% filter(species == "BC"))  # postfledging survival
summary(cox_BC_sex_fledge_tarsus)

cox_BC_sexXfledge_mass <- 
  coxme(Surv(entry, exit, event) ~ sex * fledge_mass + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings %>% filter(species == "BC"))  # postfledging survival
summary(cox_BC_sexXfledge_mass)

cox_BC_sexXfledge_tarsus <- 
  coxme(Surv(entry, exit, event) ~ sex * fledge_tarsus + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings %>% filter(species == "BC"))  # postfledging survival
summary(cox_BC_sexXfledge_tarsus)

cox_WBC_sex_hatch_order <- 
  coxme(Surv(entry, exit, event) ~ sex + relative_hatch_order + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings %>% filter(species == "WBC"))  # postfledging survival
summary(cox_WBC_sex_hatch_order)

cox_WBC_sexXhatch_order <- 
  coxme(Surv(entry, exit, event) ~ sex * relative_hatch_order + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings %>% filter(species == "WBC"))  # postfledging survival
summary(cox_WBC_sexXhatch_order)

cox_WBC_sex_fledge_mass <- 
  coxme(Surv(entry, exit, event) ~ sex + fledge_mass + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings %>% filter(species == "WBC"))  # postfledging survival
summary(cox_WBC_sex_fledge_mass)

cox_WBC_sex_fledge_tarsus <- 
  coxme(Surv(entry, exit, event) ~ sex + fledge_tarsus + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings %>% filter(species == "WBC"))  # postfledging survival
summary(cox_WBC_sex_fledge_tarsus)

cox_WBC_sexXfledge_mass <- 
  coxme(Surv(entry, exit, event) ~ sex * fledge_mass + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings %>% filter(species == "WBC"))  # postfledging survival
summary(cox_WBC_sexXfledge_mass)

cox_WBC_sexXfledge_tarsus <- 
  coxme(Surv(entry, exit, event) ~ sex * fledge_tarsus + 
          (1|year) + (1|nest_id), data =  status_dat_fledglings %>% filter(species == "WBC"))  # postfledging survival
summary(cox_WBC_sexXfledge_tarsus)