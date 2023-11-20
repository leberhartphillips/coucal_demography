# load packages
source("R/project/project_libraries.R")

# load functions
function.sources = list.files(path = "R/functions",
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
sapply(function.sources, source)

# load capture histories
data.sources = list.files(path = "cooked_data", 
                          pattern="*ch.rds$", full.names = TRUE, 
                          ignore.case = TRUE)
sapply(data.sources, load, .GlobalEnv)

#### read the data into R ----
# pre & postfleding survival combined (variable postf_age, postf_status) 
dat = read.csv("data/raw/Coucal_chick_survival_2001-2019_20200129.csv", header=T)
head(dat)
str(dat)

status_dat_nestlings <- 
  # read raw data
  read.csv("data/raw/Coucal_chick_survival_2001-2019_20200129.csv", 
           header = TRUE, stringsAsFactors = FALSE, na.strings = c("", " ", "NA")) %>% 
  
  # rename ring_ID column
  dplyr::rename(ring_ID = Ring_ID) %>% 
  
  # make all entries lower case for consistency
  mutate(Fledged_status = tolower(Fledged.),
         site = tolower(site)) %>% 
  
  # select variables of interest
  select(species, ring_ID, lab_no, sex, year, site, nest_ID, pref_age, 
         Fledged_status, postf_age, postf_status, ageC, lay_date, hatch_order) %>% 
  
  # remove all white space from data
  mutate(across(everything(), ~str_trim(.x))) %>% 
  mutate(across(everything(), ~str_replace_all(.x, fixed(" "), ""))) %>% 
  
  # specify empty data as NA
  mutate(across(everything(), ~gsub("^$|^ $", NA, .x))) %>% 
  
  # classify columns
  mutate(sex = as.factor(sex),
         ageC = as.numeric(ageC),
         postf_age = as.numeric(postf_age),
         postf_status = as.numeric(postf_status),
         hatch_order = as.numeric(hatch_order),
         pref_age = as.numeric(pref_age)) %>% 
  
  # remove rows with missing sex, age, and fledged status info
  filter(!is.na(sex) & !is.na(ageC) & !is.na(Fledged_status)) %>% 
  
  # make a unique id for each individual
  mutate(ind_ID = paste(nest_ID, lab_no, ring_ID, sep = "_"),
         
  # create the age of entry into the data (all at age 0)
         entry = 0,
  
  # specifiy the event (0 = alive, 1 = found dead)
         event = ifelse(Fledged_status == "yes", 0, 1)) %>% 
  
  # rename the pref_age column as "exit" to specify the the individual died 
  # or was last censored
  dplyr::rename(exit = pref_age) %>% 
  
  # consolidate to variables of interest
  dplyr::select(species, ind_ID, nest_ID, year, sex, entry, exit, event, hatch_order)

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
                Fledged_status, postf_age, postf_status, ageC, lay_date, hatch_order) %>% 
  
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
         pref_age = as.numeric(pref_age)) %>% 
  
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
  dplyr::rename(event = postf_status) %>% 
  
  # consolidate to variables of interest
  dplyr::select(species, ind_ID, nest_ID, year, sex, entry, exit, event, hatch_order)

status_dat_all <- 
  # read raw data
  read.csv("data/raw/Coucal_chick_survival_2001-2019_20200129.csv", 
           header = TRUE, stringsAsFactors = FALSE, na.strings = c("", " ", "NA")) %>% 
  
  # rename ring_ID column
  dplyr::rename(ring_ID = Ring_ID) %>% 
  
  # make all entries lower case for consistency
  mutate(Fledged_status = tolower(Fledged.),
         site = tolower(site)) %>% 
  
  # select variables of interest
  select(species, ring_ID, lab_no, sex, year, site, nest_ID, pref_age, 
         Fledged_status, postf_age, postf_status, ageC, statusC, 
         lay_date, hatch_order) %>% 
  
  # remove all white space from data
  mutate(across(everything(), ~str_trim(.x))) %>% 
  mutate(across(everything(), ~str_replace_all(.x, fixed(" "), ""))) %>% 
  
  # specify empty data as NA
  mutate(across(everything(), ~gsub("^$|^ $", NA, .x))) %>% 
  
  # exclude all individuals that died in the nest
  # filter(Fledged_status == "yes") %>% 
  
  # classify columns
  mutate(sex = as.factor(sex),
         ageC = as.numeric(ageC),
         postf_age = as.numeric(postf_age),
         postf_status = as.numeric(postf_status),
         hatch_order = as.numeric(hatch_order),
         pref_age = as.numeric(pref_age)) %>% 
  
  # remove rows with missing sex, age, and status info
  filter(!is.na(sex) & !is.na(ageC) & !is.na(statusC)) %>% 
  
  # make a unique id for each individual
  mutate(ind_ID = paste(nest_ID, lab_no, ring_ID, sep = "_"),
         
         # create the age of entry into the data (all at age 15)
         entry = 0,
         
         # specify the age of death or censoring
         exit = ageC,
         
         # make the event numeric and specify if 
         # the individual died (1) or was censored (0)
         event = as.numeric(statusC)) %>% 
  
  # consolidate to variables of interest
  dplyr::select(species, ind_ID, nest_ID, year, sex, entry, exit, event, hatch_order)

BC_nest_dat <- 
  status_dat_nestlings %>% 
  filter(species == "BC")

BC_fled_dat <- 
  status_dat_fledglings %>% 
  filter(species == "BC")

BC_all_dat <- 
  status_dat_all %>% 
  filter(species == "BC")

WBC_nest_dat <- 
  status_dat_nestlings %>% 
  filter(species == "WBC")

WBC_fled_dat <- 
  status_dat_fledglings %>% 
  filter(species == "WBC")

WBC_all_dat <- 
  status_dat_all %>% 
  filter(species == "WBC")

#### Hazard rate bootstrap ----

# set attempt to 0 at start of each loop
attempt <- 0

tt.a <- seq(0, 60, 1)
alpha_value <- 1.4

# store simulated estimates only if peak >= 1 and <= 10 and it's less than
# 100 attempts
while( attempt <= 100 ) {
  
  # next attempt
  attempt <- attempt + 1
  
  boot_BC_all_dat <- 
    WBC_all_dat %>% 
    dplyr::group_by(nest_ID) %>%
    dplyr::sample_n(1)
  
  try(
    BC.hz.M.a <- sshzd(Surv(exit, event, entry) ~ exit, 
                       data = filter(boot_BC_all_dat, sex == "M"), 
                       alpha = alpha_value)
  )
  # simulate an estimate
  try(
    est.M.a_boot <- hzdcurve.sshzd(object = BC.hz.M.a, time = tt.a, se = TRUE)
  )
}
# store fixed effects
mod6_sim@fixef[i, ] <- mod6_sim_try@fixef


boot_BC_all_dat <- 
  BC_all_dat %>% 
  dplyr::group_by(nest_ID) %>%
  dplyr::sample_n(1)



BC.hz.M.a <- sshzd(Surv(exit, event, entry) ~ exit, 
                   data = filter(boot_BC_all_dat, sex == "M"), 
                   alpha = alpha_value)

project(object = BC.hz.M.a)

est.M.a_boot <- hzdcurve.sshzd(object = BC.hz.M.a, time = tt.a, se = TRUE)

M_WBC_all_hazard_function_est_boot <- 
  data.frame(species = rep("BC", length(tt.a)), 
             sex = rep("M", length(tt.a)), 
             stage = rep("all", length(tt.a)),
             age = tt.a) %>% 
  mutate(estimate = 1 - est.M.a_boot$fit,
         upper = 1 - est.M.a_boot$fit * exp(1.96 * est.M.a_boot$se),
         lower = 1 - est.M.a_boot$fit / exp(1.96 * est.M.a_boot$se))

expand.grid(species = species, 
            age = time_vector,
            sex = c("M", "F")) %>% 
  mutate(estimate = c(est.M.a_boot$fit, est.M.a_boot$fit))

immigrant_pop_size = 100
k = 4
HSR = 0.4955
h = 1/2.9
egg_survival = 0.32
ISR = 0.738
immigrant_pop_size = 100
fledge_age = 15
flight_age = 36
first_year = 2001
bootstrap_name = "BC_one"
num_boot = 1
species = "BC"
iter_add = 1

# transform the daily nestling survival (DCS) to apparent fledgling success
# by calculating the product of all DCS estimates:
coucal_nestling_survival <-
  M_WBC_all_hazard_function_est_boot %>% 
  filter(age < fledge_age) %>% 
  group_by(sex) %>% 
  dplyr::summarise(value = prod(estimate), .groups = 'drop') %>%
  mutate(stage = "nestling",
         rate = "survival")

coucal_groundling_survival <-
  M_WBC_all_hazard_function_est_boot %>% 
  filter(age < flight_age) %>% 
  group_by(sex) %>% 
  dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
  mutate(stage = "groundling",
         rate = "survival")

coucal_fledgling_survival <-
  M_WBC_all_hazard_function_est_boot %>% 
  filter(age >= flight_age) %>% 
  group_by(sex) %>% 
  dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
  mutate(stage = "fledgling",
         rate = "survival")

adult_F_immigrants <- immigrant_pop_size * (1 - ISR)
adult_M_immigrants <- immigrant_pop_size * ISR

coucal_adult_immigration <- 
  data.frame(sex = c("Female", "Male"),
             value = c(adult_F_immigrants, adult_M_immigrants),
             stage = c("adult"),
             rate = c("immigration"))

coucal_egg_survival <- 
  data.frame(sex = NA,
             value = egg_survival,
             stage = c("egg"),
             rate = c("survival"))

ggplot(M_WBC_all_hazard_function_est_boot) +
  geom_line(aes(y = estimate, x = age, color = sex)) +
  geom_ribbon(aes(x = age, ymax = upper, ymin = lower, fill = sex), alpha = 0.25) +
  # geom_vline(xintercept = 15, linetype = "dashed", alpha = 0.5, color = "green") +
  # geom_vline(xintercept = 15, linetype = "dashed", alpha = 0.5, color = "green") +
  geom_vline(xintercept = 36, linetype = "dashed", alpha = 0.5, color = "grey70") +
  geom_vline(xintercept = 15, linetype = "dashed", alpha = 0.5, color = "grey70") +
  # facet_grid(species ~ .) +
  luke_theme +
  scale_fill_manual(values = plot_palette_sex) +
  scale_color_manual(values = plot_palette_sex) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0.8, 1))


#### Estimate Hazard Function using a smoother spline from the gss package
# alpha 1.4 is the default value (less may be over fitting)
alpha_value <- 1.4
BC.hz.M.a <- sshzd(Surv(exit, event, entry)~exit, 
                   data = filter(BC_all_dat, sex == "M"), 
                   alpha = alpha_value) 
BC.hz.F.a <- sshzd(Surv(exit, event, entry)~exit, 
                   data = filter(BC_all_dat, sex == "F"), 
                   alpha = alpha_value)
WBC.hz.M.a <- sshzd(Surv(exit, event, entry)~exit, 
                    data = filter(WBC_all_dat, sex == "M"), 
                    alpha = alpha_value)
WBC.hz.F.a <- sshzd(Surv(exit, event, entry)~exit, 
                    data = filter(WBC_all_dat, sex == "F"), 
                    alpha = alpha_value)

# predicted hazard
# Female all
est.F.a <- hzdrate.sshzd(BC.hz.F.a, tt.a, se = TRUE)
F_BC_all_hazard_function_est <- 
  data.frame(species = rep("BC", length(tt.a)), 
             sex = rep("F", length(tt.a)), 
             stage = rep("all", length(tt.a)),
             age = tt.a)
F_BC_all_hazard_function_est$estimate <- 1 - est.F.a$fit
F_BC_all_hazard_function_est$upper <- 1 - est.F.a$fit*exp(1.96*est.F.a$se)
F_BC_all_hazard_function_est$lower <- 1 - est.F.a$fit/exp(1.96*est.F.a$se)

# Male all
est.M.a <- hzdrate.sshzd(BC.hz.M.a, tt.a, se = TRUE)
M_BC_all_hazard_function_est <- 
  data.frame(species = rep("BC", length(tt.a)), 
             sex = rep("M", length(tt.a)), 
             stage = rep("all", length(tt.a)),
             age = tt.a)
M_BC_all_hazard_function_est$estimate <- 1 - est.M.a$fit
M_BC_all_hazard_function_est$upper <- 1 - est.M.a$fit*exp(1.96*est.M.a$se)
M_BC_all_hazard_function_est$lower <- 1 - est.M.a$fit/exp(1.96*est.M.a$se)

# predicted hazard
# Female all
est.F.a <- hzdrate.sshzd(WBC.hz.F.a, tt.a, se = TRUE)
F_WBC_all_hazard_function_est <- 
  data.frame(species = rep("WBC", length(tt.a)), 
             sex = rep("F", length(tt.a)), 
             stage = rep("all", length(tt.a)),
             age = tt.a) %>% 
  mutate(estimate = 1 - est.F.a$fit,
         upper = 1 - est.F.a$fit*exp(1.96*est.F.a$se),
         lower = 1 - est.F.a$fit/exp(1.96*est.F.a$se))

# Male all
est.M.a <- hzdrate.sshzd(WBC.hz.M.a, tt.a, se = TRUE)
M_WBC_all_hazard_function_est <- 
  data.frame(species = rep("WBC", length(tt.a)), 
             sex = rep("M", length(tt.a)), 
             stage = rep("all", length(tt.a)),
             age = tt.a)
M_WBC_all_hazard_function_est$estimate <- 1 - est.M.a$fit
M_WBC_all_hazard_function_est$upper <- 1 - est.M.a$fit*exp(1.96*est.M.a$se)
M_WBC_all_hazard_function_est$lower <- 1 - est.M.a$fit/exp(1.96*est.M.a$se)

all_hazard_function_est <- 
  bind_rows(F_WBC_all_hazard_function_est,
            M_WBC_all_hazard_function_est,
            F_BC_all_hazard_function_est,
            M_BC_all_hazard_function_est)

ggplot(all_hazard_function_est) +
  geom_line(aes(y = estimate, x = age, color = sex)) +
  geom_ribbon(aes(x = age, ymax = upper, ymin = lower, fill = sex), alpha = 0.25) +
  # geom_vline(xintercept = 15, linetype = "dashed", alpha = 0.5, color = "green") +
  # geom_vline(xintercept = 15, linetype = "dashed", alpha = 0.5, color = "green") +
  geom_vline(xintercept = 36, linetype = "dashed", alpha = 0.5, color = "grey70") +
  geom_vline(xintercept = 15, linetype = "dashed", alpha = 0.5, color = "grey70") +
  facet_grid(species ~ .) +
  luke_theme +
  scale_fill_manual(values = plot_palette_sex) +
  scale_color_manual(values = plot_palette_sex) +
  theme(legend.position = "none")

# time is exit
# status is event
# age is entry

fit <- gssanova(cbind(exit + 0.01, event) ~ entry, data = filter(BC_all_dat, sex == "M"),
                family = "weibull")

test <- sshzd(Surv(exit, event) ~ sqrt(exit) * entry, data = filter(BC_all_dat, sex == "M"))

project(BC.hz.M.a, inc = c("exit", "entry"))

sshzd(Surv(futime, status) ~ futime * age, data = stan)

#### BC hazard function smoothing ----
# hazard functions by treatment period
# note Surv() in the sshzd() function takes the arguments in a different order than Surv() in the survival package
# hazard function = hz
BC.hz.M.f <- sshzd(Surv(exit, event, entry)~exit, data = filter(BC_fled_dat, sex == "M"), alpha = 1.4) #alpha=0.25 will overfit
BC.hz.F.f <- sshzd(Surv(exit, event, entry)~exit, data = filter(BC_fled_dat, sex == "F"), alpha = 1.4)

BC.hz.M.n <- sshzd(Surv(exit, event, entry)~exit, data = filter(BC_nest_dat, sex == "M"), alpha = 1.4) #alpha=0.25 will overfit
BC.hz.F.n <- sshzd(Surv(exit, event, entry)~exit, data = filter(BC_nest_dat, sex == "F"), alpha = 1.4)

BC.hz.M.a <- sshzd(Surv(exit, event, entry)~exit, data = filter(BC_all_dat, sex == "M"), alpha = 1.4) #alpha=0.25 will overfit
BC.hz.F.a <- sshzd(Surv(exit, event, entry)~exit, data = filter(BC_all_dat, sex == "F"), alpha = 1.4)
# assess the plausibility of a proportional hazard model, or an additive model in log hazard,
# check the Kullback-Leibler projection (want the first 2 elements to be very small and the "check" to be 1) 
project(BC.hz.M.a, inc = c("exit", "entry"))
project(BC.hz.F.a, inc = c("exit", "entry"))

tt.a <- seq(0, 70, 1)
new <- data.frame(entry = c(0, 0), exit = c(15, 35.5))
est <- hzdrate.sshzd(BC.hz.M.a, new, se = TRUE)
age1 <- c(15, 35.5)
sex2 <- c(1,0)
test_df <- expand.grid(exit = seq(0, 70, 1), sex2 = c(1, 0))

BC_all_dat$sex2 <- ifelse(BC_all_dat$sex == "M", 1, 0)

# Error in `[.data.frame`(x, , object$terms[[label]]$vlist) : 
#   undefined columns selected

test <- sshzd(Surv(exit, event, entry)~exit*sex2, data = BC_all_dat, alpha = 1.4)
rates_test <- hzdcurve.sshzd(test, time = tt.a, covariates = test_df)

curve_test3 <- hzdcurve.sshzd(weibull_test, time = tt.a, covariates = sex2)
curve_test <- 
  hzdcurve.sshzd(object = weibull_test, time = tt.a)

curve_test <- 
  hzdcurve.sshzd(object = BC.hz.M.a, time = tt.a)

curve_test1 <- 
  hzdcurve.sshzd(object = BC.hz.M.a, time = tt.a, data.frame(age = age1))

curve_test2 <- 
  hzdcurve.sshzd(object = BC.hz.M.a, time = tt.a, data.frame(age = age2))

curve_test3 <- 
  hzdrate.sshzd(object = BC.hz.M.a, tt.a)

surv_test <- survexp.sshzd(BC.hz.M.a, tt.a, data.frame(age = age0))

plot(tt.a, curve_test2, type="l", lty=1, lwd=2, col="blue", xaxs="i", yaxs="i", 
     # xlim=c(0, 70), ylim=c(0.7, 1), 
     xlab=c("days since hatching"), ylab=c("hazard function"))
lines(tt.a, curve_test3, lty=2, lwd=2, col="blue")
lines(tt.a, curve_test2, lty=2, lwd=2, col="blue")


# calculate the hazard rate function for plotting for a sequence of time
# hazard = hz
tt.n <- seq(0, 15, 1)
tt.f <- seq(16, 70, 1)
tt.a <- seq(0, 70, 1)

# predicted hazard
# Female nestlings
est.M.n <- hzdrate.sshzd(BC.hz.M.n, tt.n, se = TRUE)
BC.hz.M.n.fit = 1 - est.M.n$fit
BC.hz.M.n.hi <- 1 - est.M.n$fit*exp(1.96*est.M.n$se)
BC.hz.M.n.lo <- 1 - est.M.n$fit/exp(1.96*est.M.n$se)

# Male nestlings
est.F.n <- hzdrate.sshzd(BC.hz.F.n, tt.n, se = TRUE)
BC.hz.F.n.fit = 1 - est.F.n$fit
BC.hz.F.n.hi <- 1 - est.F.n$fit*exp(1.96*est.F.n$se)
BC.hz.F.n.lo <- 1 - est.F.n$fit/exp(1.96*est.F.n$se)

# predicted hazard
# Female fledglings
est.M.f <- hzdrate.sshzd(BC.hz.M.f, tt.f, se = TRUE)
BC.hz.M.f.fit = 1 - est.M.f$fit
BC.hz.M.f.hi <- 1 - est.M.f$fit*exp(1.96*est.M.f$se)
BC.hz.M.f.lo <- 1 - est.M.f$fit/exp(1.96*est.M.f$se)

# Male fledglings
est.F.f <- hzdrate.sshzd(BC.hz.F.f, tt.f, se = TRUE)
BC.hz.F.f.fit = 1 - est.F.f$fit
BC.hz.F.f.hi <- 1 - est.F.f$fit*exp(1.96*est.F.f$se)
BC.hz.F.f.lo <- 1 - est.F.f$fit/exp(1.96*est.F.f$se)

# predicted hazard
# Female all
est.F.a <- hzdrate.sshzd(BC.hz.F.a, tt.a, se = TRUE)
F_BC_all_hazard_function_est <- 
  data.frame(species = rep("BC", length(tt.a)), 
             sex = rep("F", length(tt.a)), 
             stage = rep("all", length(tt.a)),
             age = tt.a)
F_BC_all_hazard_function_est$estimate <- 1 - est.F.a$fit
F_BC_all_hazard_function_est$upper <- 1 - est.F.a$fit*exp(1.96*est.F.a$se)
F_BC_all_hazard_function_est$lower <- 1 - est.F.a$fit/exp(1.96*est.F.a$se)

# Male all
est.M.a <- hzdrate.sshzd(BC.hz.M.a, tt.a, se = TRUE)
M_BC_all_hazard_function_est <- 
  data.frame(species = rep("BC", length(tt.a)), 
             sex = rep("M", length(tt.a)), 
             stage = rep("all", length(tt.a)),
             age = tt.a)
M_BC_all_hazard_function_est$estimate <- 1 - est.M.a$fit
M_BC_all_hazard_function_est$upper <- 1 - est.M.a$fit*exp(1.96*est.M.a$se)
M_BC_all_hazard_function_est$lower <- 1 - est.M.a$fit/exp(1.96*est.M.a$se)

# graph for hazard functions in males and females periods
plot(tt.n, BC.hz.M.n.fit, type="l", lty=1, lwd=2, col="blue", xaxs="i", yaxs="i", 
     xlim=c(0, 70), ylim=c(0.7, 1), 
     xlab=c("days since hatching"), ylab=c("hazard function"))
lines(tt.n, BC.hz.M.n.hi, lty=2, lwd=2, col="blue")
lines(tt.n, BC.hz.M.n.lo, lty=2, lwd=2, col="blue")
lines(tt.n, BC.hz.F.n.fit, lty=1, lwd=2, col="red")
lines(tt.n, BC.hz.F.n.hi, lty=2, lwd=2, col="red")
lines(tt.n, BC.hz.F.n.lo, lty=2, lwd=2, col="red")

lines(tt.f, BC.hz.M.f.fit, lty=1, lwd=2, col="blue")
lines(tt.f, BC.hz.M.f.hi, lty=2, lwd=2, col="blue")
lines(tt.f, BC.hz.M.f.lo, lty=2, lwd=2, col="blue")
lines(tt.f, BC.hz.F.f.fit, lty=1, lwd=2, col="red")
lines(tt.f, BC.hz.F.f.hi, lty=2, lwd=2, col="red")
lines(tt.f, BC.hz.F.f.lo, lty=2, lwd=2, col="red")

plot(tt.a, BC.hz.M.a.fit, type="l", lty=1, lwd=2, col="blue", xaxs="i", yaxs="i", 
     xlim=c(0, 70), ylim=c(0.7, 1),
     xlab=c("days since hatching"), ylab=c("hazard function"))
lines(tt.a, BC.hz.M.a.hi, lty=2, lwd=2, col="blue")
lines(tt.a, BC.hz.M.a.lo, lty=2, lwd=2, col="blue")
lines(tt.a, BC.hz.F.a.fit, lty=1, lwd=2, col="red")
lines(tt.a, BC.hz.F.a.hi, lty=2, lwd=2, col="red")
lines(tt.a, BC.hz.F.a.lo, lty=2, lwd=2, col="red")

#### WBC hazard function smoothing ----
# note Surv() in the sshzd() function takes the arguments in a different order than Surv() in the survival package
# hazard function = hz
WBC.hz.M.f <- sshzd(Surv(exit, event, entry)~exit, data = filter(WBC_fled_dat, sex == "M"), alpha = 1.4) #alpha=0.25 will overfit
WBC.hz.F.f <- sshzd(Surv(exit, event, entry)~exit, data = filter(WBC_fled_dat, sex == "F"), alpha = 1.4)

WBC.hz.M.n <- sshzd(Surv(exit, event, entry)~exit, data = filter(WBC_nest_dat, sex == "M"), alpha = 1.4) #alpha=0.25 will overfit
WBC.hz.F.n <- sshzd(Surv(exit, event, entry)~exit, data = filter(WBC_nest_dat, sex == "F"), alpha = 1.4)

WBC.hz.M.a <- sshzd(Surv(exit, event, entry)~exit, data = filter(WBC_all_dat, sex == "M"), alpha = 1.2) #alpha=0.25 will overfit
WBC.hz.F.a <- sshzd(Surv(exit, event, entry)~exit, data = filter(WBC_all_dat, sex == "F"), alpha = 1.2)

# calculate the hazard rate function for plotting for a sequence of time (tt = 1 to 52 weeks by 0.5)
# hazard = hz
tt.n <- seq(0, 15, 1)
tt.f <- seq(16, 70, 1)
tt.a <- seq(0, 70, 1)

# predicted hazard
# Female nestlings
est.M.n <- hzdrate.sshzd(WBC.hz.M.n, tt.n, se = TRUE)
WBC.hz.M.n.fit = 1 - est.M.n$fit
WBC.hz.M.n.hi <- 1 - est.M.n$fit*exp(1.96*est.M.n$se)
WBC.hz.M.n.lo <- 1 - est.M.n$fit/exp(1.96*est.M.n$se)

# Male nestlings
est.F.n <- hzdrate.sshzd(WBC.hz.F.n, tt.n, se = TRUE)
WBC.hz.F.n.fit = 1 - est.F.n$fit
WBC.hz.F.n.hi <- 1 - est.F.n$fit*exp(1.96*est.F.n$se)
WBC.hz.F.n.lo <- 1 - est.F.n$fit/exp(1.96*est.F.n$se)

# predicted hazard
# Female fledglings
est.M.f <- hzdrate.sshzd(WBC.hz.M.f, tt.f, se = TRUE)
WBC.hz.M.f.fit = 1 - est.M.f$fit
WBC.hz.M.f.hi <- 1 - est.M.f$fit*exp(1.96*est.M.f$se)
WBC.hz.M.f.lo <- 1 - est.M.f$fit/exp(1.96*est.M.f$se)

# Male fledglings
est.F.f <- hzdrate.sshzd(WBC.hz.F.f, tt.f, se = TRUE)
WBC.hz.F.f.fit = 1 - est.F.f$fit
WBC.hz.F.f.hi <- 1 - est.F.f$fit*exp(1.96*est.F.f$se)
WBC.hz.F.f.lo <- 1 - est.F.f$fit/exp(1.96*est.F.f$se)

# predicted hazard
# Female all
est.F.a <- hzdrate.sshzd(WBC.hz.F.a, tt.a, se = TRUE)
F_WBC_all_hazard_function_est <- 
  data.frame(species = rep("WBC", length(tt.a)), 
             sex = rep("F", length(tt.a)), 
             stage = rep("all", length(tt.a)),
             age = tt.a)
F_WBC_all_hazard_function_est$estimate <- 1 - est.F.a$fit
F_WBC_all_hazard_function_est$upper <- 1 - est.F.a$fit*exp(1.96*est.F.a$se)
F_WBC_all_hazard_function_est$lower <- 1 - est.F.a$fit/exp(1.96*est.F.a$se)

# Male all
est.M.a <- hzdrate.sshzd(WBC.hz.M.a, tt.a, se = TRUE)
M_WBC_all_hazard_function_est <- 
  data.frame(species = rep("WBC", length(tt.a)), 
             sex = rep("M", length(tt.a)), 
             stage = rep("all", length(tt.a)),
             age = tt.a)
M_WBC_all_hazard_function_est$estimate <- 1 - est.M.a$fit
M_WBC_all_hazard_function_est$upper <- 1 - est.M.a$fit*exp(1.96*est.M.a$se)
M_WBC_all_hazard_function_est$lower <- 1 - est.M.a$fit/exp(1.96*est.M.a$se)

all_hazard_function_est <- 
  bind_rows(F_WBC_all_hazard_function_est,
            M_WBC_all_hazard_function_est,
            F_BC_all_hazard_function_est,
            M_BC_all_hazard_function_est)

ggplot(all_hazard_function_est) +
  geom_line(aes(y = estimate, x = age, color = sex)) +
  geom_ribbon(aes(x = age, ymax = upper, ymin = lower, fill = sex), alpha = 0.25) +
  # geom_vline(xintercept = 15, linetype = "dashed", alpha = 0.5, color = "green") +
  # geom_vline(xintercept = 15, linetype = "dashed", alpha = 0.5, color = "green") +
  geom_vline(xintercept = 36.85, linetype = "dashed", alpha = 0.5, color = "orange") +
  annotate("rect", xmin = 34.71, xmax = 38.99, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "orange") +
  annotate("rect", xmin = 32.73, xmax = 36.94, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "green") +
  geom_vline(xintercept = 34.85, linetype = "dashed", alpha = 0.5, color = "green") +
  facet_grid(species ~ .) +
  luke_theme +
  scale_fill_manual(values = plot_palette_sex) +
  scale_color_manual(values = plot_palette_sex) +
  theme(legend.position = "none")

# graph for hazard functions in males and females periods
plot(tt.n, WBC.hz.M.n.fit, type="l", lty=1, lwd=2, col="blue", xaxs="i", yaxs="i", 
     xlim=c(0, 70), ylim=c(0.7, 1), 
     xlab=c("days since hatching"), ylab=c("hazard function"))
lines(tt.n, WBC.hz.M.n.hi, lty=2, lwd=2, col="blue")
lines(tt.n, WBC.hz.M.n.lo, lty=2, lwd=2, col="blue")
lines(tt.n, WBC.hz.F.n.fit, lty=1, lwd=2, col="red")
lines(tt.n, WBC.hz.F.n.hi, lty=2, lwd=2, col="red")
lines(tt.n, WBC.hz.F.n.lo, lty=2, lwd=2, col="red")

lines(tt.f, WBC.hz.M.f.fit, lty=1, lwd=2, col="blue")
lines(tt.f, WBC.hz.M.f.hi, lty=2, lwd=2, col="blue")
lines(tt.f, WBC.hz.M.f.lo, lty=2, lwd=2, col="blue")
lines(tt.f, WBC.hz.F.f.fit, lty=1, lwd=2, col="red")
lines(tt.f, WBC.hz.F.f.hi, lty=2, lwd=2, col="red")
lines(tt.f, WBC.hz.F.f.lo, lty=2, lwd=2, col="red")

plot(tt.a, WBC.hz.M.a.fit, type="l", lty=1, lwd=2, col="blue", xaxs="i", yaxs="i", 
     xlim=c(0, 70), ylim=c(0.7, 1),
     xlab=c("days since hatching"), ylab=c("hazard function"))
lines(tt.a, WBC.hz.M.a.hi, lty=2, lwd=2, col="blue")
lines(tt.a, WBC.hz.M.a.lo, lty=2, lwd=2, col="blue")
lines(tt.a, WBC.hz.F.a.fit, lty=1, lwd=2, col="red")
lines(tt.a, WBC.hz.F.a.hi, lty=2, lwd=2, col="red")
lines(tt.a, WBC.hz.F.a.lo, lty=2, lwd=2, col="red")

legend("topright", legend = c("Female", "Male"), cex=2, lty=c(1,1), col=c("red", "blue"), lwd=2.5, bty="n")
axis(1, labels=F, tick=T)
# axis(2, at=seq(0,0.1,0.01), labels=T, tick=T)  # tweak for GKS3, change ylim in plot too
axis(2, labels=T, tick=T)

# Create surv data set
nest_mort <- Surv(time = BC_nest_dat$entry, 
                  time2 = BC_nest_dat$exit, 
                  event = BC_nest_dat$event)

fled_mort <- Surv(time = BC_fled_dat$entry, 
                  time2 = BC_fled_dat$exit, 
                  event = BC_fled_dat$event)

cox.BC <- coxme(Surv(pref_age, pref_status) ~ sex + hatch_order + 
                  (1|year) + (1|nest_ID), data =  dat_BC) # prefledging survival
summary(cox.BC)

reg_fit_me <- coxme(formula = mort ~ sex + (1|year) + (1|nest_ID), data =  dat_BC)
summary(reg_fit)

nest_fit <- coxph(formula = nest_mort ~ BC_nest_dat$sex)
summary(nest_fit)

fled_fit <- coxph(formula = fled_mort ~ BC_fled_dat$sex)
summary(fled_fit)$coef[1]

# Now get baseline curve
baseline <- basehaz(fled_fit)

# Draw baseline hazard (that's female)
plot(baseline$time, 1 - baseline$hazard, type='l',main="Hazard rates") 

# Draw male hazard
lines(baseline$time, 1 - exp(summary(fled_fit)$coef[1])*baseline$hazard, col="blue") 

data(oldmort) #create the data

library(eha) #used for data 

# Create surv data set
mort <- Surv(time=oldmort$enter,time2=oldmort$exit,event=oldmort$event)

reg_fit <- coxph(formula=mort~oldmort$sex)
summary(reg_fit)

# Now get baseline curve
baseline <- basehaz(reg_fit)

# Draw baseline hazard (that's male)
plot(baseline$time, baseline$hazard, type='l',main="Hazard rates") 

# Draw female hazard
lines(baseline$time, exp(-0.1929)*baseline$hazard, col="blue") 

# cbind variables of interest in this analysis
X = cbind(dat$species, dat$lab_no, dat$sex, dat$pref_age, dat$pref_status, 
          dat$blank_space, dat$postf_age, dat$postf_status,  dat$ageC, 
          dat$statusC, dat$year, dat$nest_ID, dat$lay_date, dat$hatch_order, 
          dat$clutch_size, dat$brood_size, dat$No_fledgelings, dat$Fledge_age, 
          dat$Fledge_mass, dat$Fledge_tarsus, dat$Fledge_FS, dat$Fledge_FSO, 
          dat$paternity) 

# set events and states as.numeric 
# NB: model will not run if states & ages are as.factors; I had this problem 
# before with NAs defaulting states & ages as factors
dat$postf_age <- as.numeric(dat$postf_age)
dat$postf_status <- as.numeric(dat$postf_status)
dat$hatch_order <- as.numeric(dat$hatch_order)

# "Fledged_status" indicates if the individual fledged or not (no = died in nest, yes = survived nest phase)
# "pref_age" indicates when the individual died (in the case of "Fledged_status" == "no") or when it fledged (in the case of "Fledged_status" == "yes")

#### Brett's example data ----
pc <- read.table("literature/from_Brett/enchist.txt",header=T,sep=",")
# look at data
# name = bird band number plus radio frequency, year = year of study (1-5)
# trt = pre/post-construction, entry and exit = week of study period
# event = 1 for died, 0 for survived, dturb = distance of home range to turbine in km
head(pc, 20) # first 20 lines of file

dat_BC <-
  dat %>% 
  filter(species == "BC") %>% 
  filter(Fledged_status != "<NA>") %>% 
  # filter(nest_ID == "2016_33")
  select(nest_ID, lab_no, ring_ID, year, sex, pref_age, Fledged_status, hatch_order) %>% 
  mutate(ind_ID = paste(nest_ID, lab_no, ring_ID, sep = "_"),
         entry = 0,
         event = ifelse(Fledged_status == "yes", 0, 1)) %>% 
  dplyr::rename(exit = pref_age) %>% 
  select(ind_ID, nest_ID, year, sex, entry, exit, event, hatch_order)

# Kaplan-Meier comparisons
# trt comparison between pre and post, pooling across years
# mfit <- survfit(Surv(pc$entry, pc$exit, pc$event) ~ 1 )
mfit <- survfit(Surv(dat_BC$entry, dat_BC$exit, dat_BC$event) ~ dat_BC$sex )
print(survfit(Surv(dat_BC$entry, dat_BC$exit, dat_BC$event) ~ dat_BC$sex), print.rmean=TRUE)
summary(mfit)

# Plot the Kaplan-Meier cumulative survival
# xaxs, yaxs=i gets rid of rabbit holes at corner of graphs
plot(mfit, conf.int=T, lty=c(1,2,2), col=c("dark green", "red"), xaxs="i", yaxs="i", xlim=c(0, 53), lwd=2.5, xlab=c("week"), ylab=c("survival function"))
legend("topright", legend = c("Post", "Pre"), cex=2, lty=1, col=c("dark green", "red"), lwd=2.5, bty="n")

trt.prepost <- coxph(Surv(dat_BC$entry, dat_BC$exit, dat_BC$event)~ dat_BC$sex + cluster(dat_BC$ind_ID))	 # Interactive model 
summary(trt.prepost)		    # Summary of the model			
cox.zph(trt.prepost)              # Model diagnostic; testing the proportional hazard assumption (see "?cox.zph")		
windows()
plot(cox.zph(trt.prepost, transform="identity"))   # 

# test of main effects models for treatment period, distance to turbine, cluster is a random effects for individual bird
trt.prepost <- coxph(Surv(dat_BC$entry, dat_BC$exit, dat_BC$event)~ dat_BC$sex + cluster(dat_BC$ind_ID))	 # Additive model 
summary(trt.prepost)		    # Summary of the model	

# hazard functions by treatment period
# note Surv() in the sshzd() function takes the arguments in a different order than Surv() in the survival package
# hazard function = hz
hz.M <- sshzd(Surv(exit, event, entry)~exit, data=filter(dat_BC, sex == "M")) #alpha=0.25 will overfit
hz.F <- sshzd(Surv(exit, event, entry)~exit, data=filter(dat_BC, sex == "F"))
hz.all <- sshzd(Surv(exit, event, entry)~exit, data=pc, alpha=1.2)

# calculate the hazard rate function for plotting for a sequence of time (tt = 1 to 52 weeks by 0.5)
# hazard = hz
tt <- seq(1,15,0.5)

# predicted hazard
# pre-construction
est <- hzdrate.sshzd(hz.M,tt,se=TRUE)
hz.M.fit = 1 - est$fit
hz.M.hi <- 1 - est$fit*exp(1.96*est$se)
hz.M.lo <- 1 - est$fit/exp(1.96*est$se)

# post-construction
est <- hzdrate.sshzd(hz.F,tt,se=TRUE)
hz.F.fit = 1 - est$fit
hz.F.hi <- 1 - est$fit*exp(1.96*est$se)
hz.F.lo <- 1 - est$fit/exp(1.96*est$se)

# pooled
est <- hzdrate.sshzd(hz.all,tt,se=TRUE)
hzall.fit = 1 - est$fit
hzall.hi <- 1 - est$fit*exp(1.96*est$se)
hzall.lo <- 1 - est$fit/exp(1.96*est$se)

# graph for hazard functions in pre and post periods
windows()  # new graphics window				
plot(tt, hz.M.fit, type="l", lty=1, lwd=2, col="red", xaxs="i", yaxs="i", 
     xlim=c(1, 15), ylim=c(0.7, 1), 
     xlab=c("days since hatching"), ylab=c("hazard function"))
lines(tt, hz.M.hi, lty=2, lwd=2, col="red")
lines(tt, hz.M.lo, lty=2, lwd=2, col="red")
lines(tt, hz.F.fit, lty=1, lwd=2, col="dark green")
lines(tt, hz.F.hi, lty=2, lwd=2, col="dark green")
lines(tt, hz.F.lo, lty=2, lwd=2, col="dark green")
legend("topright", legend = c("Post", "Pre"), cex=2, lty=c(1,1), col=c("dark green", "red"), lwd=2.5, bty="n")
axis(1, labels=F, tick=T)
# axis(2, at=seq(0,0.1,0.01), labels=T, tick=T)  # tweak for GKS3, change ylim in plot too
axis(2, labels=T, tick=T)


##------------------------------------------------------------------
# subset data by species (i.e. BC, WBC & CTC) for survival analysis
#-------------------------------------------------------------------
dat_BC = subset(dat, species=="BC")   ## subset data for BC 
str(dat_BC)
summary(dat_BC)

dat_WBC = subset(dat, species=="WBC")   ## subset data for WBC 
str(dat_WBC)
summary(dat_WBC)

dat_BC <- mutate

# prefledging survival
hzd_BC_nest <- sshzd(Surv(pref_age, pref_status) ~ sex + hatch_order + 
                       (1|year) + (1|nest_ID), 
                     data = dat_BC, alpha=1.2)
summary(hzd_BC_nest)

# postfledging survival
hzd_BC_fled <- sshzd(Surv(postf_age, postf_status) ~ sex  + hatch_order + 
                       (1|year) + (1|nest_ID), data =  dat_BC)
summary(hzd_BC_fled)

# combined pre-&postfledging survival
hzd_BC_all <- sshzd(Surv(ageC, statusC) ~ sex + hatch_order + 
                       (1|year) + (1|nest_ID), data =  dat_BC)
summary(hzd_BC_fled)

sshzd

cox.BC <- sshzd(Surv(pref_age, pref_status) ~ sex + hatch_order + 
                  (1|year) + (1|nest_ID), data =  dat_BC)

# calculate the daily survival as 1 - the daily hazard rate
hzdrate.sshzd()