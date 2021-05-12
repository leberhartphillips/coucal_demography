
#--------------------------
# Clear working environment
#--------------------------

rm(list=ls())

#--------------------------
## Load packages
#--------------------------
library(survival)
library(ggplot2)
library(magrittr)
library(ggpubr)
library(survminer) 
library(bdsmatrix)
library(coxme) 
library(epitools) 

#-----------------------------------------
# set working directory
#-----------------------------------------
# getwd()
# setwd("C:/Users/isafari/86B3F2D6AC2B00178AE1F2F77F781B0C.TMP/Documents/MY STATISTICAL ANALYSES/2019/Coucal_survival_2020v") ##change as required

#------------------------------------------------------------------------
## read the data into R
## pre & postfleding survival combined (variable postf_age, postf_status) 
#------------------------------------------------------------------------
list.files()
dat = read.csv("data/raw/Coucal_chick_survival_2001-2019_20200129.csv", header=T)
head(dat)
str(dat)
#-------------------------------------------------------------------
####cbind variables of interest in this analysis
#-------------------------------------------------------------------
X = cbind(dat$species, dat$lab_no, dat$sex, dat$pref_age, dat$pref_status, 
          dat$blank_space, dat$postf_age, dat$postf_status,  dat$ageC, 
          dat$statusC, dat$year, dat$nest_ID, dat$lay_date, dat$hatch_order, 
          dat$clutch_size, dat$brood_size, dat$No_fledgelings, dat$Fledge_age, 
          dat$Fledge_mass, dat$Fledge_tarsus, dat$Fledge_FS, dat$Fledge_FSO, 
          dat$paternity) 

#-------------------------------------------------------------------
# set events and states as.numeric 
# NB: model will not run if states & ages are as.factors; I had this problem 
# before with NAs defaulting states & ages as factors
#--------------------------------------------------------------------
dat$postf_age <- as.numeric(dat$postf_age)
dat$postf_status <- as.numeric(dat$postf_status)
dat$hatch_order <- as.numeric(dat$hatch_order)

##------------------------------------------------------------------
# subset data by species (i.e. BC, WBC & CTC) for survival analysis
#-------------------------------------------------------------------
dat_BC = subset(dat, species=="BC")   ## subset data for BC 
str(dat_BC)
summary(dat_BC)

dat_WBC = subset(dat, species=="WBC")   ## subset data for WBC 
str(dat_WBC)
summary(dat_WBC)

dat_CTC = subset(dat, species=="CTC")   ## subset data for CTC 
str(dat_CTC)
summary(dat_CTC)

## -----------------------------------------------------------------------
## Analysis of sex specific survival (pre-, post-fledging and combined pre+postf)
# cox model with random effects 
## -----------------------------------------------------------------------
# BC data
cox.BC <- coxme(Surv(pref_age, pref_status) ~ sex + hatch_order + 
                  (1|year) + (1|nest_ID), data =  dat_BC) # prefledging survival
summary(cox.BC)

cox.BC <- coxme(Surv(postf_age, postf_status) ~ sex  + hatch_order + 
                  (1|year) + (1|nest_ID), data =  dat_BC)  # postfledging survival
summary(cox.BC)

cox.BC <- coxme(Surv(ageC, statusC) ~ sex + hatch_order + 
                  (1|year) + (1|nest_ID), data =  dat_BC)  # combined pre-&postfledging survival
summary(cox.BC)

# WBC data
cox.WBC <- coxme(Surv(pref_age, pref_status) ~ sex + hatch_order + 
                   (1|year) + (1|nest_ID), data =  dat_WBC) # prefledging survival
summary(cox.WBC)

cox.WBC <- coxme(Surv(postf_age, postf_status) ~ sex + hatch_order + 
                   (1|year) + (1|nest_ID), data =  dat_WBC)  # postfledging survival
summary(cox.WBC)

cox.WBC <- coxme(Surv(ageC, statusC) ~ sex + hatch_order + 
                   (1|year) + (1|nest_ID), data =  dat_WBC)  # combined pre-& postfledging survival
summary(cox.WBC)

# # CTC data
# cox.CTC <- coxme(Surv(pref_age, pref_status) ~ sex + hatch_order  + (1|year) + (1|nest_ID), data =  dat_CTC) # prefledging survival
# summary(cox.CTC)

#------------------------------------------------------
# looking more into the survival parameters
#------------------------------------------------------
cox.BC # print(cox.BC)
fixef(cox.BC) # fixed efects coefficients
ranef(cox.BC) # random effects coefficients
vcov(cox.BC)  # fixed efects variance matrix
VarCorr(cox.BC) # random effects parameters (theta)

cox.WBC # print(cox.WBC)
fixef(cox.WBC) # fixed efects coefficients
ranef(cox.WBC) # random effects coefficients
vcov(cox.WBC)  # fixed efects variance matrix
VarCorr(cox.WBC) # random effects parameters (theta)

cox.CTC # print(cox.CTC)
fixef(cox.CTC) # fixed efects coefficients
ranef(cox.CTC) # random effects coefficients
vcov(cox.CTC)  # fixed efects variance matrix
VarCorr(cox.CTC) # random effects parameters (theta)

#--------------------------------------------------------------
######### create a survival object by sexes in BC
#--------------------------------------------------------------

bcfit.by.sex <- survfit(Surv(postf_age, postf_status) ~ sex, data = dat_BC) 
bcfit.by.sex

bcfit.by.sex <- survfit(Surv(ageC, statusC) ~ sex, data = dat_BC) 
bcfit.by.sex

prob_BC <- 
  summary(bcfit.by.sex, times = c(0:109)) # survival probabilities for males & females from day 1 until day 109 (the last age in BC data)

test <- 
  cbind(prob_BC$time, prob_BC$surv, prob_BC$lower, prob_BC$upper, 
        c(rep("F", 110), rep("M", 95))) %>% 
  as_data_frame() %>% 
  dplyr::rename(age = V1,
         survival = V2,
         lowerCI = V3,
         upperCI = V4,
         sex = V5) %>% 
  dplyr::mutate(age = as.numeric(age),
         survival = as.numeric(survival),
         lowerCI = as.numeric(lowerCI),
         upperCI = as.numeric(upperCI),
         sex = as.factor(sex))

ggplot(data = filter(test, age != 0), aes(x = age, y = survival, group = sex)) +
  geom_point(aes(color = sex))  +
  # stat_smooth(method = lm, formula = y ~ poly(x, 3, raw = TRUE)) +
  geom_smooth(method = lm, formula = y ~ log(x)) +
  # stat_smooth(method = gam, formula = y ~ s(x)) +
  facet_wrap(sex ~ .)

age_data <- data.frame(age = c(1:110))

model_F <- lm(survival ~ log(age), data = filter(test, sex == "F" & age != 0))
model_M <- lm(survival ~ log(age), data = filter(test, sex == "M" & age != 0))
predictions_F <- 
  model_F %>% 
  predict(age_data) %>% 
  as_data_frame() %>% 
  dplyr::rename(cum_surv = value) %>% 
  dplyr::mutate(cum_surv_trans = cum_surv - 0.136,
         sex = "F")

predictions_M <- 
  model_M %>% 
  predict(age_data) %>% 
  as_data_frame() %>% 
  dplyr::rename(cum_surv = value) %>% 
  dplyr::mutate(cum_surv_trans = cum_surv - 0.113,
         sex = "M")

predictions <- 
  bind_rows(predictions_F,
            predictions_M)

summary(bcfit.by.sex)$table # gives summaries of events,
prob_BC <- summary(bcfit.by.sex) ### gives detailed summaries of events, number at risk, survival estimates & confidence intervals

ggplot()

surv_diff_BC_sexes <- survdiff(Surv(postf_age, postf_status) ~ sex , data = dat_BC) # use Log-Rank test to compare survival for males and females BC
surv_diff_BC_sexes

surv_diff_BC_sexes_strata_hatch_order <- survdiff(Surv(postf_age, postf_status) ~ sex + strata(hatch_order), data = dat_BC) # use Log-Rank test to compare survival for males and females BC statified by hatching order
surv_diff_BC_sexes_strata_hatch_order

#-----------------------------------------------------------------------
# create survival object by sex in WBC
#-----------------------------------------------------------------------
# wbcfit.by.sex <- survfit(Surv(postf_age, postf_status) ~ sex + hatch_order + (1|year), data = dat_WBC)
# wbcfit.by.sex

wbcfit.by.sex <- survfit(Surv(postf_age, postf_status) ~ sex + (1|year), data = dat_WBC)
wbcfit.by.sex

prob_WBC <- summary(wbcfit.by.sex, times = c(0:125)) # survival probabilities for males & females from day 1 until day 125 (the last age in WBC data)
prob_WBC

summary(wbcfit.by.sex)$table # gives summaries of events,
summary(wbcfit.by.sex) ### gives detailed summaries of events, number at risk, survival estimates & confidence intervals

surv_diff_WBC_sexes <- survdiff(Surv(postf_age, postf_status) ~ sex , data = dat_WBC) # use Log-Rank test to compare survival for males and females WBC
surv_diff_WBC_sexes

surv_diff_WBC_sexes_strata_hatch_order <- survdiff(Surv(postf_age, postf_status) ~ sex + strata(hatch_order), data = dat_WBC) # use Log-Rank test to compare survival for males and female WBC statified by hatching order
surv_diff_WBC_sexes_strata_hatch_order


# ------------------------------------------------------------------------
## Plotting sex-specific survival curves
#-------------------------------------------------------------------------
## plotting sex specific survival curves for BC using ggsurvplot 
# ------------------------------------------------------------------------
bcfit.by.sex <- survfit(Surv(postf_age, postf_status) ~ sex , data = dat_BC)
bcfit.by.sex
bcfit.by.sex <- survfit(Surv(ageC, statusC) ~ sex , data = dat_BC)

ggsurvplot(bcfit.by.sex, 
           break.time.by = 10, 
           risk.table = TRUE,
           ylab = "Cumulative Survival (± 95% CI)",
           xlab = "Time since hatching (days)",
           conf.int = TRUE) #  plot the survival## ok



BC_plot1 <- ggsurvplot(bcfit.by.sex,     # survfit object with calculated statistics.
                       risk.table = F,  # "abs_pct"show number and percent of the individuals at at risk
                       ncensor.plot = F,
                       pval = TRUE,  # show p-value of log-rank test.
                       pval.coord = c(0, 0.2),
                       pval.size = 7, # numeric value specifying the p-value text size. Default is 5.
                       pval.method = TRUE,
                       pval.method.coord = c(0, 0.250),
                       pval.method.size = 7,
                       conf.int = T,  # show confidence intervals for point estimaes of survival curves.
                       ylab = "Cumulative Survival",
                       xlab = "Time (days)",
                       xlim = c(0,110),    # present narrower X axis, but not affect survival estimates.
                       break.time.by = 10, 
                       legend =  c(0.11, 0.4), # legend = "top", # legend =  c(x, y)
                       font.legend = c(20, "plain", "black"),
                       legend.title = " ", # "Sex",
                       legend.labs = c("Female (N = 66)", "Male     (N = 59)")) 
BC_plot1

BC_plot1 <- 
  ggpar(BC_plot1,
  font.title    = c(20, "bold", "black"),         
  font.subtitle = c(20, "bold", "black"),
  font.caption  = c(20, "plain", "black"),        
  font.x        = c(20, "bold", "black"), 
  font.y        = c(20, "bold", "black"), 
  font.xtickslab = c(20, "plain", "black"),
  font.ytickslab = c(20, "plain", "black"),
  legend =  c(0.11, 0.4), # legend = "top", # legend =  c(x, y)
  font.legend = c(20, "plain", "black"))

BC_plot1

# -------------------------------------------------------------------------------
# plotting survival curves for all WBC by sex using ggsurvplot 
#--------------------------------------------------------------------------------
wbcfit.by.sex <- survfit(Surv(postf_age, postf_status) ~ sex , data = dat_WBC)
wbcfit.by.sex
wbcfit.by.sex <- survfit(Surv(ageC, statusC) ~ sex , data = dat_WBC)

ggsurvplot(wbcfit.by.sex, break.time.by = 10, risk.table = TRUE) #  plot the survival## ok

WBC_plot1 <- ggsurvplot(wbcfit.by.sex,     # survfit object with calculated statistics.
                       risk.table = F,  # "abs_pct"show number and percent of the individuals at at risk
                       ncensor.plot = F,
                       pval = TRUE,  # show p-value of log-rank test.
                       pval.coord = c(0, 0.2),
                       pval.size = 7, # numeric value specifying the p-value text size. Default is 5.
                       pval.method = TRUE,
                       pval.method.coord = c(0, 0.250),
                       pval.method.size = 7,
                       conf.int = T,  # show confidence intervals for point estimaes of survival curves.
                       ylab = "Survival probability",
                       xlab = "Time (days)",
                       xlim = c(0,126),    # present narrower X axis, but not affect survival estimates.
                       break.time.by = 10, 
                       legend =  c(0.11, 0.4), # legend = "top", # legend =  c(x, y)
                       font.legend = c(20, "plain", "black"),
                       legend.title = " ", # "Sex",
                       legend.labs = c("Female (N = 50)", "Male     (N = 57)")) 
WBC_plot1

WBC_plot1 <- ggpar(WBC_plot1,
  font.title    = c(20, "bold", "black"),         
  font.subtitle = c(20, "bold", "black"),
  font.caption  = c(20, "plain", "black"),        
  font.x        = c(20, "bold", "black"), 
  font.y        = c(20, "bold", "black"), 
  font.xtickslab = c(20, "plain", "black"),
  font.ytickslab = c(20, "plain", "black"),
  legend =  c(0.11, 0.4), # legend = "top", # legend =  c(x, y)
  font.legend = c(20, "plain", "black"))

WBC_plot1
# -------------------------------------------------------------------------------
## plotting survival curves for all CTC nestlings by sex using ggsurvplot 
# -------------------------------------------------------------------------------
ctcfit.by.sex <- survfit(Surv(pref_age, pref_status) ~ sex , data = dat_CTC)
ctcfit.by.sex

ggsurvplot(ctcfit.by.sex, break.time.by = 2, risk.table = TRUE) #  plot the survival## ok

CTC_plot1 <- ggsurvplot(ctcfit.by.sex,     # survfit object with calculated statistics.
                        risk.table = F,  # "abs_pct"show number and percent of the individuals at at risk
                        ncensor.plot = F,
                        pval = TRUE,  # show p-value of log-rank test.
                        pval.coord = c(0, 0.2),
                        pval.size = 7, # numeric value specifying the p-value text size. Default is 5.
                        pval.method = TRUE,
                        pval.method.coord = c(0, 0.250),
                        pval.method.size = 7,
                        conf.int = T,  # show confidence intervals for point estimaes of survival curves.
                        ylab = "Survival probability",
                        xlab = "Time (days)",
                        xlim = c(0,21),    # present narrower X axis, but not affect survival estimates.
                        break.time.by = 2, 
                        legend =  c(0.095, 0.4), # legend = "top", # legend =  c(x, y)
                        font.legend = c(20, "plain", "black"),
                        legend.title = " ", # "Sex",
                        legend.labs = c("Female (N = 11)", "Male     (N = 10)")) 
CTC_plot1

CTC_plot1 <- ggpar(CTC_plot1,
                   font.title    = c(20, "bold", "black"),         
                   font.subtitle = c(20, "bold", "black"),
                   font.caption  = c(20, "plain", "black"),        
                   font.x        = c(20, "bold", "black"), 
                   font.y        = c(20, "bold", "black"), 
                   font.xtickslab = c(20, "plain", "black"),
                   font.ytickslab = c(20, "plain", "black"),
                   legend =  c(0.095, 0.4), # legend = "top", # legend =  c(x, y)
                   font.legend = c(20, "plain", "black"))

CTC_plot1
##########################################
# END
#----------------------------------------






