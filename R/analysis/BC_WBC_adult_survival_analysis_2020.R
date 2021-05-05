
rm(list=ls())

## Load survival package
library(survival)
library(ggplot2)# install.packages("ggplot2")
library(magrittr)
library(ggpubr)
library(survminer) #install.packages("survminer")
library(bdsmatrix)
library(coxme) ##install.packages("coxme")
library(epitools)#install.packages("epitools")


getwd()
setwd("C:/Users/isafari/86B3F2D6AC2B00178AE1F2F77F781B0C.TMP/Documents/MY STATISTICAL ANALYSES/2019/Coucal_survival_2020v") ##change as required

list.files()
dat = read.csv("data/raw/Coucal_adults_survival_2001-2019_20200129.csv", header=T)
head(dat)
str(dat)


#-----------------------------------------------
#cbind the variables of interest in this analysis
#-----------------------------------------------
X = cbind(dat$year, dat$species, dat$Alu, dat$sex, dat$body_mass, dat$right_tarsus, dat$bill_length, dat$wing_length, dat$tail_body, dat$age1, dat$status1, dat$ageC, dat$statusC) # bind sex and species


#--------------------------------------------------
# subset data by species i.e. BC, WBC & CTC
#--------------------------------------------------
dat_BC = subset(dat, species== "BC")   ## subset data for BC 
str(dat_BC)
summary(dat_BC)

dat_WBC = subset(dat, species=="WBC")   ## subset data for WBC 
str(dat_WBC)
summary(dat_WBC)

dat_CTC = subset(dat, species=="CTC")   ## subset data for CTC 
str(dat_CTC)
summary(dat_CTC)

## ------------------------------------------------
## Analysis of sex specific mortality
## ------------------------------------------------
# cox model with random effects (year)

# BC adults
cox.BC <- coxme(Surv(ageC, statusC) ~ sex + (1|year), data =  dat_BC)
summary(cox.BC)


# WBC adults
cox.WBC <- coxme(Surv(ageC, statusC) ~ sex + (1|year), data =  dat_WBC)
summary(cox.WBC)


# CTC adults (note that this model will not run because there is no any mortality observed in CTCs) 
# cox.CTC <- coxme(Surv(ageC, statusC) ~ sex + (1|year), data =  dat_CTC)
# summary(cox.CTC)

#----------------------------------------------------
# looking more into the effects sizes and  (co)variances
cox.BC # print(cox.BC)
fixef(cox.BC) # fixed efects coefficients
ranef(cox.BC) # random effects coefficients
vcov(cox.BC)  # fixed efects variance matrix
vcov(cox.WBC)  # fixed efects variance matrix
VarCorr(cox.BC) # random effects parameters (theta)
VarCorr(cox.WBC)  # random effects parameters (theta)

#--------------------------------------------------------------
######### create a survival object by sexes in BC
#--------------------------------------------------------------
bcfit.by.sex <- survfit(Surv(ageC, statusC) ~ sex + (1|year), data = dat_BC)
bcfit.by.sex

bcfit.by.sex <- survfit(Surv(ageC, statusC) ~ sex, data = dat_BC) 
bcfit.by.sex

 prob_BC <- summary(bcfit.by.sex, times = c(0:1127)) # survival probabilities for males & females from day 1 until day 1127 (the last age in BC data)
 prob_BC

 summary(bcfit.by.sex)$table # gives summaries of events,
 summary(bcfit.by.sex) ### gives detailed summaries of events, number at risk, survival estimates & confidence intervals

surv_diff_BC_sexes <- survdiff(Surv(ageC, statusC) ~ sex , data = dat_BC) # use Log-Rank test to compare survival for males and females BC
surv_diff_BC_sexes

surv_diff_BC_body_mass <- survdiff(Surv(ageC, statusC) ~ body_mass, data = dat_BC) # use Log-Rank test to compare survival by body_mass BC
surv_diff_BC_body_mass

surv_diff_BC_sexes_body_mass <- survdiff(Surv(ageC, statusC) ~ sex + body_mass, data = dat_BC) # use Log-Rank test to compare survival by sex and body_mass
surv_diff_BC_sexes_body_mass


surv_diff_BC_sexes_strata_body_mass <- survdiff(Surv(ageC, statusC) ~ sex + strata(body_mass), data = dat_BC) # use Log-Rank test to compare survival for males and females BC statified by hatching order
surv_diff_BC_sexes_strata_body_mass


#-----------------------------------------------------------------------
# create survival object by sex in WBC
#-----------------------------------------------------------------------
wbcfit.by.sex <- survfit(Surv(ageC, statusC) ~ sex + (1|year), data = dat_WBC)
wbcfit.by.sex

# wbcfit.by.sex <- survfit(Surv(ageC, statusC) ~ sex, data = dat_WBC) 
# wbcfit.by.sex

prob_WBC <- summary(wbcfit.by.sex, times = c(0:2055)) # survival probabilities for males & females from day 1 until day 2055 (the last age in WBC data)
prob_WBC

summary(wbcfit.by.sex)$table # gives summaries of events,
summary(wbcfit.by.sex) ### gives detailed summaries of events, number at risk, survival estimates & confidence intervals

surv_diff_WBC_sexes <- survdiff(Surv(ageC, statusC) ~ sex , data = dat_WBC) # use Log-Rank test to compare survival for males and females BC
surv_diff_WBC_sexes


surv_diff_WBC_body_mass <- survdiff(Surv(ageC, statusC) ~ body_mass, data = dat_WBC) # use Log-Rank test to compare survival by body_mass BC
surv_diff_WBC_body_mass

surv_diff_WBC_sexes_body_mass <- survdiff(Surv(ageC, statusC) ~ sex + body_mass, data = dat_WBC) # use Log-Rank test to compare survival by sex and body_mass
surv_diff_WBC_sexes_body_mass


surv_diff_WBC_sexes_strata_body_mass <- survdiff(Surv(ageC, statusC) ~ sex + strata(body_mass), data = dat_WBC) # use Log-Rank test to compare survival for males and females BC statified by hatching order
surv_diff_WBC_sexes_strata_body_mass

#--------------------------------------------------------------------------------
## create survival object by sex in CTC
#This model will not run properly because there are no mortality in the CTC data
#---------------------------------------------------------------------------------
ctcfit.by.sex <- survfit(Surv(ageC, statusC) ~ sex + (1|year), data = dat_CTC)
ctcfit.by.sex

# ctcfit.by.sex <- survfit(Surv(ageC, statusC) ~ sex, data = dat_CTC) 
# ctcfit.by.sex

prob_CTC <- summary(ctcfit.by.sex, times = c(0:668)) # survival probabilities for males & females from day 1 until day 668 (the last age in CTC data)
prob_CTC

summary(ctcfit.by.sex)$table # gives summaries of events,
summary(ctcfit.by.sex) ### gives detailed summaries of events, number at risk, survival estimates & confidence intervals

surv_diff_CTC_sexes <- survdiff(Surv(ageC, statusC) ~ sex , data = dat_CTC) # use Log-Rank test to compare survival for males and females BC
surv_diff_CTC_sexes


surv_diff_CTC_body_mass <- survdiff(Surv(ageC, statusC) ~ body_mass, data = dat_CTC) # use Log-Rank test to compare survival by body_mass BC
surv_diff_CTC_body_mass

surv_diff_CTC_sexes_body_mass <- survdiff(Surv(ageC, statusC) ~ sex + body_mass, data = dat_CTC) # use Log-Rank test to compare survival by sex and body_mass
surv_diff_CTC_sexes_body_mass

surv_diff_CTC_sexes_strata_body_mass <- survdiff(Surv(ageC, statusC) ~ sex + strata(body_mass), data = dat_CTC) # use Log-Rank test to compare survival for males and females BC statified by hatching order
surv_diff_CTC_sexes_strata_body_mass

# --------------------------------------------------------------------------
## 
## plotting sex specific adult survival curves for BC using ggsurvplot 
## -------------------------------------------------------------------------
bcfit.by.sex <- survfit(Surv(ageC, statusC) ~ sex + (1|year), data = dat_BC)
bcfit.by.sex

ggsurvplot(bcfit.by.sex, break.time.by = 100, risk.table = TRUE) #  plot the survival## ok


BC_plot1 <- ggsurvplot(bcfit.by.sex,     # survfit object with calculated statistics.
 risk.table = F,  # "abs_pct"show number and percent of the individuals at at risk
 ncensor.plot = F,
 pval = TRUE,  # show p-value of log-rank test.
 pval.coord = c(0, 0.45),
 pval.size = 7, # numeric value specifying the p-value text size. Default is 5.
 pval.method = TRUE,
 pval.method.coord = c(0, 0.50),
 pval.method.size = 7,
 conf.int = T,  # show confidence intervals for point estimaes of survival curves.
 ylab = "Survival probability",
 xlab = "Time (days)",
 xlim = c(0,1130),    # present narrower X axis, but not affect survival estimates.
 break.time.by = 100, 
 legend =  c(0.15, 0.8), # legend = "top", # legend =  c(x, y)
 font.legend = c(20, "plain", "black"),
 legend.title = " ", # "Sex",
legend.labs = c("Female (N = 196)", "Male     (N = 167)")) 
BC_plot1

BC_plot1 <- ggpar(BC_plot1,
  font.title    = c(20, "bold", "black"),         
  font.subtitle = c(20, "bold", "black"),
  font.caption  = c(20, "plain", "black"),        
  font.x        = c(20, "bold", "black"), 
  font.y        = c(20, "bold", "black"), 
  font.xtickslab = c(20, "plain", "black"),
  font.ytickslab = c(20, "plain", "black"),
    legend =  c(0.15, 0.8), # legend = "top", # legend =  c(x, y)
  font.legend = c(20, "plain", "black"))

BC_plot1


# -------------------------------------------------------------------------------
# WBC
###### plotting survival curves for all WBC adults based on sex using ggsurvplot 
#---------------------------------------------------------------------------------
wbcfit.by.sex <- survfit(Surv(ageC, statusC) ~ sex + (1|year), data = dat_WBC)
wbcfit.by.sex

ggsurvplot(wbcfit.by.sex, break.time.by = 100, risk.table = TRUE)  ## ok

WBC_plot1 <- ggsurvplot(wbcfit.by.sex,               # survfit object with calculated statistics.
  risk.table = F,   #"abs_pct" show number and percent of individuals at risk
  #risk.table = TRUE,       # show risk table.
  ncensor.plot = F, 
  pval = TRUE,    # show p-value of log-rank test.
  pval.coord = c(0, 0.45),
  pval.size = 7, # numeric value specifying the p-value text size. Default is 5.
  pval.method = TRUE,
  pval.method.coord = c(0, 0.50),
  pval.method.size = 7,
  conf.int = T, # show confidence intervals for point estimaes of survival curves.
  ylab = "Survival probability",
  xlab = "Time (days)",
  xlim = c(0,2060),        # present narrower X axis, but not affect survival estimates.
  break.time.by = 100, 
  legend =  c(0.15, 0.8), # legend = "top", # legend =  c(x, y)
  font.legend = c(20, "plain", "black"),
  legend.title = " ", # "Sex",
  legend.labs = c("Female (N = 60)", "Male     (N = 93)")) 
WBC_plot1

WBC_plot1 <- ggpar(WBC_plot1,
  font.title    = c(20, "bold", "black"),         
  font.subtitle = c(20, "bold", "black"), 
  font.caption  = c(20, "plain", "black"),        
  font.x        = c(20, "bold", "black"),          
  font.y        = c(20, "bold", "black"),      
  font.xtickslab = c(20, "plain", "black"),
  font.ytickslab = c(20, "plain", "black"),
  legend.title = " ", # "Sex",
  legend =  c(0.15, 0.8), # legend = "top", # legend =  c(x, y)
  font.legend = c(20, "plain", "black"))

WBC_plot1

# -------------------------------------------------------------------------------
## plotting survival curves for all CTC nestlings based on sex using ggsurvplot 
# -------------------------------------------------------------------------------


ggsurvplot(pref_ctcfit.by.sex, break.time.by = 2, risk.table = TRUE)  ## ok

CTC_plot1 <- ggsurvplot(pref_ctcfit.by.sex,      # survfit object with calculated statistics.
  risk.table = F,   #"abs_pct",     #show number and percent of the individuals at begining who are at risk
  #risk.table = TRUE,       # show risk table.
  ncensor.plot = F,
  pval = TRUE,      # show p-value of log-rank test.
  pval.coord = c(0, 0.20),
  pval.size = 7, # numeric value specifying the p-value text size. Default is 5.
  pval.method = TRUE,
  pval.method.coord = c(0, 0.25),
  pval.method.size = 7,
  conf.int = T,  # show confidence intervals for point estimaes of survival curves.
  ylab = "Survival probability",
  xlab = "Time (days)",
  xlim = c(0, 22),        # present narrower X axis, but not affect survival estimates.
  break.time.by = 2,
  font.legend = c(20, "plain", "black"),
  legend.title = " ", # "Sex",
  legend.labs = c("Female (N = 07)", "Male     (N = 10)")) 

CTC_plot1 <- ggpar(CTC_plot1,
    font.title    = c(20, "bold", "black"),         
    font.subtitle = c(20, "bold", "black"), 
    font.caption  = c(20, "plain", "black"),        
    font.x        = c(20, "bold", "black"),          
    font.y        = c(20, "bold", "black"),      
    font.xtickslab = c(20, "plain", "black"),
    font.ytickslab = c(20, "plain", "black"),
    legend =  c(0.15, 0.5), # legend = "top", # legend =  c(x, y)
    legend.title = " ", # "Sex",
    font.legend = c(20, "plain", "black"))

CTC_plot1




##########################################
##########################################
# ---------------------------------------------
## computing daily survival probabilities
prob_prefBC <- summary(pref_bcfit.by.sex, times = c(0:22)) # survival probabilities for males & females from day 1 until day 22
prob_prefBC


prob_prefWBC <- summary(pref_wbcfit.by.sex, times = c(0:22)) # survival probabilities for males & females from day 1 until day 22
prob_prefWBC


prob_prefCTC <- summary(pref_ctcfit.by.sex, times = c(0:22)) # survival probabilities for males & females from day 1 until day 22
prob_prefCTC

# ---------------------------------------------------------










## --------------------------------------------
# plotting sex ratios and 95% CIs
# ------------------------------------------------
# pre-fleding SRs
# ----------------------------
list.files()
#dat = read.csv("Plotting_pref_juve_SR_BC_WBC_mean.csv", header=T)
dat = read.csv("Plotting_postf_juve_SR_BC_WBC_mean.csv", header=T)
head(dat)
str(dat)

####cbind the variables of interest in this analysis
X = cbind(dat$age, dat$meanSR, dat$lowerSR, dat$upperSR, dat$spp) # bind sex and species

######### subset the dat that need to be analysed only
## subset 
dat_BC = subset(dat, spp=="Black coucal") 
head(dat_BC)
str(dat_BC)

## subset of the data for postfledging survival
dat_WBC = subset(dat, spp=="White-browed coucal") 
head(dat_WBC)
str(dat_WBC)

######### Plotting SR with the upper and lower 95% CI

par( mfrow = c( 1, 2 ))

## plot for BC 
plot(meanSR ~ age, 
     data = dat_BC, 
     type = "l",
     col="black",
     lwd=5,
     ylim = c(0, 1),   # c(0.45, 0.55), 
     ylab = "sex ratio (proportion of males)",
     xlim =c(0, 93),
     xlab=" Age (days)",
     main= " ", # "mean sex ratios of BC nestlings before leaving the nest",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
 lines(dat_BC$age, dat_BC$lowerSR, col = "black", lty = 2, lwd=2) ##
 lines(dat_BC$age, dat_BC$upperSR, col = "black", lty = 2, lwd=2)
abline(h=0.5, col="green", lty = 1, lwd=3) # 

# Add a legend
legend("topright", legend=c("mean", "95% CI", "Parity"),
       col=c("black", "black", "green"), lty=1:2, lwd= 3, cex=1.5)


# ---------------------------
#### separate plot for WBC
plot(meanSR ~ age, 
     data = dat_WBC, 
     type = "l",
     col="blue",
     lwd=5,
     ylim = c(0, 1),   # c(0.45, 0.55), 
     ylab = "sex ratio (proportion of males)",
     xlim =c(0, 20),
     xlab=" Age (days)",
     main= " ", # "mean sex ratios of BC nestlings before leaving the nest",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
lines(dat_WBC$age, dat_WBC$lowerSR, col = "blue", lty = 2, lwd=2) ##
lines(dat_WBC$age, dat_WBC$upperSR, col = "blue", lty = 2, lwd=2)
abline(h=0.5, col="green", lty = 1, lwd=3) # 

# Add a legend
legend("topright", legend=c("mean", "95% CI", "Parity"),
       col=c("blue", "blue", "green"), lty=1:2, lwd= 3, cex=1.5)

#---------------------------------------
### combined pre-fledging SR
plot(meanSR ~ age, 
     data = dat_BC, 
     type = "l",
     col="black",
     lwd=5,
     ylim = c(0, 1),   # c(0.45, 0.55), 
     ylab = "sex ratio (proportion of males)",
     xlim =c(0, 20),
     xlab=" Age (days)",
     main="mean sex ratios of BC & WBC nestlings before leaving the nest",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
lines(dat_BC$age, dat_BC$lowerSR, col = "black", lty = 2, lwd=2) ##
lines(dat_BC$age, dat_BC$upperSR, col = "black", lty = 2, lwd=2)
abline(h=0.5, col="green", lty = 1, lwd=3) # 

### plot WBC data on the same plot with BC

lines(dat_WBC$age, dat_WBC$meanSR, col = "blue", lty = 1, lwd=5) ##
lines(dat_WBC$age, dat_WBC$upperSR, col = "blue", lty = 2, lwd=1.5)
lines(dat_WBC$age, dat_WBC$lowerSR, col = "blue", lty = 2, lwd=1.5) ##

# Add a legend
legend("topright", legend=c("Black coucal", "White-browed coucal"),
       col=c("black", "blue"), lty=1:1, lwd= 4, cex=2)

## ----------------------------------------------------------------
# plotting postfledge SR 








##
# -------------------------
# calculate binary means & 95% conf intervals
##
list.files()
dat2 = read.csv("Pref_ageWBC.csv", header=T)
head(dat2)
str(dat2)
# Get CI from binomial distribution
wbc_CI <- binom.exact(dat2$males, dat2$total)
wbc_CI

# ----
dat3 = read.csv("Pref_ageBC.csv", header=T)
head(dat3)
str(dat3)
# Get CI from binomial distribution
bc_CI <- binom.exact(dat3$males, dat3$total)
bc_CI
#------------------------------------------------ 















## ------------------------------------------------
## Analysis of sex specific postfledging survival
## ------------------------------------------------
### modelling postfledging survival
# BC postfledge
cox.BC.postf <- coxme(Surv(postf_age, postf_status==1) ~ sex + hatch_order + brood_size + sex*hatch_order + (1|year) + (1|nest_ID), data =  dat_BC_postf)
summary(cox.BC.postf)


# BC postfledge
cox.WBC.postf <- coxme(Surv(postf_age, postf_status==1) ~ sex + hatch_order + brood_size + sex*hatch_order + fledging_age + fledging_weight + (1|year) + (1|nest_ID), data =  dat_WBC_postf)
summary(cox.WBC.postf)

# -----------------------------------------------
### create postfleding survival object by sexes in BC & WBC
###  BC prefedge
postf_bcfit.by.sex <- survfit(Surv(postf_age, postf_status) ~ sex , data = dat_BC_postf) 
postf_bcfit.by.sex

prob_postfBC <- summary(postf_bcfit.by.sex, times = c(0:130)) # survival probabilities for males & females from day 1 until day 130
prob_postfBC

pref_bcfit.by.sex <- survfit(Surv(pref_age, pref_status) ~ sex + hatch_order, data = dat_BC_pref) 
pref_bcfit.by.sex

summary(pref_bcfit.by.sex) ### gives summaries of events, number at risk, survival estimates & confidence intervals

pref_surv_diff_BC_sexes <- survdiff(Surv(pref_age, pref_status) ~ sex + hatch_order, data = dat_BC_pref) # use Log-Rank test to compare survival for males and females BC
pref_surv_diff_BC_sexes

pref_surv_diff_BC_sexes_strataHO  <- survdiff(Surv(pref_age, pref_status) ~ sex + strata(hatch_order), data = dat_BC_pref) # use Log-Rank test to compare survival for males and females BC statified by hatching order
pref_surv_diff_BC_sexes_strataHO

##-----------------------------------------------


