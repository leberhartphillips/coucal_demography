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
  select(species, ring_ID, lab_no, sex, year, site, nest_ID, pref_age, 
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

BC_nest_dat <- 
  status_dat_nestlings %>% 
  filter(species == "BC")

BC_fled_dat <- 
  status_dat_fledglings %>% 
  filter(species == "BC")

# hazard functions by treatment period
# note Surv() in the sshzd() function takes the arguments in a different order than Surv() in the survival package
# hazard function = hz
hz.M.f <- sshzd(Surv(exit, event, entry)~exit, data = filter(BC_fled_dat, sex == "M"), alpha = 1.4) #alpha=0.25 will overfit
hz.F.f <- sshzd(Surv(exit, event, entry)~exit, data = filter(BC_fled_dat, sex == "F"), alpha = 1.4)

hz.M.n <- sshzd(Surv(exit, event, entry)~exit, data = filter(BC_nest_dat, sex == "M"), alpha = 1.4) #alpha=0.25 will overfit
hz.F.n <- sshzd(Surv(exit, event, entry)~exit, data = filter(BC_nest_dat, sex == "F"), alpha = 1.4)

# calculate the hazard rate function for plotting for a sequence of time (tt = 1 to 52 weeks by 0.5)
# hazard = hz
tt.n <- seq(0, 15, 1)
tt.f <- seq(16, 70, 1)

# predicted hazard
# Female nestlings
est.M.n <- hzdrate.sshzd(hz.M.n, tt.n, se = TRUE)
hz.M.n.fit = 1 - est.M.n$fit
hz.M.n.hi <- 1 - est.M.n$fit*exp(1.96*est.M.n$se)
hz.M.n.lo <- 1 - est.M.n$fit/exp(1.96*est.M.n$se)

# Male nestlings
est.F.n <- hzdrate.sshzd(hz.F.n, tt.n, se = TRUE)
hz.F.n.fit = 1 - est.F.n$fit
hz.F.n.hi <- 1 - est.F.n$fit*exp(1.96*est.F.n$se)
hz.F.n.lo <- 1 - est.F.n$fit/exp(1.96*est.F.n$se)

# predicted hazard
# Female fledglings
est.M.f <- hzdrate.sshzd(hz.M.f, tt.f, se = TRUE)
hz.M.f.fit = 1 - est.M.f$fit
hz.M.f.hi <- 1 - est.M.f$fit*exp(1.96*est.M.f$se)
hz.M.f.lo <- 1 - est.M.f$fit/exp(1.96*est.M.f$se)

# Male fledglings
est.F.f <- hzdrate.sshzd(hz.F.f, tt.f, se = TRUE)
hz.F.f.fit = 1 - est.F.f$fit
hz.F.f.hi <- 1 - est.F.f$fit*exp(1.96*est.F.f$se)
hz.F.f.lo <- 1 - est.F.f$fit/exp(1.96*est.F.f$se)

# graph for hazard functions in males and females periods
plot(tt.n, hz.M.n.fit, type="l", lty=1, lwd=2, col="blue", xaxs="i", yaxs="i", 
     xlim=c(0, 70), ylim=c(0.7, 1), 
     xlab=c("days since hatching"), ylab=c("hazard function"))
lines(tt.n, hz.M.n.hi, lty=2, lwd=2, col="blue")
lines(tt.n, hz.M.n.lo, lty=2, lwd=2, col="blue")
lines(tt.n, hz.F.n.fit, lty=1, lwd=2, col="red")
lines(tt.n, hz.F.n.hi, lty=2, lwd=2, col="red")
lines(tt.n, hz.F.n.lo, lty=2, lwd=2, col="red")

lines(tt.f, hz.M.f.fit, lty=1, lwd=2, col="blue")
lines(tt.f, hz.M.f.hi, lty=2, lwd=2, col="blue")
lines(tt.f, hz.M.f.lo, lty=2, lwd=2, col="blue")
lines(tt.f, hz.F.f.fit, lty=1, lwd=2, col="red")
lines(tt.f, hz.F.f.hi, lty=2, lwd=2, col="red")
lines(tt.f, hz.F.f.lo, lty=2, lwd=2, col="red")

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