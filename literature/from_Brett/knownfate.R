# knownfate.R
# time to event analysis for known fate data
# survival of female prairie chickens pre/post wind power
# Winder et al. 2014 Journal of Applied Ecology 51:395-405
# Updated 23 August 2019

# bookkeeping
rm(list = ls())

# INITIALIZE PACKAGES

# install packages if needed
# install.packages("survival")
# install.packages("gss")

# load package survival, needed for Kaplan-Meier analyses
library(survival)
# load package gss, needed for hazard functions
library(gss)

# DATA PROCESSING

# get datafile enchist.txt from folder
#pc <- read.table(file.choose(),header=T,sep=",")
# setwd("C:/Users/brett.sandercock/OneDrive - NINA/Documents/Teaching/Rworkshop/Mod3survival")
pc <- read.table("literature/from_Brett/enchist.txt",header=T,sep=",")
pc.pre <- subset(pc, trt=="pre") # subset data
pc.post <- subset(pc, trt=="post")

# look at data
# name = bird band number plus radio frequency, year = year of study (1-5)
# trt = pre/post-construction, entry and exit = week of study period
# event = 1 for died, 0 for survived, dturb = distance of home range to turbine in km
# head(pc, 20) # first 20 lines of file

# Kaplan-Meier comparisons
# trt comparison between pre and post, pooling across years
# mfit <- survfit(Surv(pc$entry, pc$exit, pc$event) ~ 1 )
mfit <- survfit(Surv(pc$entry, pc$exit, pc$event) ~ pc$trt )
print(survfit(Surv(pc$entry, pc$exit, pc$event) ~ pc$trt), print.rmean=TRUE)
summary(mfit)

# Plot the Kaplan-Meier cumulative survival
# xaxs, yaxs=i gets rid of rabbit holes at corner of graphs
plot(mfit, conf.int=T, lty=c(1,2,2), col=c("dark green", "red"), xaxs="i", yaxs="i", xlim=c(0, 53), lwd=2.5, xlab=c("week"), ylab=c("survival function"))
legend("topright", legend = c("Post", "Pre"), cex=2, lty=1, col=c("dark green", "red"), lwd=2.5, bty="n")

# Cox proportional hazards
# test of interactive model for treatment period, distance to turbine
# cluster is a random effects for individual bird, robust se in output are adjusted values
trt.prepost <- coxph(Surv(pc$entry, pc$exit, pc$event)~ pc$trt + pc$dturb + cluster(pc$name))	 # Interactive model 
summary(trt.prepost)		    # Summary of the model			
cox.zph(trt.prepost)              # Model diagnostic; testing the proportional hazard assumption (see "?cox.zph")		
windows()
plot(cox.zph(trt.prepost, transform="identity"))   # Plotting the model diagnostics, identity puts time on same scale as KM plots

# test of main effects models for treatment period, distance to turbine, cluster is a random effects for individual bird
trt.prepost <- coxph(Surv(pc$entry, pc$exit, pc$event)~ pc$trt + pc$dturb + cluster(pc$name))	 # Additive model 
summary(trt.prepost)		    # Summary of the model			

# Pooling across both periods
#mfit <- survfit(Surv(pc$entry, pc$exit, pc$event) ~ 1)
#print(survfit(Surv(pc$entry, pc$exit, pc$event) ~ 1), print.rmean=TRUE)
#summary(mfit)
# xaxs, yaxs=i gets rid of rabbit holes at corner of graphs
#plot(mfit, lty=1, col=c("dark green"), xaxs="i", yaxs="i", xlim=c(0, 53), lwd=2.5, xlab=c("week"), ylab=c("survival function"))
#legend("topright", legend = c("All years"), cex=2, lty=1, col=c("dark green"), lwd=2.5, bty="n")

# HAZARD FUNCTION NONCYCLIC

# smoothing factor, default of 1, higher values are more smoothed
# Gu 2014 J. Stat Software 58:art5 recommends alpha=1.4 
alpha.smoothie = 1.2

# hazard functions by treatment period
# note Surv() in the sshzd() function takes the arguments in a different order than Surv() in the survival package
# hazard function = hz
hz.pre <- sshzd(Surv(exit, event, entry)~exit, data=pc.pre, alpha=1.2) #alpha=0.25 will overfit
hz.post <- sshzd(Surv(exit, event, entry)~exit, data=pc.post, alpha=1.2)
hz.all <- sshzd(Surv(exit, event, entry)~exit, data=pc, alpha=1.2)

# calculate the hazard rate function for plotting for a sequence of time (tt = 1 to 52 weeks by 0.5)
# hazard = hz
tt <- seq(1,52,0.5)

# predicted hazard
# pre-construction
est <- hzdrate.sshzd(hz.pre,tt,se=TRUE)
hzpre.fit = est$fit
hzpre.hi <- est$fit*exp(1.96*est$se)
hzpre.lo <- est$fit/exp(1.96*est$se)

# post-construction
est <- hzdrate.sshzd(hz.post,tt,se=TRUE)
hzpost.fit = est$fit
hzpost.hi <- est$fit*exp(1.96*est$se)
hzpost.lo <- est$fit/exp(1.96*est$se)

# pooled
est <- hzdrate.sshzd(hz.all,tt,se=TRUE)
hzall.fit = est$fit
hzall.hi <- est$fit*exp(1.96*est$se)
hzall.lo <- est$fit/exp(1.96*est$se)

# graph for hazard functions in pre and post periods
windows()  # new graphics window				
plot(tt, hzpre.fit, type="l", lty=1, lwd=2, col="red", xaxs="i", yaxs="i", xlim=c(0, 53), ylim=c(0, 0.12), xlab=c("week"), ylab=c("hazard function"))
lines(tt, hzpre.hi, lty=2, lwd=2, col="red")
lines(tt, hzpre.lo, lty=2, lwd=2, col="red")
lines(tt, hzpost.fit, lty=1, lwd=2, col="dark green")
lines(tt, hzpost.hi, lty=2, lwd=2, col="dark green")
lines(tt, hzpost.lo, lty=2, lwd=2, col="dark green")
legend("topright", legend = c("Post", "Pre"), cex=2, lty=c(1,1), col=c("dark green", "red"), lwd=2.5, bty="n")
axis(1, labels=F, tick=T)
# axis(2, at=seq(0,0.1,0.01), labels=T, tick=T)  # tweak for GKS3, change ylim in plot too
axis(2, labels=T, tick=T)

# graph for hazard functions pooling years
#windows()  # new graphics window				
#plot(tt, hh3, type="l", lwd=2, lty=2, col="dark green", xaxs="i", yaxs="i", xlim=c(0, 53), ylim=c(0, 0.04), xlab=c("week"), ylab=c("hazard function"))
#legend("topleft", legend = c("All years"), cex=2, lty=1:2, col=c("dark green", "red"), lwd=2.5, bty="n")
#axis(1, labels=F, tick=T)
#axis(2, labels=T, tick=T)