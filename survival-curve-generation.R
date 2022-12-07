library(flexsurv)
library(survival)
library(readxl)
library(tidyr)

#This code assumes the user has already digitized a Kaplan-Meier curve and has generated sample #independent patient data (IPD) from the digitized curve
#Replace “directory pathway” with user’s PC directory containing the IPD file 
setwd(“directory pathway”) 

#data <- read.csv("pt_level_data.csv")           #Import the patient-level datapoints 
filename <- "pt_level_data.csv"                   

#Load data, remove NA rows and replace all 10000 with Inf
filename <- read_excel(filename, sheet="R data") %>% 
  drop_na() %>% 
  replace(((.)==10000), Inf)

#Create vector of start and end times, total number of persons in the study 

#The second parameter specified in the rep() function instructs R how many times it should replicate the first parameter specified 
times_start <- c(rep(filename$start_time_censor, filename$n_censors), rep(filename$start_time_event, filename$n_events))
times_end <- c(rep(filename$end_time_censor, filename$n_censors), rep(filename$end_time_event, filename$n_events))

n <- sum(n_events) + sum(n_censors)

#Add times for patients at risk at the last time point 
#repl1 <- unname(unlist(read_excel(filename, sheet="Number events & censored", range="O16", col_names = FALSE)))
#repl2 <- unname(unlist(read_excel(filename, sheet="Number events & censored", range="O15", col_names=FALSE)))

#Manually enter these numbers from the excel sheet
repl1 <- 30
repl2 <- 60

times_start <- c(times_start, rep(repl1, repl2))
times_end <- c(times_end, rep(Inf,repl2))

#Run Models 
#The Surv function must return the response variable 
#interval2 specifies interval-censored observations 
exp <- survreg(Surv(times_start, times_end, type="interval2")~1, dist="exponential")
weib <- survreg(Surv(times_start, times_end, type="interval2")~1, dist="weibull")
lnorm <- survreg(Surv(times_start, times_end, type="interval2")~1, dist="lognormal")
llogis <- survreg(Surv(times_start, times_end, type="interval2")~1, dist="loglogistic")


gamma <- flexsurvreg(Surv(times_start, times_end, type="interval2")~1, dist="gamma")
gomp <- flexsurvreg(Surv(times_start, times_end, type="interval2")~1, dist="gompertz")
gengamma1 <- flexsurvreg(Surv(times_start, times_end, type="interval2")~1, dist="gengamma")


#Plot the Curves 
survplot <- survfit(Surv(times_start, times_end, type="interval2")~1)
plot(survplot, conf.int=TRUE) #xlim=c(0,40))


#Add predicted curves for surv objects 
lines(gomp, col="brown", ci=FALSE)
lines(gamma,col="yellow",ci=FALSE)
lines(gengamma1,col="violet",ci=TRUE)

lines(predict(exp, type="quantile", p=seq(.01,.99, by=.01))[1,], seq(.99,.01,by=-.01),col="red", type="l")
lines(predict(weib, type="quantile", p=seq(.01,.99, by=.01))[1,], seq(.99,.01,by=-.01),col="blue", type="l")
lines(predict(lnorm, type="quantile", p=seq(.01,.99, by=.01))[1,], seq(.99,.01,by=-.01),col="green", type="l")
lines(predict(llogis, type="quantile", p=seq(.01,.99, by=.01))[1,], seq(.99,.01,by=-.01),col="cyan", type="l")

#Manual Construction of Survival Curves (mimics excel): 
survplot <- survfit(Surv(times_start, times_end, type="interval2")~1)
plot(survplot,conf.int = FALSE, xlim=c(0,60))

#Exponential 
summary(exp)$table
t(chol(vcov(exp)))
sqrt(diag(vcov(exp)))
lambda <- exp(-exp$icoef)
lines(x=seq(0,100,by=1),y=exp(-lambda*seq(0,100,by=1)), type="l")
logLik(exp)

#Weibull 
summary(weib)$table
t(chol(vcov(weib)))
c <- 1/exp(weib$icoef[2])
lines(x=seq(0,100,by=1),y=exp(-(seq(0,100,by=1)/b)**c), type="l")
logLik(weib)

# lognormal
summary(lnorm)$table
t(chol(vcov(lnorm)))
a <- lnorm$icoef[1]
b <- exp(lnorm$icoef[2])
lines(x=seq(0,100,by=1),y=1-pnorm((log(seq(0,100,by=1))-a)/b), type="l")
logLik(lnorm)

# loglogistic
summary(llogis)$table
t(chol(vcov(llogis)))
b <- exp(llogis$icoef[1])
c <- 1/exp(llogis$icoef[2])
lines(x=seq(0,100,by=1),y=(1+((seq(0,100,by=1)/b)**c))**-1, type="l")
logLik(llogis)

# gamma
gamma$res.t
t(chol(vcov(gamma)))
shape <- exp(gamma$coefficients[1])
rate <- exp(gamma$coefficients[2])
lines(x=seq(0,100,by=1),y=1-pgamma(seq(0,100,by=1),shape,rate), type="l")
logLik(gamma)


# generalized gamma 1 (three parameters)
summary(gengamma1)
gengamma1$res.t
t(chol(vcov(gengamma1)))
mu <- gengamma1$coefficients["mu"]
sigma <- exp(gengamma1$coefficients["sigma"])
Q <- gengamma1$coefficients["Q"]
lines(x=seq(0,100,by=0.1),y=1-pgengamma(seq(0,100,by=0.1),mu,sigma,Q), type="l")
logLik(gengamma1)

# gompertz (uses flexsurv's rgompertz parameterization, a=shape, b=rate)
gomp$res.t
t(chol(vcov(gomp)))
a <- gomp$coefficients[1]
b <- exp(gomp$coefficients[2])
lines(x=seq(0,100,by=1),y=exp((-b/a)*(exp(a*seq(0,100,by=1))-1)), type="l",col="green")
logLik(gomp)

#Plot to fit the original kaplan-meier curve 
km <- survfit(Surv(times_start, times_end, type="interval2")~1)
summary(km)
plot(km, xlim=c(0,50))

lines(x=seq(0,100,by=0.1),y=1-pgengamma(seq(0,100,by=0.1),mu,sigma,Q), type="l")

