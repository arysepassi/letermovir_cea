#This code snippet assumes the user has generated a list of times (months) and survival probabilities from a piecewise curve (e.g. combination of data from Saullo et al and US life table data)

#Note: This code was used twice – once to generate a smoothed spline curve for letermovir survival data, and again for PET-only survival data. 

#Replace “directory pathway” with user’s directory for use: 
Setwd("directory pathway")

#Import file containing the survival probabilities and times: 
data<-read.csv("test-cubic-spline.csv")

library(splines)  
library(ISLR)

attach(data)

monthlims <- range(MONTH)
month.grid <- seq(from=monthlims[1], to=monthlims[2])

#Fit a cubic spline model with three knots (cutpoints) at three percentiles to test for fit: 
fit <- lm(P.SURVIVAL.PET~bs(MONTH, knots=c(30,50,65)),data=data)
summary(fit)

plot(MONTH, data$P.SURVIVAL.PET,col="grey",xlab="Month",ylab="Survival")
points(month.grid,predict(fit,newdata=list(MONTH=month.grid)),col="darkgreen",lwd=2, type="l")

#Test using smoothing splines: 
fit1 <- smooth.spline(MONTH,P.SURVIVAL.PET,df=25) #25 degrees of freedom for splines fit best 

plot(MONTH,P.SURVIVAL.PET,col="grey", xlab="Months", ylab="P(Survival)")
points(month.grid,predict(fit,newdata=list(MONTH=month.grid)),col="darkgreen",lwd=2,type='l')
abline(v=c(30,50,65),lty=2,col="darkgreen")
lines(fit1,col="red",lwd=2)
legend("topright",c("Smoothing Spline with 25 df","Cubic Spline"),col=c("red","darkgreen"),lwd=2)


fit2 <- smooth.spline(MONTH, P.SURVIVAL.PET,cv=TRUE)
fit2

plot(MONTH,P.SURVIVAL.PET,col="grey")
lines(fit2, lwd=2,col='purple')
legend("topright",("Smoothing Splines with 102.7 df selected by CV"),col="purple",lwd=2)

#At this point, the user can select which model provided the best fit for the original data and export the data for use in the model: 
#Extract survival data: 
spline <- data.frame(fit2$x,
                     fit2$y)

write.csv(spline, "write to user's path")
