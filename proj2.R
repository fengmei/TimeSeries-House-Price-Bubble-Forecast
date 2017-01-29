install.packages("tseries")
install.packages("TTR")
install.packages("forecast")
library("forecast")
require("tseries") 
## Read data2
data2 <-scan(file="~/SJSU_ClASSES/Math265_TimeSeries/FinalProject/proj2.txt")

## explore data
head(data2)
summary(data2)
sapply(data2,class)

#names(data1) <- c("Proj1")
length(data2)

## plot data
plot.ts(data2,type="b",main="PROJ2 Time Series", xlab="Time",ylab="Value")
plot.ts(data2[1000:1500],type="b")  ## looks stationary with some peridicity
# stationary test: pass
adf.test(data2)

## acf, pacf
par(mfrow=c(1,2))
   acf(data2, type = c("correlation"),lag.max=50,main="Plot of Proj2 Sample Autocorrelation")
   acf(data2, type = c("partial"),lag.max=50,main="Plot of Proj2 Sample Partial Autocorrelation")
par(mfrow=c(1,1))  ## acf sinusodical decay to 0, pacf tells that this maybe a ARMA process

## Periodogram
spec.pgram(data2,taper=0.1, main="Proj2 Periodogram") ## maybe ARMA(2,2),ARMA(2,3),ARMA(3,2), ARMA(3,3) or try larger one to ARMA(4,4).....

##
f <- auto.arima(data2,max.p=5,max.q=6,start.p=2,start.q=0,stationary=TRUE,seasonal=FALSE)
f
## difference? Doesn't look better.

differenced.data2 <- diff(data2, differences = 1)
plot(differenced.data2, type="b")

par(mfrow=c(2,1))
   acf(differenced.data2, type = c("correlation"))
   acf(differenced.data2, type = c("partial"))
par(mfrow=c(1,1))  

## log? not better
log.data2 <- log(data2)
plot(log.data2,type="b")

## check and subtract the mean
##data1_zeromean <- data1 - mean(data1)
mean(data2)

## so candidates models are ARMA(2,2),ARMA(2,3), ARMA(3,2),ARMA(3,3), ARMA(4,2),ARMA(2,4), ARMA(4,3),ARMA(3,4) fit the candidates model

## Split data to 90% for training and 10% for test 

fitting.data2 <- data2[1:1400]
test.data2 <- data2[1401:1500]


## Fit the models

data2_fit1 <- arima(fitting.data2,order=c(2,0,2),method="ML")
data2_fit2 <- arima(fitting.data2,order=c(2,0,3),method="ML")
data2_fit3 <- arima(fitting.data2,order=c(3,0,2),method="ML")
data2_fit4 <- arima(fitting.data2,order=c(3,0,3),method="ML")
data2_fit5 <- arima(fitting.data2,order=c(4,0,2),method="ML")
data2_fit6 <- arima(fitting.data2,order=c(2,0,4),method="ML")
data2_fit7 <- arima(fitting.data2,order=c(4,0,3),method="ML")
data2_fit8 <- arima(fitting.data2,order=c(3,0,4),method="ML")

data2_fit1
data2_fit2
data2_fit3
data2_fit4
data2_fit5
data2_fit6
data2_fit7
data2_fit8


names(data2_fit1)
### check MSE

sum(((data2_fit1$res)^2))/(length(fitting.data2) -(2+2) -(2+2+1))
sum(((data2_fit2$res)^2))/(length(fitting.data2) -(2+3) -(2+3+1))
sum(((data2_fit3$res)^2))/(length(fitting.data2) -(3+2) -(3+2+1))
sum(((data2_fit4$res)^2))/(length(fitting.data2) -(3+3) -(3+3+1))
sum(((data2_fit5$res)^2))/(length(fitting.data2) -(4+2) -(4+2+1))
sum(((data2_fit6$res)^2))/(length(fitting.data2) -(2+4) -(2+4+1))
sum(((data2_fit7$res)^2))/(length(fitting.data2) -(4+3) -(4+3+1))
sum(((data2_fit8$res)^2))/(length(fitting.data2) -(3+4) -(3+4+1))

## coefficient significance test show that ARMA(3,2),ARMA(4,2),ARMA(4,3),ARMA(3,4) are not significant
## Choose In the rest models ARMA(2,2), ARMA(2,3), ARMA(3,3), ARMA(2,4)


## diagnostics

tsdiag(data2_fit1)
tsdiag(data2_fit2)
tsdiag(data2_fit3)
tsdiag(data2_fit4)
tsdiag(data2_fit5)
tsdiag(data2_fit6)
tsdiag(data2_fit7)
tsdiag(data2_fit8)


qqnorm(data2_fit4$resi, main="ARMA(3,3) Normal QQ plot")


## plot spectrum against the periodogram

my.spectrum <- function(phi.of.b, theta.of.b, variance=1)
{
   p <- length(phi.of.b)
   q <- length(theta.of.b)

   omega <- seq(from=0, to=pi, by=.001)

   phi.of.e.minus.i.omega <- 1
   phi.of.e.i.omega <- 1
  
   if(p>1)
   {   for(i in 2:p)
      {
          phi.of.e.minus.i.omega <-  phi.of.e.minus.i.omega + phi.of.b[i]*exp(complex(imaginary = -(i-1))*omega)
          phi.of.e.i.omega <-  phi.of.e.i.omega + phi.of.b[i]*exp(complex(imaginary = (i-1))*omega)
      }
   }

   theta.of.e.minus.i.omega <- 1
   theta.of.e.i.omega <- 1

   if(q>1)
   {
      for(i in 2:q)
      {
          theta.of.e.minus.i.omega <-  theta.of.e.minus.i.omega + theta.of.b[i]*exp(complex(imaginary = -(i-1))*omega)
          theta.of.e.i.omega <-  theta.of.e.i.omega + theta.of.b[i]*exp(complex(imaginary = (i-1))*omega)
      }
   }

    my.spectrum <- (variance/(2*pi))*Re(theta.of.e.minus.i.omega*theta.of.e.i.omega/(phi.of.e.minus.i.omega*phi.of.e.i.omega))

   plot(omega, 10*log10(my.spectrum), ylab="spectrum (in decibels)", type="l")   
}


par(mfrow=c(2,1))
   spec.pgram(fitting.data2)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data2_fit3$coef[1:3])), theta.of.b=c(1,data2_fit3$coef[4:5]), variance=data2_fit3$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data2, main="Original Training Data Periodogram")  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data2_fit4$coef[1:3])), theta.of.b=c(1,data2_fit4$coef[4:6]), variance=data2_fit4$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data2)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data2_fit5$coef[1:4])), theta.of.b=c(1,data2_fit5$coef[5:6]), variance=data2_fit5$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data2)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data2_fit7$coef[1:4])), theta.of.b=c(1,data2_fit7$coef[5:7]), variance=data2_fit7$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data2)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data2_fit5$coef[1:4])), theta.of.b=c(1,data2_fit5$coef[5:6]), variance=data2_fit5$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data2)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data2_fit4$coef[1:3])), theta.of.b=c(1,data2_fit4$coef[4:6]), variance=data2_fit4$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data2)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data2_fit4$coef[1:3])), theta.of.b=c(1,data2_fit4$coef[4:6]), variance=data2_fit4$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data2)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data2_fit4$coef[1:3])), theta.of.b=c(1,data2_fit4$coef[4:6]), variance=data2_fit4$sigma2)
par(mfrow=c(1,1))


## plot theoretical acf against the sample acf

par(mfrow=c(2,2))
acf(fitting.data2, type = c("correlation"))
acf(fitting.data2, type = c("partial"))
plot(ARMAacf(ar=c(data2_fit4$coef[1:5]),lag.max=25),type="h")
abline(h=0)
plot(ARMAacf(ar=c(data2_fit3$coef[1:5]),lag.max=25, pacf=TRUE),type="h")
abline(h=0)

par(mfrow=c(2,2))
acf(fitting.data1, type = c("correlation"))
acf(fitting.data1, type = c("partial"))
plot(ARMAacf(ar=c(data1_fit7$coef[1:9]),lag.max=25),type="h")
abline(h=0)
plot(ARMAacf(ar=c(data1_fit7$coef[1:9]),lag.max=25, pacf=TRUE),type="h")
abline(h=0)

par(mfrow=c(2,2))
acf(fitting.data2, type = c("correlation"),main="Original training data ACF")
acf(fitting.data2, type = c("partial"),main="Original training data pACF")
plot(ARMAacf(ar=c(data2_fit4$coef[1:3]),ma=c(data2_fit4$coef[4:6]),lag.max=25),type="h",ylab="acf",main="Model ARMA(3,3) acf")
abline(h=0)
plot(ARMAacf(ar=c(data2_fit4$coef[1:3]),ma=c(data2_fit4$coef[4:6]),lag.max=25, pacf=TRUE),type="h",ylab="pacf",main="Model ARMA(3,3) pacf")
abline(h=0)

par(mfrow=c(1,1))

par(mfrow=c(2,2))
acf(fitting.data2, type = c("correlation"),main="Original training data ACF")
acf(fitting.data2, type = c("partial"),main="Original training data pACF")
plot(ARMAacf(ar=c(data2_fit5$coef[1:4]),ma=c(data2_fit5$coef[5:6]),lag.max=25),type="h",ylab="acf",main="Model ARMA(3,3) acf")
abline(h=0)
plot(ARMAacf(ar=c(data2_fit5$coef[1:4]),ma=c(data2_fit5$coef[5:6]),lag.max=25, pacf=TRUE),type="h",ylab="pacf",main="Model ARMA(3,3) pacf")
abline(h=0)

par(mfrow=c(1,1))

## predict using Test data 

data2_pred1<- predict(data2_fit1, n.ahead=100, se.fit=FALSE)
data2_pred2<- predict(data2_fit2, n.ahead=100, se.fit=FALSE)
data2_pred3<- predict(data2_fit3, n.ahead=100, se.fit=FALSE)
data2_pred4<- predict(data2_fit4, n.ahead=100, se.fit=FALSE)
data2_pred5<- predict(data2_fit5, n.ahead=100, se.fit=FALSE)
data2_pred6<- predict(data2_fit6, n.ahead=100, se.fit=FALSE)
data2_pred7<- predict(data2_fit7, n.ahead=100, se.fit=FALSE)
data2_pred8<- predict(data2_fit8, n.ahead=100, se.fit=FALSE)

data2_sse1 <- sum((test.data2 - data2_pred1)^2)
data2_sse2 <- sum((test.data2 - data2_pred2)^2)
data2_sse3 <- sum((test.data2 - data2_pred3)^2)
data2_sse4 <- sum((test.data2 - data2_pred4)^2)
data2_sse5 <- sum((test.data2 - data2_pred5)^2)
data2_sse6 <- sum((test.data2 - data2_pred6)^2)
data2_sse7 <- sum((test.data2 - data2_pred7)^2)
data2_sse8 <- sum((test.data2 - data2_pred8)^2)

data2_sse1
data2_sse2
data2_sse3
data2_sse4
data2_sse5
data2_sse6
data2_sse7
data2_sse8


## Plot and see the prediction
par(mfrow=c(2,1))
plot(data2_pred5,type="b")
plot(test.data2,type="b")
par(mfrow=c(1,1))

## COMPARE Models         AR(2)   AR(3)   AR(4)  AR(5)  AR(6)   AR(7)   AR(8)    AR(9)    ARMA(2,1)  ARMA(3,1) 
## Coef Signif Test	   |   Y		Y		Y	   Y	  Y		  Y		  Y		  Y          Y		    Y			
##  AIC                |5075.41 4924.98 4850.83 4824.11 4802.03 4796.08 4783.13  4782.34   5009.06    4758.6 
## MSE				   |4572.86 3266.52 2766.15 2605.04 2478.72 2445.71 2375.2  2371.08   3942.97    2243.53 
## Periodogram         |   L		L		H	   H	  M		  M		  M		  M			 L          H					
## test prediction sse |        3200350  3061494 3014249 3091009 2999213 3127040 3047744   3312724    2975883
## Parsinomy		   |   Y		Y		Y	   Y	  Y		  N		  N		  N			  Y 		Y				
## Diagnostics		   |   N		N		N	   M	  M		  Y       Y       Y           N         Y

## Fit the model to the entire dataset
data2.fit <- arima(data2, order=c(3,0,3),method="ML")
sum(((data2.fit$res)^2))/(length(data2) -(3+3) -(3+3+1))

tsdiag(data2.fit)
qqnorm(data2.fit$resi)

## predict 13 points ahead
data2.pred <- predict(data2.fit, n.ahead=13, se.fit=TRUE)


data2.prediction <- data2.pred$pred
data2.pred.lowerbound <- data2.prediction - 2*data2.pred$se
data2.pred.upperbound <- data2.prediction + 2*data2.pred$se

cbind(data2.pred.lowerbound,data2.prediction,data2.pred.upperbound)
save <- rbind(data2.pred.lowerbound,data2.prediction,data2.pred.upperbound)
write.table(save, file="proj2.csv",sep=",")

## plot
plot(1400:1500, data2[1400:1500], xlim=c(1400, 1515),  type="b",xlab="Time", ylab="Value", main="Proj 2 Forecast of Future 13 Points")

lines(1501:1513, data2.prediction, type="b", col=2)
lines(1501:1513, data2.pred.lowerbound, type="l", col=2)
lines(1501:1513, data2.pred.upperbound, type="l", col=2)
plot(forecast(data2.fit,h=13),xlim=c(1400,1515))



