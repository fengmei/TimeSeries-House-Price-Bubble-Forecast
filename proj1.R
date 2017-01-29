install.packages("tseries")
install.packages("TTR")
install.packages("forecast")
library("forecast")
require("tseries") 
## Read data1
data1 <-scan(file="~/SJSU_ClASSES/Math265_TimeSeries/FinalProject/proj1.txt")

## explore data
head(data1)
summary(data1)
sapply(data1,class)

#names(data1) <- c("Proj1")
length(data1)

## plot data
plot.ts(data1,type="b",main="Data1 Time Series", xlab="Time",ylab="Value")  ## looks stationary

## acf, pacf
par(mfrow=c(1,2))
   acf(data1, type = c("correlation"), main="Plot of Proj1 Sample Autocorrelation")
   acf(data1, type = c("partial"),main="Plot of Proj1 Sample Partial Autocorrelation")
par(mfrow=c(1,1))  ## maybe AR(4), AR(5), AR(6), AR(7), AR(8)

## Periodogram
spec.pgram(data1,taper=0.1, main="Proj 1 Periodogram") ## maybe AR(2), AR(3), AR(4), ARMA(3,1), ARMA(3,2),ARMA(4,2), ARMA(4,4)

## using auto.arima to test a model
##f <- auto.arima(data1,max.p=10,max.q=5,start.p=2,start.q=0,stationary=TRUE,seasonal=FALSE)
##f <- auto.arima(data1,max.p=10,max.q=5,start.p=2,start.q=1,stationary=TRUE,seasonal=FALSE)

## difference? Doesn't look better.
adf.test(data1)
differenced.data1 <- diff(data1, differences = 1)
plot(differenced.data1, type="b")

par(mfrow=c(2,1))
   acf(differenced.data1, type = c("correlation"))
   acf(differenced.data1, type = c("partial"))
par(mfrow=c(1,1))  

## log? non-stationary now.No!
log.data1 <- log(data1)
plot(log.data1,type="b")

## check and subtract the mean
##data1_zeromean <- data1 - mean(data1)
mean(data1)

## so candidates models are AR(3), AR(4), AR(5),AR(6), AR(7), AR(8),ARMA(3,2),ARMA(4,2), ARMA(4,4)fit the candidates model
#ar(data1,aic=FALSE, order.max=4, method=c("mle"))
#ar(data1,aic=FALSE, order.max=4, method=c("ols"))

## Split data to 90% for training and 10% for test 

fitting.data1 <- data1[1:450]
test.data1 <- data1[451:501]



data1_fit1_CSS <- arima(fitting.data1,order=c(3,0,0))

data1_fit0 <- arima(fitting.data1,order=c(2,0,0),method="ML")
data1_fit1 <- arima(fitting.data1,order=c(3,0,0),method="ML")
data1_fit2 <- arima(fitting.data1,order=c(4,0,0),method="ML")
data1_fit3 <- arima(fitting.data1,order=c(5,0,0),method="ML")
data1_fit4 <- arima(fitting.data1,order=c(6,0,0),method="ML")
data1_fit5 <- arima(fitting.data1,order=c(7,0,0),method="ML")
data1_fit6 <- arima(fitting.data1,order=c(8,0,0),method="ML")
data1_fit7 <- arima(fitting.data1,order=c(9,0,0),method="ML")
data1_fit8 <- arima(fitting.data1,order=c(2,0,1))
data1_fit9 <- arima(fitting.data1,order=c(3,0,1),method="ML")

data1_fit15 <- arima(fitting.data1,order=c(2,0,2))
data1_fit10 <- arima(fitting.data1,order=c(3,0,2),method="ML")
data1_fit11 <- arima(fitting.data1,order=c(10,0,0),method="ML")
data1_fit12 <- arima(fitting.data1,order=c(4,0,1),method="ML")
data1_fit13 <- arima(fitting.data1,order=c(4,0,2),method="ML")
data1_fit14 <- arima(fitting.data1,order=c(4,0,4),method="ML")

data1_fit1_CSS
data1_fit0
data1_fit1  ## CSS and ML got the very close results, so choose ML method
data1_fit2
data1_fit3
data1_fit4
data1_fit5
data1_fit6
data1_fit7
data1_fit8
data1_fit9
data1_fit10
data1_fit11
data1_fit12
data1_fit13
data1_fit14
data1_fit15

names(data1_fit1)
### ???? ## check MSE
sum(((data1_fit0$res)^2))/(length(fitting.data1) -(2) -(2+1))
sum(((data1_fit1$res)^2))/(length(fitting.data1) -(3) -(3+1))
sum(((data1_fit2$res)^2))/(length(fitting.data1) -(4) -(4+1))
sum(((data1_fit3$res)^2))/(length(fitting.data1) -(5) -(5+1))
sum(((data1_fit4$res)^2))/(length(fitting.data1) -(6) -(6+1))
sum(((data1_fit5$res)^2))/(length(fitting.data1) -(7) -(7+1))
sum(((data1_fit6$res)^2))/(length(fitting.data1) -(8) -(8+1))
sum(((data1_fit7$res)^2))/(length(fitting.data1) -(9) -(9+1))
sum(((data1_fit8$res)^2))/(length(fitting.data1) -(2+1) -(2+1+1))
sum(((data1_fit9$res)^2))/(length(fitting.data1) -(3+1) -(3+1+1))
#sum(((data1_fit10$res)^2))/(length(fitting.data1) -(4+2) -(4+2+1))
#sum(((data1_fit11$res)^2))/(length(fitting.data1) -(4+4) -(4+4+1))


## coefficient significance test show that ARMA(3,2),ARMA(4,2),ARMA(4,4) are not significant
## Choose In the rest models


## diagnostics
tsdiag(data1_fit0)
tsdiag(data1_fit1)
tsdiag(data1_fit2)
tsdiag(data1_fit3)
tsdiag(data1_fit4)
tsdiag(data1_fit5)
tsdiag(data1_fit6)
tsdiag(data1_fit7)
tsdiag(data1_fit8)
tsdiag(data1_fit9)
#tsdiag(data1_fit10)
#tsdiag(data1_fit11)

qqnorm(data1_fit9$resi,main="Proj1 Model ARMA(3,1) Q-Q plot ")


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

   plot(omega, 10*log10(my.spectrum), ylab="spectrum (in decibels)", type="l", main="Fitted Model Spectrum")   
}

par(mfrow=c(2,1))
   spec.pgram(fitting.data1)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data1_fit0$coef[1:2])), theta.of.b=c(1), variance=data1_fit0$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data1)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data1_fit1$coef[1:3])), theta.of.b=c(1), variance=data1_fit1$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data1)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data1_fit2$coef[1:4])), theta.of.b=c(1), variance=data1_fit2$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data1)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data1_fit3$coef[1:5])), theta.of.b=c(1), variance=data1_fit3$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data1)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data1_fit4$coef[1:6])), theta.of.b=c(1), variance=data1_fit4$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data1)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data1_fit5$coef[1:7])), theta.of.b=c(1), variance=data1_fit5$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data1)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data1_fit6$coef[1:8])), theta.of.b=c(1), variance=data1_fit6$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data1)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data1_fit7$coef[1:9])), theta.of.b=c(1), variance=data1_fit7$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data1)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data1_fit8$coef[1:2])), theta.of.b=c(1,data1_fit8$coef[3]), variance=data1_fit8$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data1, main="Original training data Periodogram")  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data1_fit9$coef[1:3])), theta.of.b=c(1,data1_fit9$coef[4]), variance=data1_fit9$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data1)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data1_fit15$coef[1:2])), theta.of.b=c(1,data1_fit15$coef[3:4]), variance=data1_fit15$sigma2)
par(mfrow=c(1,1))

## plot theoretical acf against the sample acf

par(mfrow=c(2,2))
acf(fitting.data1, type = c("correlation"))
acf(fitting.data1, type = c("partial"))
plot(ARMAacf(ar=c(data1_fit1$coef[1:3]),lag.max=25),type="h")
abline(h=0)
plot(ARMAacf(ar=c(data1_fit1$coef[1:3]),lag.max=25, pacf=TRUE),type="h")
abline(h=0)

par(mfrow=c(2,2))
acf(fitting.data1, type = c("correlation"))
acf(fitting.data1, type = c("partial"))
plot(ARMAacf(ar=c(data1_fit7$coef[1:9]),lag.max=25),type="h")
abline(h=0)
plot(ARMAacf(ar=c(data1_fit7$coef[1:9]),lag.max=25, pacf=TRUE),type="h")
abline(h=0)

par(mfrow=c(2,2))
acf(fitting.data1, type = c("correlation"),main="Original training data ACF")
acf(fitting.data1, type = c("partial"),main="Original training data pACF")
plot(ARMAacf(ar=c(data1_fit9$coef[1:3]),ma=c(data1_fit9$coef[4]),lag.max=25),type="h",ylab="acf",main="Model ARMA(3,1) acf")
abline(h=0)
plot(ARMAacf(ar=c(data1_fit9$coef[1:3]),ma=c(data1_fit9$coef[4]),lag.max=25, pacf=TRUE),type="h",main="Model ARMA(3,1) pacf")
abline(h=0)

par(mfrow=c(1,1))


## predict using Test data 
data1_pred0<- predict(data1_fit0, n.ahead=51, se.fit=FALSE)
data1_pred1<- predict(data1_fit1, n.ahead=51, se.fit=FALSE)
data1_pred2<- predict(data1_fit2, n.ahead=51, se.fit=FALSE)
data1_pred3<- predict(data1_fit3, n.ahead=51, se.fit=FALSE)
data1_pred4<- predict(data1_fit4, n.ahead=51, se.fit=FALSE)
data1_pred5<- predict(data1_fit5, n.ahead=51, se.fit=FALSE)
data1_pred6<- predict(data1_fit6, n.ahead=51, se.fit=FALSE)
data1_pred7<- predict(data1_fit7, n.ahead=51, se.fit=FALSE)
data1_pred8<- predict(data1_fit8, n.ahead=51, se.fit=FALSE)
data1_pred9<- predict(data1_fit9, n.ahead=51, se.fit=FALSE)
data1_pred15<- predict(data1_fit15, n.ahead=51, se.fit=FALSE)

data1_sse0 <- sum((test.data1 - data1_pred0)^2)
data1_sse1 <- sum((test.data1 - data1_pred1)^2)
data1_sse2 <- sum((test.data1 - data1_pred2)^2)
data1_sse3 <- sum((test.data1 - data1_pred3)^2)
data1_sse4 <- sum((test.data1 - data1_pred4)^2)
data1_sse5 <- sum((test.data1 - data1_pred5)^2)
data1_sse6 <- sum((test.data1 - data1_pred6)^2)
data1_sse7 <- sum((test.data1 - data1_pred7)^2)
data1_sse8 <- sum((test.data1 - data1_pred8)^2)
data1_sse9 <- sum((test.data1 - data1_pred9)^2)
data1_sse15 <- sum((test.data1 - data1_pred15)^2)
data1_sse0
data1_sse1
data1_sse2
data1_sse3
data1_sse4
data1_sse5
data1_sse6
data1_sse7
data1_sse8
data1_sse9
data1_sse15

## Plot and see the prediction
par(mfrow=c(2,1))
plot(data1_pred9,type="b")
plot(test.data1,type="b")
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
data1.fit <- arima(data1, order=c(3,0,1),method="ML")
data1.fit
sum(((data1.fit$res)^2))/(length(data1) -(3+1) -(3+1+1))

tsdiag(data1.fit)
qqnorm(data1.fit$resi)

## predict 13 points ahead
data1.pred <- predict(data1.fit, n.ahead=13, se.fit=TRUE)

data1.prediction <- data1.pred$pred
data1.pred.lowerbound <- data1.prediction - 2*data1.pred$se
data1.pred.upperbound <- data1.prediction + 2*data1.pred$se

cbind(data1.pred.lowerbound,data1.prediction,data1.pred.upperbound)
save <- rbind(data1.pred.lowerbound,data1.prediction,data1.pred.upperbound)
write.table(save, file="proj1.csv",sep=",")
## plot
plot(400:501, data1[400:501], xlim=c(400, 520),  type="b", xlab="Time", ylab="Value", main="Proj_1 Forecast of Future 13 Points")

lines(502:514, data1.prediction, type="b", col=2)
lines(502:514, data1.pred.lowerbound, type="l", col=2)
lines(502:514, data1.pred.upperbound, type="l", col=2)
plot(forecast(data1.fit,h=13),xlim=c(400,520))



