install.packages("tseries")
install.packages("TTR")
install.packages("forecast")
library("forecast")
require("tseries") 
## Read data3
data3 <-scan(file="~/SJSU_ClASSES/Math265_TimeSeries/FinalProject/proj3.txt")

## explore data
head(data3)
summary(data3)
sapply(data3,class)
length(data3)

## plot data
plot.ts(data3,type="b",main="PROJ 3 Time Series")  ## looks not stationary, with trend
## stationary test, can't pass - unit root test
adf.test(data3)

## check differenced data, looks stationary now

differenced.data3 <- diff(data3, differences = 1)
plot(differenced.data3, type="b")
adf.test(differenced.data3)

## acf, pacf
par(mfrow=c(1,2))
   acf(data3, type = c("correlation"),main="Data3 Sample ACF")
   acf(data3, type = c("partial"),main="Data3 Sample pACF")
par(mfrow=c(1,1))  ## not stationary

## Periodogram
spec.pgram(data3,taper=0.1) ## not accurate now

## differenced acf, pacf, pediodogram  ; much better, like AR, or ARMA
par(mfrow=c(1,2))
   acf(differenced.data3, type = c("correlation"), main="Differenciated Data3 ACF")
   acf(differenced.data3, type = c("partial"),main="Differenciated Data3 pACF")
par(mfrow=c(1,1))  

par(mfrow=c(1,2))
spec.pgram(data3,taper=0.1,main="Data3 Periodogram")
spec.pgram(differenced.data3,taper=0.1,main="Differenciated Data3 Periodogram") ## looks like the data1....
par(mfrow=c(1,1))  


## auto.arima fit a model
f <- auto.arima(data3,max.p=10,max.q=6,start.p=2,start.q=1)
f


## log? No!
log.data3 <- log(data3)
plot(log.data3,type="b")

## check and subtract the mean
##data1_zeromean <- data1 - mean(data1)
mean(differenced.data3)

## so candidates models are AR(), AR(),ARMA(),ARMA(), ARMA(), fit the candidates model
#ar(data1,aic=FALSE, order.max=4, method=c("mle"))
#ar(data1,aic=FALSE, order.max=4, method=c("ols"))

## Split data to 90% for training and 10% for test 

fitting.data3 <- data3[1:520]
test.data3 <- data3[521:584]
fitting.train.data3 <- differenced.data3[1:520]
test.data3 <- differenced.data3[521:584]

## using Arima(), to include contant. https://www.otexts.org/fpp/8/7
data3_fit1 <- Arima(fitting.data3,order=c(2,1,0),method="ML",include.drift=TRUE)
data3_fit2 <- Arima(fitting.data3,order=c(3,1,0),method="ML",include.drift=TRUE)
data3_fit3 <- Arima(fitting.data3,order=c(5,1,0),method="ML",include.drift=TRUE)
data3_fit4 <- Arima(fitting.data3,order=c(8,1,0),method="ML",include.drift=TRUE)
data3_fit5 <- Arima(fitting.data3,order=c(10,1,0),method="ML",include.drift=TRUE)
data3_fit6 <- Arima(fitting.data3,order=c(2,1,1),method="ML",include.drift=TRUE)
data3_fit7 <- Arima(fitting.data3,order=c(3,1,1),method="ML",include.drift=TRUE)
data3_fit8 <- Arima(fitting.data3,order=c(2,1,2),method="ML",include.drift=TRUE)
data3_fit9 <- Arima(fitting.data3,order=c(3,1,2),method="ML",include.drift=TRUE)




data3_fit1  ## CSS and ML got the very close results, so choose ML method
data3_fit2
data3_fit3
data3_fit4
data3_fit5
data3_fit6
data3_fit7
data3_fit8
data3_fit9

## coefficient significance test show that ARMA(3,2),ARMA(2,2),AR(10) are not significant
## Choose In the rest models

names(data3_fit1)
## check MSE

sum(((data3_fit1$res)^2))/(length(fitting.data3) -(2) -(2+1))
sum(((data3_fit2$res)^2))/(length(fitting.data3) -(3) -(3+1))
sum(((data3_fit3$res)^2))/(length(fitting.data3) -(5) -(5+1))
sum(((data3_fit4$res)^2))/(length(fitting.data3) -(8) -(8+1))
sum(((data3_fit5$res)^2))/(length(fitting.data3) -(10) -(10+1))
sum(((data3_fit6$res)^2))/(length(fitting.data3) -(2+1) -(2+1+1))
sum(((data3_fit7$res)^2))/(length(fitting.data3) -(3+1) -(3+1+1))
sum(((data3_fit8$res)^2))/(length(fitting.data3) -(2+2) -(2+2+1))
sum(((data3_fit9$res)^2))/(length(fitting.data3) -(3+2) -(3+2+1))



## diagnostics
tsdiag(data3_fit1) #N
tsdiag(data3_fit2) #N
tsdiag(data3_fit3) #N
tsdiag(data3_fit4) #Y
tsdiag(data3_fit5) #Y
tsdiag(data3_fit6) #N
tsdiag(data3_fit7) #Y
tsdiag(data3_fit8) #N
tsdiag(data3_fit9) #Y

qqnorm(data3_fit4$resi)
qqnorm(data3_fit7$resi, main="ARIMA(3,1,1) QQ Plot")

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
   spec.pgram(fitting.data3)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data3_fit1$coef[1:2])), theta.of.b=c(1), variance=data3_fit1$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.train.data3)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data3_fit2$coef[1:3])), theta.of.b=c(1), variance=data3_fit2$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(diff(fitting.train.data3))  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data3_fit3$coef[1:5])), theta.of.b=c(1), variance=data3_fit3$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.train.data3)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data3_fit4$coef[1:8])), theta.of.b=c(1), variance=data3_fit4$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data3)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data3_fit5$coef[1:10])), theta.of.b=c(1), variance=data3_fit5$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.train.data3)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data3_fit6$coef[1:2])), theta.of.b=c(1,data3_fit6$coef[3]), variance=data3_fit6$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.train.data3,main="Original training data Periodogram")  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data3_fit7$coef[1:3])), theta.of.b=c(1,data3_fit7$coef[4]), variance=data3_fit7$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.train.data3)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data3_fit8$coef[1:2])), theta.of.b=c(1,data3_fit8$coef[3:4]), variance=data3_fit8$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data3)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data3_fit9$coef[1:3])), theta.of.b=c(1,data3_fit9$coef[4:5]), variance=data3_fit9$sigma2)
par(mfrow=c(1,1))

## !!!NOT HERE.  plot theoretical acf against the sample acf

par(mfrow=c(2,2))
acf(fitting.train.data3, type = c("correlation"))
acf(fitting.train.data3, type = c("partial"))
plot(ARMAacf(ar=c(data3_fit4$coef[1:8]),lag.max=25),type="h")
abline(h=0)
plot(ARMAacf(ar=c(data3_fit4$coef[1:8]),lag.max=25, pacf=TRUE),type="h")
abline(h=0)

par(mfrow=c(2,2))
acf(fitting.train.data3, type = c("correlation"),main="Original training data ACF")
acf(fitting.train.data3, type = c("partial"),main="Original training data pACF")
plot(ARMAacf(ar=c(data3_fit8$coef[1:2]),ma=c(data3_fit8$coef[3:4]),lag.max=25),type="h",ylab="acf",main="Model ARIMA(3,1,1) acf")
abline(h=0)
plot(ARMAacf(ar=c(data3_fit8$coef[1:2]),ma=c(data3_fit8$coef[3:4]),lag.max=25, pacf=TRUE),type="h",ylab="acf",main="Model ARMA(3,1,1) pacf")
abline(h=0)
par(mfrow=c(1,1))


## predict using Test data 

data3_pred1 <- forecast(data3_fit1, h=64)
data3_pred2 <- forecast(data3_fit2, h=64)
data3_pred3 <- forecast(data3_fit3, h=64)
data3_pred4 <- forecast(data3_fit4, h=64)
data3_pred5 <- forecast(data3_fit5, h=64)
data3_pred6 <- forecast(data3_fit6, h=64)
data3_pred7 <- forecast(data3_fit7, h=64)
data3_pred8 <- forecast(data3_fit8, h=64)
data3_pred9 <- forecast(data3_fit9, h=64)

data3_sse1 <- sum((test.data3 - data3_pred1$mean)^2)
data3_sse2 <- sum((test.data3 - data3_pred2$mean)^2)
data3_sse3 <- sum((test.data3 - data3_pred3$mean)^2)
data3_sse4 <- sum((test.data3 - data3_pred4$mean)^2)
data3_sse5 <- sum((test.data3 - data3_pred5$mean)^2)
data3_sse6 <- sum((test.data3 - data3_pred6$mean)^2)
data3_sse7 <- sum((test.data3 - data3_pred7$mean)^2)
data3_sse8 <- sum((test.data3 - data3_pred8$mean)^2)
data3_sse9 <- sum((test.data3 - data3_pred9$mean)^2)

data3_sse1
data3_sse2
data3_sse3
data3_sse4
data3_sse5
data3_sse6
data3_sse7
data3_sse8
data3_sse9

## Plot and see the prediction
par(mfrow=c(2,1))
plot(data3_pred4$mean,type="b")
plot(test.data3,type="b")
par(mfrow=c(1,1))

## COMPARE Models         AR(2)   AR(3)   AR(5)  AR(8)  AR(10)   ARMA(2,1)   ARMA(3,1)    ARMA(2,2)    ARMA(3,2)  
## Coef Signif Test	   |   Y		Y		Y	   Y	  Y		  Y		  Y		  Y          Y		    Y			
##  AIC                |5075.41 4924.98 4850.83 4824.11 4802.03 4796.08 4783.13  4782.34   5009.06    4758.6 
## MSE				   |4572.86 3266.52 2766.15 2605.04 2478.72 2445.71 2375.2  2371.08   3942.97    2243.53 
## Periodogram         |   L		L		H	   H	  M		  M		  M		  M			 L          H					
## test prediction sse |        3200350  3061494 3014249 3091009 2999213 3127040 3047744   3312724    2975883
## Parsinomy		   |   Y		Y		Y	   Y	  Y		  N		  N		  N			  Y 		Y				
## Diagnostics		   |   N		N		N	   M	  M		  Y       Y       Y           N         Y

## Fit the model to the entire dataset
data3.fit <- Arima(data3, order=c(3,1,1),method="ML",include.drift=TRUE)

sum(((data3.fit$res)^2))/(length(data3) -(3+1) -(3+1+1))
proj3forecast <- rbind(data3forecast$lower[,2],data3forecast$mean,data3forecast$upper[,2])
write.table(proj3forecast, file="proj3.csv",sep=",")
tsdiag(data3.fit)
qqnorm(data3.fit$resi)

## predict 13 points ahead
forecast(data3.fit,h=13)
plot(forecast(data3.fit,h=13),type="b",xlim=c(550,600),main="PROJ 3 13-points Forecast")



