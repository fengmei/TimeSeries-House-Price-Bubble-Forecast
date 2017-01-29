install.packages("tseries")
install.packages("TTR")
install.packages("forecast")
library("forecast")
require("tseries") 
## Read data4
data4 <-scan(file="~/SJSU_ClASSES/Math265_TimeSeries/FinalProject/deposits.txt")

## explore data
head(data4)
summary(data4)
sapply(data4,class)

#names(data1) <- c("Proj1")
length(data4)

## plot data
plot.ts(data4,type="b",xlab="Year",ylab="Value",main="Proj 4 Sedimentary Layer Thickness Time Series")
plot.ts(data4[200:400],type="b")  ## looks not stationary
# stationary test: pass
adf.test(data4)

## check differenced data, looks stationary now

differenced.data4 <- diff(data4, differences = 1)
plot(differenced.data4, type="b")
adf.test(differenced.data4)

## acf, pacf
par(mfrow=c(2,1))
   acf(data4, type = c("correlation"),lag.max=100)
   acf(data4, type = c("partial"),lag.max=100)
par(mfrow=c(1,1))  ## not stationary

## Periodogram
spec.pgram(data4,taper=0.1) ## not accurate now

## differenced acf, pacf, pediodogram  ; much better, like AR, or ARMA
par(mfrow=c(2,1))
   acf(differenced.data4, type = c("correlation"),lag.max=100)
   acf(differenced.data4, type = c("partial"),lag.max=100)
par(mfrow=c(1,1))  

spec.pgram(differenced.data4,taper=0.1) ## 

## auto.arima fit a model
f <- auto.arima(data4,max.p=10,max.q=6,start.p=0,start.q=1, seasonal=TRUE)
f


## log? becomes seasonal
log.data4 <- log(data4)
plot(log.data4,type="b",xlab="Year",ylab="Value",main="Proj 4 Logged Sedimentary Layer Thickness")
adf.test(log.data4)
diff.log.data4 <- diff(log.data4)
plot(diff.log.data4, type="b",xlab="Year",ylab="Value",main="Proj 4 Differientiated Log Sedimentary Layer Thickness Time Series")
adf.test(diff.log.data4)

par(mfrow=c(2,2))
   acf(log.data4, type = c("correlation"),lag.max=100,,main="Log-Data4 pACF")
   acf(log.data4, type = c("partial"),lag.max=100,,main="Log-Data4 pACF")

   acf(diff.log.data4, type = c("correlation"),lag.max=100,,main="Differenciated log-Data4 ACF")
   acf(diff.log.data4, type = c("partial"),lag.max=100,,main="Differenciated log-Data4 pACF")
par(mfrow=c(1,1))  
par(mfrow=c(1,2))
spec.pgram(log.data4,taper=0.1,main="Log-Data4 Periodogram") ## 
spec.pgram(diff.log.data4,taper=0.1,main="Differenced Log-Data4 Periodogram") ## 
par(mfrow=c(1,1))  
## check the mean
mean(differenced.data4)




## so candidates models are ARIMA(0,1,1), ARMA(1,1,1),ARMA(2,1,1), ARMA(1,1,2),ARMA(2,1,2) fit the candidates model
## Split data to 90% for training and 10% for test 

fitting.data4 <- log.data4[1:560]
test.data4 <- log.data4[561:624]



fitting.data4.diff <-diff.log.data4[1:560]


## using Arima(), to include contant. https://www.otexts.org/fpp/8/7
data4_fit1 <- Arima(fitting.data4,order=c(0,1,1),method="ML",include.drift=TRUE)
data4_fit2 <- Arima(fitting.data4,order=c(1,1,1),method="ML",include.drift=TRUE)
data4_fit3 <- Arima(fitting.data4,order=c(2,1,1),method="ML",include.drift=TRUE)
data4_fit4 <- Arima(fitting.data4,order=c(1,1,2),method="ML",include.drift=TRUE)
data4_fit5 <- Arima(fitting.data4,order=c(2,1,2),method="ML",include.drift=TRUE)
data4_fit6 <- Arima(fitting.data4,order=c(2,1,0),method="ML",include.drift=TRUE)

data4_fit1
data4_fit2
data4_fit3
data4_fit4
data4_fit5
data4_fit6

## coefficient significance test show that ARIMA(2,1,2),ARIMA(2,1,1),ARIMA(1,1,2) are not significant
## Choose In the rest models

names(data4_fit1)
## check MSE

sum(((data4_fit1$res)^2))/(length(fitting.data4) -(1) -(1+1))
sum(((data4_fit2$res)^2))/(length(fitting.data4) -(2) -(2+1))
sum(((data4_fit3$res)^2))/(length(fitting.data4) -(3) -(3+1))
sum(((data4_fit4$res)^2))/(length(fitting.data4) -(3) -(3+1))
sum(((data4_fit5$res)^2))/(length(fitting.data4) -(4) -(4+1))


## diagnostics
tsdiag(data4_fit1) #N
tsdiag(data4_fit2) #Y
tsdiag(data4_fit3) #Y
tsdiag(data4_fit4) #Y
tsdiag(data4_fit5) #Y

qqnorm(data4_fit1$resi)
qqnorm(data4_fit2$resi,main="ARIMA(1,1,1) Normal QQ Plot")
qqnorm(data4_fit3$resi)
qqnorm(data4_fit4$resi)
qqnorm(data4_fit5$resi)
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
   spec.pgram(fitting.data4.diff)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1), theta.of.b=c(1,data4_fit1$coef[1]), variance=data4_fit1$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data4.diff,main="Original training data Periodogram")  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data4_fit2$coef[1])), theta.of.b=c(1,data4_fit2$coef[2]), variance=data4_fit2$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data4.diff)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data4_fit3$coef[1:2])), theta.of.b=c(1,data4_fit3$coef[3]), variance=data4_fit3$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data4.diff)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data4_fit4$coef[1])), theta.of.b=c(1,data4_fit4$coef[2:3]), variance=data4_fit4$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
   spec.pgram(fitting.data4.diff)  ##  tapr = .1 by default
   my.spectrum(phi.of.b=c(1, -(data4_fit5$coef[1:2])), theta.of.b=c(1,data4_fit5$coef[3:4]), variance=data4_fit5$sigma2)
par(mfrow=c(1,1))

## !!!plot theoretical acf against the sample acf

par(mfrow=c(2,2))
acf(fitting.data4.diff, type = c("correlation"),main="Original training data ACF")
acf(fitting.data4.diff, type = c("partial"),main="Original training data pACF")
plot(ARMAacf(ma=c(data4_fit1$coef[1]),lag.max=25),type="h")
abline(h=0)
plot(ARMAacf(ma=c(data4_fit1$coef[1]),lag.max=25, pacf=TRUE),type="h",ylab="acf",main="Model ARIMA(0,1,1) acf")
)
abline(h=0)

par(mfrow=c(2,2))
acf(fitting.data4.diff, type = c("correlation"),main="Original training data ACF")
acf(fitting.data4.diff, type = c("partial"),main="Original training data ACF")
plot(ARMAacf(ar=c(data4_fit2$coef[1]),ma=c(data4_fit2$coef[2]),lag.max=25),type="h",ylab="acf",main="Model ARIMA(1,1,1) acf")
abline(h=0)
plot(ARMAacf(ar=c(data4_fit2$coef[1]),ma=c(data4_fit2$coef[2]),lag.max=25, pacf=TRUE),type="h",ylab="acf",main="Model ARIMA(1,1,1) pacf")
abline(h=0)

par(mfrow=c(2,2))
acf(fitting.data4.diff, type = c("correlation"))
acf(fitting.data4.diff, type = c("partial"))
plot(ARMAacf(ar=c(data4_fit3$coef[1:2]),ma=c(data4_fit3$coef[3]),lag.max=25),type="h")
abline(h=0)
plot(ARMAacf(ar=c(data4_fit3$coef[1:2]),ma=c(data4_fit3$coef[3]),lag.max=25, pacf=TRUE),type="h")
abline(h=0)

par(mfrow=c(2,2))
acf(fitting.data4.diff, type = c("correlation"))
acf(fitting.data4.diff, type = c("partial"))
plot(ARMAacf(ar=c(data4_fit4$coef[1]),ma=c(data4_fit4$coef[2:3]),lag.max=25),type="h")
abline(h=0)
plot(ARMAacf(ar=c(data4_fit4$coef[1:2]),ma=c(data4_fit4$coef[2:3]),lag.max=25, pacf=TRUE),type="h")
abline(h=0)

par(mfrow=c(1,1))


## predict using Test data 

data4_pred1 <- forecast(data4_fit1, h=64)
data4_pred2 <- forecast(data4_fit2, h=64)
data4_pred3 <- forecast(data4_fit3, h=64)
data4_pred4 <- forecast(data4_fit4, h=64)
data4_pred5 <- forecast(data4_fit5, h=64)


data4_sse1 <- sum((test.data4 - data4_pred1$mean)^2)
data4_sse2 <- sum((test.data4 - data4_pred2$mean)^2)
data4_sse3 <- sum((test.data4 - data4_pred3$mean)^2)
data4_sse4 <- sum((test.data4 - data4_pred4$mean)^2)
data4_sse5 <- sum((test.data4 - data4_pred5$mean)^2)

data4_sse1
data4_sse2
data4_sse3
data4_sse4
data4_sse5

## Plot and see the prediction
par(mfrow=c(2,1))
plot(data4_pred2$mean,type="b")
plot(test.data4,type="b")
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
data4.fit <- Arima(log.data4, order=c(1,1,1), method="ML",include.drift=TRUE)
tsdiag(data4.fit)
qqnorm(data4.fit$resi)

sum(((data4.fit$res)^2))/(length(data3) -(3+1) -(3+1+1))
proj4forecast <- data4.prediction
proj4forecast <- rbind(exp(proj4forecast$lower[,2]),exp(proj4forecast$mean),exp(proj4forecast$upper[,2]))
write.table(proj4forecast, file="proj4.csv",sep=",")

## predict 13 points ahead
data4.prediction<- forecast(data4.fit,h=13)
plot(forecast(data4.fit,h=13),xlim=c(1,640))
exp(data4.prediction$mean)
exp(data4.prediction$lower[,2])
exp(data4.prediction$upper[,2])


plot(560:624, data4[560:624], xlim=c(560,640),type="b",xlab="Year",ylab="Value",main="Proj 4 Sedimentary Layer Thickness Forecast")
lines(625:637,exp(data4.prediction$mean),type="b",col=2)
lines(625:637,exp(data4.prediction$lower[,2]),type="l", col=2)
lines(625:637,exp(data4.prediction$upper[,2]),type="l", col=2)


