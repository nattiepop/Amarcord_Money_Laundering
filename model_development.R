# Set directory

setwd("~/Master/Time Serie Analysis/Case 2")
library(forecast)
library(TSA)

## Try hourly incoming wires ##

data <- read.csv("amarcord_1hr.csv")
hourly <- ts(data$Trans_Amt, frequency = 7)
plot(hourly)

par(mfrow=c(2,1))
acf(hourly, main="")
pacf(hourly, main="")

decompose(hourly)
plot(decompose(hourly))

library(tseries)
adf.test(hourly)

# ETS

model.ets = ets(hourly)
par(mfrow=c(1,1))
plot(forecast(model.ets, h = 49))

summary(model.ets)

plot(forecast(model.ets, h=49)$fitted)

fit <- forecast(model.ets, h=49)$fitted

plot(hourly)
lines(fit, col='red')

accuracy(hourly, fit)

upper <- fitted(model.ets) + 1.96*(model.ets$sigma)
lower <- fitted(model.ets) - 1.96*(model.ets$sigma)
plot(hourly, type="n")
polygon(c(time(hourly),rev(time(hourly))), c(upper,rev(lower)), 
        col=rgb(0,0,0.6,0.2), border=FALSE)
lines(hourly)
lines(fitted(model.ets),col='red')
out <- (hourly < lower | hourly > upper)
points(time(hourly)[out], hourly[out], pch=19, col='blue')

# sigma2 <- var(hourly)
# 
# upper <- fitted(model.ets) + 1.96*sqrt(sigma2)
# lower <- fitted(model.ets) - 1.96*sqrt(sigma2)
# plot(hourly, type="n")
# polygon(c(time(hourly),rev(time(hourly))), c(upper,rev(lower)),
#         col=rgb(0,0,0.6,0.2), border=FALSE)
# lines(hourly)
# lines(fitted(model.ets),col='red')
# out <- (hourly < lower | hourly > upper)
# points(time(hourly)[out], hourly[out], pch=19, col='blue')

## Try deseasonal hourly incoming wires ##

data <- read.csv("amarcord_1hr.csv")
hourly <- ts(data$Trans_Amt, frequency = 7)
plot(hourly)

par(mfrow=c(2,1))
plot(hourly)
acf(hourly)

decompose(hourly)
plot(decompose(hourly))

# seasonal <- decompose(hourly)$season
# deseasonal_hourly <- hourly - seasonal
# 
# hourly <- ts(deseasonal_hourly, frequency = 1)

# ARIMA

model.test <- auto.arima(hourly, stepwise = T, approximation = F, trace = T)

summary(model.test)

# model1 <- Arima(hourly,order=c(2,0,2),method='ML')
# model2 <- Arima(hourly,order=c(1,0,0),method='ML')
# model3 <- Arima(hourly,order=c(0,0,1),method='ML')
# model4 <- Arima(hourly,order=c(1,0,1),method='ML')

model1 <- Arima(hourly,order=c(0,0,0),seasonal=c(2,0,0),method='ML')
model2 <- Arima(hourly,order=c(1,0,0),seasonal=c(2,0,0),method='ML')
model3 <- Arima(hourly,order=c(0,0,1),seasonal=c(2,0,0),method='ML')
model4 <- Arima(hourly,order=c(1,0,1),seasonal=c(2,0,0),method='ML')

list.model <- list(model1,model2,model3,model4)

for (i in list.model){
  print(summary(i))
}

for (i in list.model){
  par(mfrow=c(1,2))
  hist(rstandard(i),xlab="Standardised residuals",main="")
  qqnorm(rstandard(i))
  qqline(rstandard(i))
}

for (i in list.model){
  forec = forecast(i, h=49)
  par(mfrow=c(1,1))
  plot(forec)
}

tsdiag(model4,gof.lag=20)

for (i in list.model){
  plot(forecast(i, h=49)$fitted)
}

for (i in list.model){
  fit <- forecast(i)$fitted
  par(mfrow=c(1,1))
  plot(hourly)
  lines(fit, col='red')
}

for (i in list.model){
  print(accuracy(hourly, forecast(i)$fitted))
}

upper <- fitted(model4) + 1.96*sqrt(model4$sigma2)
lower <- fitted(model4) - 1.96*sqrt(model4$sigma2)
plot(hourly, type="n", ylim=range(lower,upper))
polygon(c(time(hourly),rev(time(hourly))), c(upper,rev(lower)), 
        col=rgb(0,0,0.6,0.2), border=FALSE)
lines(hourly)
lines(fitted(model4),col='red')
out <- (hourly < lower | hourly > upper)
points(time(hourly)[out], hourly[out], pch=19, col='blue')

# Decomposition

data <- read.csv("amarcord_1hr.csv")
data$inc <- data$Trans_Amt
data

lm <- lm(inc~t, data=data)
summary(lm)

data$inc.trend <- (lm$coefficients[1] + (lm$coefficients[2] * data$t))
data$inc.detrend <- data$inc/data$inc.trend

df <- data

seasonal <- aggregate(df$inc.detrend, list(df$dow), FUN=mean) 
seasonal <- seasonal[,2]+(1-mean(seasonal[,2]))

data$seasonal <- rep(seasonal,30)

data$inc.multidecompose <- data$inc.trend*data$seasonal

accuracy(data$inc.multidecompose, data$inc)

plot(data$inc, type = 'l')
lines(data$inc.multidecompose, col='red')

sigma2 <- var(data$inc.multidecompose)

upper <- data$inc.multidecompose + 1.96*sqrt(sigma2)
lower <- data$inc.multidecompose - 1.96*sqrt(sigma2)
plot(data$inc, type="n")
polygon(c(time(data$inc),rev(time(data$inc))), c(upper,rev(lower)), 
        col=rgb(0,0,0.6,0.2), border=FALSE)
lines(data$inc)
lines(data$inc.multidecompose,col='red')
out <- (data$inc < lower | data$inc > upper)
points(time(data$inc)[out], data$inc[out], pch=19, col='blue')

# NNAR

data <- read.csv("amarcord_1hr.csv")
hourly <- ts(data$Trans_Amt, frequency = 7)
plot(hourly)

model.nnar = nnetar(hourly)
summary(model.nnar)
nnetforecast <- forecast(model.nnar, h = 35, PI = F) 
plot(nnetforecast)

fit <- fitted(model.nnar)

plot(hourly)
lines(fit, col='red')

accuracy(hourly, fit)

sigma2 <- var(hourly)

upper <- fitted(model.nnar) + 1.96*sqrt(sigma2)
lower <- fitted(model.nnar) - 1.96*sqrt(sigma2)
plot(hourly, type="n")
polygon(c(time(hourly),rev(time(hourly))), c(upper,rev(lower)), 
        col=rgb(0,0,0.6,0.2), border=FALSE)
lines(hourly)
lines(fitted(model.nnar),col='red')
out <- (hourly < lower | hourly > upper)
points(time(hourly)[out], hourly[out], pch=19, col='blue')


