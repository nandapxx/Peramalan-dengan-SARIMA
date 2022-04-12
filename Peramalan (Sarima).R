## Install Packages
install.packages("tseries")
install.packages("TSA")
install.packages("timeseries")
install.packages("starma")
install.packages("MTS")
install.packages("fGarch")
install.packages("portes")
install.packages("MVN")
install.packages("readr")
install.packages("tsoutliers")
install.packages("FitAR")
library(tseries)
library(TSA)
library(timeSeries)
library(starma)
library(MTS)
library(fGarch)
library(portes)
library(MVN)
library(readr)
library(tsoutliers)
library(FitAR)

## Input Data
data= read.csv(file.choose(),header=T,sep=";")
attach(data)
X1=data$WNPMI

## Plot
data_X1=ts(X1,start=c(2011,1),frequency=12)
ts.plot(data_X1,ylab="",main="PLOT DATA WNPMI",col="magenta",lwd=2)

# Cek stasioner dalam varian
BoxCox.lambda(data_X1)

# ADF Test X1
adf.test(data_X1)
acf(data_X1)
pacf(data_X1)

# Differencing musiman
X1d=diff(data_X1,lag=12,differences = 1)
ts.plot(X1d)
adf.test(X1d)

## Uji ACF PACF WNPMI
# ACF PACF (X1)
acf(X1d, main="ACF WNPMI (D=1)",lag=40)
axis(1,at=seq(1,40,by=1))
pacf(X1d, main="PACF WNPMI (D=1)",lag=40)
axis(1,at=seq(5,40,by=1))

# Differencing 1
dX1=diff(data_X1)
plot.ts(dX1)
adf.test(dX1)
acf(dX1)
pacf(dX1)

## Differencing musiman 
X1t=diff(dX1,lag=12,differences = 1)
ts.plot(X1t)
adf.test(X1t)

## Uji ACF PACF WNPMI
# ACF PACF (X1)
acf(X1t, main="ACF WNPMI (D=1)",lag=40)
axis(1,at=seq(1,40,by=1))
pacf(X1t, main="PACF WNPMI (D=1)",lag=40)
axis(1,at=seq(5,40,by=1))

## Estimasi Parameter WNPMI
b1=Arima(data_X1,order=c(1,1,1),seasonal=list(order=c(1,1,0),period=12));b1
b2=Arima(data_X1,order=c(2,1,1),seasonal=list(order=c(1,1,0),period=12));b2
b3=Arima(data_X1,order=c(1,1,1),seasonal=list(order=c(1,1,1),period=12));b3
b4=Arima(data_X1,order=c(2,1,1),seasonal=list(order=c(1,1,1),period=12));b4
b5=Arima(data_X1,order=c(1,1,1),seasonal=list(order=c(2,1,1),period=12));b5
b6=Arima(data_X1,order=c(2,1,1),seasonal=list(order=c(2,1,1),period=12));b6
b7=Arima(data_X1,order=c(1,1,1),seasonal=list(order=c(0,1,1),period=12));b7
b8=Arima(data_X1,order=c(2,1,1),seasonal=list(order=c(0,1,1),period=12));b8
b9=Arima(data_X1,order=c(1,1,1),seasonal=list(order=c(2,1,0),period=12));b9
b10=Arima(data_X1,order=c(2,1,1),seasonal=list(order=c(2,1,0),period=12));b10
b11=Arima(data_X1,order=c(0,1,1),seasonal=list(order=c(1,1,0),period=12));b11
b12=Arima(data_X1,order=c(2,1,0),seasonal=list(order=c(1,1,0),period=12));b12
b13=Arima(data_X1,order=c(0,1,1),seasonal=list(order=c(1,1,1),period=12));b13
b14=Arima(data_X1,order=c(2,1,0),seasonal=list(order=c(1,1,1),period=12));b14
b15=Arima(data_X1,order=c(0,1,1),seasonal=list(order=c(2,1,1),period=12));b15
b16=Arima(data_X1,order=c(2,1,0),seasonal=list(order=c(2,1,1),period=12));b16
b17=Arima(data_X1,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12));b17
b18=Arima(data_X1,order=c(2,1,0),seasonal=list(order=c(0,1,1),period=12));b18
b19=Arima(data_X1,order=c(0,1,1),seasonal=list(order=c(2,1,0),period=12));b19
b20=Arima(data_X1,order=c(2,1,0),seasonal=list(order=c(2,1,0),period=12));b20
b21=Arima(data_X1,order=c(1,1,0),seasonal=list(order=c(1,1,0),period=12));b21
b22=Arima(data_X1,order=c(1,1,0),seasonal=list(order=c(1,1,1),period=12));b22
b23=Arima(data_X1,order=c(1,1,0),seasonal=list(order=c(2,1,1),period=12));b23
b24=Arima(data_X1,order=c(1,1,0),seasonal=list(order=c(0,1,1),period=12));b24
b25=Arima(data_X1,order=c(1,1,0),seasonal=list(order=c(2,1,0),period=12));b25
#hasil: b7 dengan nilai aic terkecil

## UJI SIGNIFIKANSI PARAMETER
printstatarima <- function (x, digits = 4,se=T,...){
  if (length(x$coef) > 0) {
    cat("\nCoefficients:\n")
    coef <- round(x$coef, digits = digits)
    if (se && nrow(x$var.coef)) {
      ses <- rep(0, length(coef))
      ses[x$mask] <- round(sqrt(diag(x$var.coef)), digits = digits)
      coef <- matrix(coef, 1, dimnames = list(NULL, names(coef)))
      coef <- rbind(coef, s.e. = ses)
      statt <- coef[1,]/ses
      pval  <- 2*pt(abs(statt), df=length(x$residuals)-1, lower.tail = F)
      coef <- rbind(coef, t=round(statt,digits=digits),sign.=round(pval,digits=digits))
      coef <- t(coef)
    }
    print.default(coef, print.gap = 2)
  }
}

## Uji Signifikansi
printstatarima(b7)