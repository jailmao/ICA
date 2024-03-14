#install.packages("tidyverse") 
#uncomment if you dont have it installed
library(tidyverse)
library(VineCopula)
library(fGarch)
library(KScorrect)
library(stats)
library(ADGofTest)
library(tseries)

n225 <- read.csv("~/STAT/ICA/N225.csv", stringsAsFactors=FALSE)
fchi <- read.csv("~/STAT/ICA/FCHI.csv", stringsAsFactors=FALSE)

n225_r <- function(t) {
  return (log(as.numeric(n225[t+1,6]))-log(as.numeric(n225[t,6])))
}
fchi_r <- function(t) {
  return (log(as.numeric(fchi[t+1,6]))-log(as.numeric(fchi[t,6])))
}

numOfWeeks = nrow(n225)
returnsDataSet <- data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c("Week", "n225_returns", "fchi_returns", "n225_adjclose", "fchi_adjclose"))))

for (i in 1:numOfWeeks-1){
  returnsDataSet[i,1]=i
  returnsDataSet[i,2]=n225_r(i)
  returnsDataSet[i,3]=fchi_r(i)
  returnsDataSet[i,4]=as.numeric(n225[i,6])
  returnsDataSet[i,5]=as.numeric(fchi[i,6])
 
}
returnsDataSet <- na.omit(returnsDataSet)

print(head(returnsDataSet))



#the identification step to find the order p 
print(jarque.bera.test(returnsDataSet$n225_returns))

par(mfrow=c(3,2))
plot(returnsDataSet$Week, returnsDataSet$n225_returns, type= "l")
plot(returnsDataSet$n225_adjclose, type="l")
acf(returnsDataSet$n225_returns, col="green", lwd=2)
pacf(returnsDataSet$n225_returns, col="green", lwd=2)
acf(returnsDataSet$n225_returns^2, col="red", lwd=2)
pacf(returnsDataSet$n225_returns^2, col="red", lwd=2)
#so we will use 3 lags.
par(mfrow=c(3,2))
plot(returnsDataSet$Week, returnsDataSet$fchi_returns, type= "l")
plot(returnsDataSet$fchi_adjclose, type="l")
acf(returnsDataSet$fchi_returns, col="green", lwd=2)
pacf(returnsDataSet$fchi_returns, col="green", lwd=2)
acf(returnsDataSet$fchi_returns^2, col="red", lwd=2)
pacf(returnsDataSet$fchi_returns^2, col="red", lwd=2)
#so we will use 12 lags. (?)

modeln225=garchFit(formula=~arma(3,0)+garch(1,1),data=returnsDataSet$n225_returns,trace=F,cond.dist="norm")
modelfchi=garchFit(formula=~arma(12,0)+garch(1,1),data=returnsDataSet$fchi_returns,trace=F,cond.dist="norm")
print(modeln225)


# Step 3: Model checking
res1 <- residuals(modeln225, standardize=TRUE)
par(mfrow=c(2,1))
acf(res1, col="green", lwd=2)
acf(res1^2, col="red", lwd=2)
par(mfrow=c(1,1))
Box.test(res1, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res1^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)

u1<-pnorm(res1, mean=0, sd=1)[4:numOfWeeks]
hist(u1)
