#import our libraries
library(tidyverse)
library(VineCopula)
library(fGarch)
library(KScorrect)
library(stats)
library(ADGofTest)
library(tseries)
library(fBasics)
library(MASS)

#import the historical data
n225 <- read.csv("~/STAT/ICA/N225.csv", stringsAsFactors=FALSE)
fchi <- read.csv("~/STAT/ICA/FCHI.csv", stringsAsFactors=FALSE)

#creates a function returning log returns
n225_r <- function(t) {
  return (log(as.numeric(n225[t+1,6]))-log(as.numeric(n225[t,6])))
}
fchi_r <- function(t) {
  return (log(as.numeric(fchi[t+1,6]))-log(as.numeric(fchi[t,6])))
}
#construct a data set with adjclose values as well as log returns
numOfWeeks = nrow(n225)
returnsDataSet <- data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c("Week", "n225_returns", "fchi_returns", "n225_adjclose", "fchi_adjclose"))))

for (i in 1:numOfWeeks-1){
  returnsDataSet[i,1]=i
  returnsDataSet[i,2]=n225_r(i)
  returnsDataSet[i,3]=fchi_r(i)
  returnsDataSet[i,4]=as.numeric(n225[i+1,6])
  returnsDataSet[i,5]=as.numeric(fchi[i+1,6])
 
}
returnsDataSet <- na.omit(returnsDataSet)

print(head(returnsDataSet))



#the identification step to find the order p 
print(jarque.bera.test(returnsDataSet$n225_returns))
#create plots to check for autocorrelation
par(mfrow=c(3,2))
plot(returnsDataSet$Week, returnsDataSet$n225_returns, type= "l")
plot(returnsDataSet$n225_adjclose, type="l")
acf(returnsDataSet$n225_returns, col="green", lwd=2)
pacf(returnsDataSet$n225_returns, col="green", lwd=2)
acf(returnsDataSet$n225_returns^2, col="red", lwd=2)
pacf(returnsDataSet$n225_returns^2, col="red", lwd=2)
#so we will use 0 lags.
#we do the same for the fchi
par(mfrow=c(3,2))
plot(returnsDataSet$Week, returnsDataSet$fchi_returns, type= "l")
plot(returnsDataSet$fchi_adjclose, type="l")
acf(returnsDataSet$fchi_returns, col="green", lwd=2)
pacf(returnsDataSet$fchi_returns, col="green", lwd=2)
acf(returnsDataSet$fchi_returns^2, col="red", lwd=2)
pacf(returnsDataSet$fchi_returns^2, col="red", lwd=2)
#so we will use 1 lag. 
#which conditional distribution should we use for the garch model?
distributions<- c("norm", "std", "sstd", "snorm", "sged", "ged")
cond_dist_table1 <- data.frame(matrix(ncol=3, nrow=0, dimnames=list(NULL, c("dist", "aic", "bic"))))

for (k in 1:6){  
  model1=garchFit(formula=~arma(3,0)+garch(1,1),data=returnsDataSet$n225_returns,trace=F,cond.dist=distributions[k])
  cond_dist_table1[k,1]=distributions[k]
  
  cond_dist_table1[k,2]=(summary(model1))$ics["AIC"]
  cond_dist_table1[k,3]=(summary(model1))$ics["BIC"]
  }

cond_dist_table1
#so we use the sstd
cond_dist_table2 <- data.frame(matrix(ncol=3, nrow=0, dimnames=list(NULL, c("dist", "aic", "bic"))))
for (k in 1:6){  
  model1=garchFit(formula=~garch(1,1),data=returnsDataSet$fchi_returns,trace=F,cond.dist=distributions[k])
  cond_dist_table2[k,1]=distributions[k]
  
  cond_dist_table2[k,2]=(summary(model1))$ics["AIC"]
  cond_dist_table2[k,3]=(summary(model1))$ics["BIC"]
}

cond_dist_table2

#we choose also sstd


#we conclude by fitting an garch(1,1) model n225, and ar(1)garch(1,1) to model fchi
model1=garchFit(formula=~garch(1,1),data=returnsDataSet$n225_returns,trace=F,cond.dist="sstd")
model2=garchFit(formula=~arma(1,0)+garch(1,1),data=returnsDataSet$fchi_returns,trace=F,cond.dist="sstd")
model1coef <- model1@fit$coef
model2coef <- model2@fit$coef


# Step 3: Model checking
res1 <- residuals(model1, standardize=TRUE)
par(mfrow=c(2,2))
acf(res1, col="green", lwd=2)
acf(res1^2, col="red", lwd=2)
pacf(res1, col="green", lwd=2)
pacf(res1^2, col="red", lwd=2)
par(mfrow=c(1,1))

Box.test(res1, lag = 10, type = c("Ljung-Box"), fitdf = 1)#passed
Box.test(res1^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)#passed


u1<-psstd(res1, nu= model1coef["shape"], xi= model1coef["skew"])
#Kolmogorov-Smirnov test
KStest1<-LcKS(u1, cdf = "punif")
KStest1$p.value#passed
#Anderson-Darling test
ADtest1<-ad.test(u1, null="punif")
ADtest1$p.value#passed


res2 <- residuals(model2, standardize=TRUE)
par(mfrow=c(2,2))
acf(res2, col="green", lwd=2)
acf(res2^2, col="red", lwd=2)
pacf(res2, col="green", lwd=2)
pacf(res2^2, col="red", lwd=2)
par(mfrow=c(1,1))
Box.test(res2, lag = 10, type = c("Ljung-Box"), fitdf = 1)#passed
Box.test(res2^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)#passed

#suggests we have captured autocorrelation

u2<-psstd(res2, nu= model2coef["shape"], xi= model2coef["skew"])
#Kolmogorov-Smirnov test
KStest2<-LcKS(u2, cdf = "punif")
KStest2$p.value#passed
#Anderson-Darling testhttp://127.0.0.1:16565/graphics/7bfc49bb-419c-4527-be45-5dea6b62c5ca.png
ADtest2<-ad.test(u2, null="punif")
ADtest2$p.value#passed
#histograms
par(mfrow=c(2,2))
hist1_x <- seq(min(res1, res2), max(res1, res2), length.out = 100)

hist(u1, breaks=50)
hist(u2, breaks=50)

hist(res1, prob=TRUE,  breaks=50)
lines(hist1_x, dsstd(hist1_x, nu = model1coef["shape"], xi = model1coef["skew"]), add = TRUE, col = "red")
hist(res2, prob= TRUE, breaks=50)
lines(hist1_x, dsstd(hist1_x, nu = model2coef["shape"], xi = model2coef["skew"]), add = TRUE, col = "red")





#generating epsilon using copulas
copModel=BiCopSelect(u1, u2, familyset=NA, selectioncrit="AIC", indeptest=TRUE, level=0.05,se = TRUE)
copModel
copModel2=BiCopSelect(u1, u2, familyset=NA, selectioncrit="BIC", indeptest=TRUE, level=0.05,se = TRUE)
copModel2

#the best model to use is Bivariate copula: BB1 (par = 0.54, par2 = 1.26, tau = 0.37)  so we will use BB1
N=10000


u_sim=BiCopSim(N, family=copModel$family, copModel$par,  copModel$par2)
#check if plots of observed vs simulated copula
par(mfrow=c(2,1))
plot(u1,u2)
plot(u_sim[,1], u_sim[,2], col=rgb(0,100,0,50,maxColorValue=255), pch=16)
res1_sim=qsstd(u_sim[,1], nu=model1coef["shape"], xi=model1coef["skew"])
res2_sim=qsstd(u_sim[,2], nu=model2coef["shape"], xi=model2coef["skew"])

print(tail(u1, n=1))
#introduce AR GARCH effects
varEstimate1 <- slot(model1, "sigma.t")
#generate unstandardized residuals
uns_res1<-residuals(model1, standardize=FALSE)
uns_res2<-residuals(model2, standardize=FALSE)
sigmaSquaredHat1=model1coef["omega"] + model1coef["alpha1"] * (tail(uns_res1, n=1))^2+ model1coef["beta1"] * (tail(varEstimate1, n=1)^2)
tail=as.numeric(tail(returnsDataSet$fchi_returns, n=1))

sim_values1 = (model1coef["mu"]+sqrt(sigmaSquaredHat1)*res1_sim)
#check histograms of simulated returns against past returns
par(mfrow=c(2,1))
hist(sim_values1)
hist(returnsDataSet$n225_returns)

#and GARCH EFFECTS for the other model
varEstimate2 <- slot(model2, "sigma.t")

sigmaSquaredHat2=model2coef["omega"] + model2coef["alpha1"] * (tail(uns_res2, n=1))^2+ model2coef["beta1"] * (tail(varEstimate2, n=1)^2)
sqrt(sigmaSquaredHat1)

sim_values2 = (model2coef["mu"]+model2coef["ar1"]*tail+sqrt(sigmaSquaredHat2)*res2_sim)
#histograms of simulated returns vs past
par(mfrow=c(2,1))
hist(sim_values2)
hist(returnsDataSet$fchi_returns)
#check plots of simulated values vs past values:
par(mfrow=c(2,1))
plot(sim_values1, sim_values2)
plot(returnsDataSet$n225_returns, returnsDataSet$fchi_returns)
#find var
VaR=quantile(0.5*(sim_values1+sim_values2), c(0.05,0.01))

VaR