#install.packages("fBasics") 
#uncomment if you dont have it installed
library(tidyverse)
library(VineCopula)
library(fGarch)
library(KScorrect)
library(stats)
library(ADGofTest)
library(tseries)
library(fBasics)
library(MASS)
library(quantmod)
library(tidyquant)
#import the historical data
options("getSymbols.warning4.0"=FALSE)
options("getSymbols.yahoo.warning"=FALSE)
# Downloading ftse price using quantmod

getSymbols("^FTSE", from = '2000-01-01',
           to = "2023-01-01",warnings = FALSE,
           auto.assign = TRUE, periodicity="weekly")
getSymbols("^N225", from = '2000-01-01',
           to = "2023-01-01",warnings = FALSE,
           auto.assign = TRUE, periodicity="weekly")

#creates a function returning log returns
returns <- function(input, t){
  return (log(as.numeric(input[t+1,6]))-log(as.numeric(input[t,6])))
}

#construct a data set with adjclose values as well as log returns
numOfWeeks = nrow(N225)
returnsDataSet <- data.frame(matrix(ncol=3, nrow=0, dimnames=list(NULL, c("Week", "n225_returns", "ftse_returns"))))

for (i in 1:numOfWeeks-1){
  returnsDataSet[i,1]=i
  returnsDataSet[i,2]=returns(FTSE,i)
  returnsDataSet[i,3]=returns(N225,i)

 
}
returnsDataSet <- na.omit(returnsDataSet)

print(head(returnsDataSet))



#the identification step to find the order p 
print(jarque.bera.test(returnsDataSet$n225_returns))

par(mfrow=c(3,2))
plot(returnsDataSet$Week, returnsDataSet$n225_returns, type= "l")
plot(N225$N225.Adjusted, type="l")
acf(returnsDataSet$n225_returns, col="green", lwd=2)
pacf(returnsDataSet$n225_returns, col="green", lwd=2)
acf(returnsDataSet$n225_returns^2, col="red", lwd=2)
pacf(returnsDataSet$n225_returns^2, col="red", lwd=2)
#so we will use 3 lags.
par(mfrow=c(3,2))
plot(returnsDataSet$Week, returnsDataSet$ftse_returns, type= "l")
plot(FTSE$FTSE.Adjusted, type="l")
acf(returnsDataSet$ftse_returns, col="green", lwd=2)
pacf(returnsDataSet$ftse_returns, col="green", lwd=2)
acf(returnsDataSet$ftse_returns^2, col="red", lwd=2)
pacf(returnsDataSet$ftse_returns^2, col="red", lwd=2)
#so we will use 3 lags. (?)
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
  model1=garchFit(formula=~arma(3,0)+garch(1,1),data=returnsDataSet$ftse_returns,trace=F,cond.dist=distributions[k])
  cond_dist_table2[k,1]=distributions[k]
  
  cond_dist_table2[k,2]=(summary(model1))$ics["AIC"]
  cond_dist_table2[k,3]=(summary(model1))$ics["BIC"]
}

cond_dist_table2
#we choose also sstd
#we conclude by fitting an ar(3)garch(1,1) model to both data
model1=garchFit(formula=~arma(3,0)+garch(1,1),data=returnsDataSet$n225_returns,trace=F,cond.dist="sstd")
model2=garchFit(formula=~arma(3,0)+garch(1,1),data=returnsDataSet$ftse_returns,trace=F,cond.dist="sstd")





# Step 3: Model checking
res1 <- residuals(model1, standardize=TRUE)
par(mfrow=c(2,2))
acf(res1, col="green", lwd=2)
acf(res1^2, col="red", lwd=2)
pacf(res1, col="green", lwd=2)
pacf(res1^2, col="red", lwd=2)
par(mfrow=c(1,1))
Box.test(res1, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res1^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
u1<-pnorm(res1, mean=0, sd=1)[4:numOfWeeks-4]

head(rev(u1))
res2 <- residuals(model2, standardize=TRUE)
par(mfrow=c(2,2))
acf(res2, col="green", lwd=2)
acf(res2^2, col="red", lwd=2)
pacf(res2, col="green", lwd=2)
pacf(res2^2, col="red", lwd=2)
par(mfrow=c(1,1))
Box.test(res2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res2^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)

#suggests we have captured autocorrelation

u2<-pnorm(res2, mean=0, sd=1)[4:numOfWeeks-4]
#needs checking
par(mfrow=c(2,1))
hist(u1)
hist(u2)

#generating epsilon using copulas
copModel=BiCopSelect(u1, u2, familyset=NA, selectioncrit="AIC", indeptest=TRUE, level=0.05,se = TRUE)
copModel
copModel2=BiCopSelect(u1, u2, familyset=NA, selectioncrit="BIC", indeptest=TRUE, level=0.05,se = TRUE)
copModel2

#the best model to use is the Bivariate copula: t (par = 0.6, par2 = 17.31, tau = 0.41) so we will use students t
N=10000
varEstimate1 <- slot(model1, "sigma.t")
varEstimate1
u_sim=BiCopSim(N, family=copModel$family, copModel$par,  copModel$par2)
res1_sim=qnorm(u_sim[,1], mean = 0, sd = 1) 
res2_sim=qnorm(u_sim[,2], mean = 0, sd = 1) 
res1_sim
print(tail(u1, n=1))
#introduce GARCH effects
sigmaSquaredHat1=model1@fit$coef["omega"] + model1@fit$coef["alpha1"] * (tail(u1, n=1))^2+ model1@fit$coef["beta1"] * (tail(varEstimate1, n=1))^2
tail=as.numeric(tail(returnsDataSet$n225_returns, n=3))

sim_values1 = (model1@fit$coef["mu"]+model1@fit$coef["ar1"]*tail[3] +model1@fit$coef["ar2"]*tail[2] +model1@fit$coef["ar3"]*tail[1])+sigmaSquaredHat1*res1_sim


#and for the other model
varEstimate2 <- slot(model2, "sigma.t")
varEstimate2
sigmaSquaredHat2=model2@fit$coef["omega"] + model2@fit$coef["alpha1"] * (tail(u2, n=1))^2+ model2@fit$coef["beta1"] * (tail(varEstimate2, n=1))^2
tail2=as.numeric(tail(returnsDataSet$ftse_returns, n=3))

sim_values2 = (model2@fit$coef["mu"]+model2@fit$coef["ar1"]*tail2[3] +model2@fit$coef["ar2"]*tail2[2] +model2@fit$coef["ar3"]*tail2[1])+sigmaSquaredHat2*res2_sim
sim_values1
sim_values2
#find var
hist(sim_values1)
#this looks like a normal distribution  so lets fit one
hist(sim_values2)
# Fit normal distribution to log-transformed values


fit1 <- fitdistr(sim_values1, densfun = "normal")
fit2 <- fitdistr(sim_values_2, densfun = "normal")
portfolio_mean=0.5*(fit1$estimate[1]+fit2$estimate[1])
portfolio_sd=sqrt(0.25*((fit1$estimate[2])^2+(fit2$estimate[2])^2)+0.5*(fit1$estimate[2])*(fit2$estimate[2])*cor(returnsDataSet$n225_returns, returnsDataSet$ftse_returns))

VaR99=2.326*portfolio_sd-portfolio_mean
VaR95=1.645*portfolio_sd-portfolio_mean

VaR99
VaR95