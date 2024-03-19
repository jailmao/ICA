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
  returnsDataSet[i,4]=as.numeric(n225[i+1,6])
  returnsDataSet[i,5]=as.numeric(fchi[i+1,6])
 
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
#so we will use 3 lags. (?)

modeln225=garchFit(formula=~arma(3,0)+garch(1,1),data=returnsDataSet$n225_returns,trace=F,cond.dist="norm")
modelfchi=garchFit(formula=~arma(3,0)+garch(1,1),data=returnsDataSet$fchi_returns,trace=F,cond.dist="norm")
modeln225


# Step 3: Model checking
res1 <- residuals(modeln225, standardize=TRUE)
par(mfrow=c(2,1))
acf(res1, col="green", lwd=2)
acf(res1^2, col="red", lwd=2)
par(mfrow=c(1,1))
Box.test(res1, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res1^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
u1<-pnorm(res1, mean=0, sd=1)[4:numOfWeeks-3]
hist(u1)
print (rev(u1))
res2 <- residuals(modelfchi, standardize=TRUE)
par(mfrow=c(2,1))
acf(res2, col="green", lwd=2)
acf(res2^2, col="red", lwd=2)
par(mfrow=c(1,1))
Box.test(res2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res2^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)



u2<-pnorm(res2, mean=0, sd=1)[4:numOfWeeks-3]

hist(u2)

#generating epsilon using copulas
copModel=BiCopSelect(u1, u2, familyset=NA, selectioncrit="AIC", indeptest=TRUE, level=0.05,se = TRUE)
copModel
#the best model to use is the Bivariate copula: t (par = 0.6, par2 = 17.31, tau = 0.41) so we will use students t
N=10000
varEstimate1 <- slot(modeln225, "sigma.t")
varEstimate1
u_sim=BiCopSim(N, family=copModel$family, copModel$par,  copModel$par2)
res1_sim=qnorm(u_sim[,1], mean = 0, sd = 1) 
res2_sim=qnorm(u_sim[,2], mean = 0, sd = 1) 
res1_sim
print(tail(u1, n=1))
#introduce GARCH effects
sigmaSquaredHat1=modeln225@fit$coef["omega"] + modeln225@fit$coef["alpha1"] * (tail(u1, n=1))^2+ modeln225@fit$coef["beta1"] * (tail(varEstimate1, n=1))^2
tail=as.numeric(tail(returnsDataSet$n225_returns, n=3))

sim_values_n225 = (modeln225@fit$coef["mu"]+modeln225@fit$coef["ar1"]*tail[3] +modeln225@fit$coef["ar2"]*tail[2] +modeln225@fit$coef["ar3"]*tail[1])+sigmaSquaredHat1*res1_sim


#and for the other model
varEstimate2 <- slot(modelfchi, "sigma.t")
varEstimate2
sigmaSquaredHat2=modelfchi@fit$coef["omega"] + modelfchi@fit$coef["alpha1"] * (tail(u2, n=1))^2+ modelfchi@fit$coef["beta1"] * (tail(varEstimate2, n=1))^2
tail2=as.numeric(tail(returnsDataSet$fchi_returns, n=3))

sim_values_fchi = (modelfchi@fit$coef["mu"]+modelfchi@fit$coef["ar1"]*tail2[3] +modelfchi@fit$coef["ar2"]*tail2[2] +modelfchi@fit$coef["ar3"]*tail2[1])+sigmaSquaredHat2*res2_sim
sim_values_n225
sim_values_fchi
#find var
hist(sim_values_n225)
#this looks like a normal distribution  so lets fit one
# Fit normal distribution to log-transformed values


fit1 <- fitdistr(sim_values_n225, densfun = "normal")
fit2 <- fitdistr(sim_values_fchi, densfun = "normal")
portfolio_mean=0.5*(fit1$estimate[1]+fit2$estimate[1])
portfolio_sd=sqrt(0.25*((fit1$estimate[2])^2+(fit2$estimate[2])^2)+0.5*(fit1$estimate[2])*(fit2$estimate[2])*cor(returnsDataSet$n225_returns, returnsDataSet$fchi_returns))

VaR99=2.326*portfolio_sd-portfolio_mean
VaR95=1.645*portfolio_sd-portfolio_mean

VaR99
VaR95