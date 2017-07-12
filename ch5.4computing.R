# ch5.4 analyses

#time series (hourly) over two years of neutron levels both fast and thermal
#associated time series for humidity, soil moisture, and rainfall
#we do NOT have measurements of biomass, even though neutron levels are associated
#with humidity, soil moisture, and biomass with jumps caused by rainfall

#modeling Water Vapor Corrected Level 2 Neutron Counts (cph)
#the idea is that we predict N_(t+1) using N(1,...,t), input \hat{N_(t+1)} to
#the equation using a fixed value for N_0, and get \hat{\Theta_(t+1)}. Then compare
#the estimated soil moisture at (t+1) with the observed soil moisture at (t+1).

#The research goal (overall) is to deconvolve neutron levels to get separate component
#time series for humidity, soil moisture, rainfall, and biomass (so as to solve for
#biomass). 

# classical detrend/deseasonal, ARIMA, bayesTS

rm(list=ls())
library(TTR)
data<-read.csv('STER_daily_obs.csv', header=T, row.names=NULL)
# log series after bring min close to zero,, try to reduce spikes
WVClog<-log10(data$WVC_FNC-2100)
WVClogts<-ts(WVClog, frequency=365, start=c(2011, 180))
plot(ts(WVClogts))
# decomposition assuming additivity (moving average)
WVCtscomp<-decompose(WVClogts)
plot(WVCtscomp, cex=1.4, cex.lab=1.4, cex.axis=1.4)
# save as decompose_v2.pdf

# differencing
WVCtsdiff1<-diff(WVClogts, differences=1)
WVCtsdiff2<-diff(WVCtsdiff1, differences=1)
WVCtsdiff3<-diff(WVCtsdiff2, differences=1)
WVCtsdiff4<-diff(WVCtsdiff3, differences=1)
par(mfrow=c(2,2))
plot.ts(WVCtsdiff1,ylim=c(-2.8,1.8))
plot.ts(WVCtsdiff2,ylim=c(-2.8,1.8))
plot.ts(WVCtsdiff3,ylim=c(-2.8,1.8))
plot.ts(WVCtsdiff4,ylim=c(-2.8,1.8))

# smooth by MA
WVCts<-ts(data$WVC_FNC, frequency=365, start=c(2011, 180))
WVCsma20<-SMA(WVCts, n=20)
# decomposition assuming additivity (moving average)
WVCtsma20<-ts(WVCsma20, frequency=365, start=c(2011,180))
WVCts20comp<-decompose(WVCtsma20)
# remove seasonal
WVCsma20sadj<-WVCtsma20-WVCts20comp$seasonal
# extend trend to cover entire series
WVCtr<-c(WVCts20comp$trend[413:594],WVCts20comp$trend[202:594],WVCts20comp$trend[202:402])
# remove seasonal and trend
WVCsma20adj<-WVCts20comp$seasonal-WVCtr
WVCsmad1<-diff(WVCsma20adj, differences=1)
WVCsmad2<-diff(WVCsmad1, differences=1)
WVCsmad3<-diff(WVCsmad2, differences=1)
# examine second difference
par(mfrow=c(3,1))
plot.ts(WVCsma20adj)
plot.ts(WVCsmad2)
plot(1:90,WVCsmad2[161:250],type='l')

# hmm? take difference of two years, model it as ARIMA
## take 2nd yr - 1st yr to remove season then difference (like SARIMA)
WVCts2<-ts(data$WVC_FNC[1:730],frequency=365, start=c(2011, 1))
WVCts2d<-WVCts2[366:730]-WVCts2[1:365]
WVCts2dd<-diff(WVCts2d, differences=1)
WVCts2dd2<-diff(WVCts2dd, differences=1)
WVCts2dd3<-diff(WVCts2dd2, differences=1)
# note WVCts2dd looks best, others have expanded range
 
# fit ARMA
par(mfrow=c(2,1))
acf(WVCts2dd,lag.max=20, main="", cex.lab=1.2, cex.axis=1.2)
pacf(WVCts2dd,lag.max=20, main="", cex.lab=1.2, cex.axis=1.2)
# save as 2yrdiff_acf_pacf_v2.pdf

# both have significant terms up to lag 5
mset<-vector("list", 36)
maic<-c(rep(0,36))
for(i in 1:6){
for(j in 1:6){
mset[[6*(i-1)+j]]<-arima(WVCts2dd, order=c(i,0,j), method="ML", include.mean=FALSE)
maic[6*(i-1)+j]<-mset[[6*(i-1)+j]]$aic
}
}
# maic chooses 16th model or ARMA(3,4)
# shortcut using auto.arima() assume non-seasonal (using AICc)
library(forecast)
mauto<-auto.arima(WVCts2dd,max.order=20, stepwise=FALSE)
# chooses arima(3,0,5)

## assess ARMA(3,4) and ARMA(3,5)
m1<-arima(WVCts2dd,order=c(3,0,4), include.mean=FALSE)
m2<-arima(WVCts2dd,order=c(3,0,5), include.mean=FALSE)
Box.test(m1$residuals, type="Ljung")
	Box-Ljung test
data:  m1$residuals
X-squared = 9e-04, df = 1, p-value = 0.976
Box.test(m2$residuals, type="Ljung")
	Box-Ljung test
data:  m2$residuals
X-squared = 1e-04, df = 1, p-value = 0.9937
plotForecastErrors2 <- function(forecasterrors,main)
   {
      # make a histogram of the forecast errors:
      mybinsize <- IQR(forecasterrors)/4
      mysd   <- sd(forecasterrors)
      mymin  <- min(forecasterrors) - mysd*5
      mymax  <- max(forecasterrors) + mysd*3
      # generate normally distributed data with mean 0 and standard deviation mysd
      mynorm <- rnorm(10000, mean=0, sd=mysd)
      mymin2 <- min(mynorm)
      mymax2 <- max(mynorm)
      if (mymin2 < mymin) { mymin <- mymin2 }
      if (mymax2 > mymax) { mymax <- mymax2 }
      # make a bw histogram of the forecast errors, with the normally distributed data overlaid:
      mybins <- seq(mymin, mymax, mybinsize)
      hist(forecasterrors, col="grey90", freq=FALSE, breaks=mybins, main=main, xlab="forecast errors")
      # freq=FALSE ensures the area under the histogram = 1
      # generate normally distributed data with mean 0 and standard deviation mysd
      myhist <- hist(mynorm, plot=FALSE, breaks=mybins, xlab="forecast errors")
      # plot the normal curve as a black line on top of the histogram of forecast errors:
      points(myhist$mids, myhist$density, type="l", col="black", lwd=2)
   }
par(mfrow=c(2,1))
plotForecastErrors2(m1$residuals, main="residuals from ARMA(3,4)")
plotForecastErrors2(m2$residuals, main="residuals from ARMA(3,5)")
# save as forecast_errors_SARIMA_v2.pdf

## classical decomposition, several ways
WVClog<-log10(data$WVC_FNC-2100)
WVClogts<-ts(WVClog, frequency=365, start=c(2011, 180))
WVCtscomp<-decompose(WVClogts)
plot(WVCtscomp)
WVCtscomp2<-stl(WVClogts,s.window="periodic")
plot(WVCtscomp2, set.pars=list(mar = c(0, 6, 0, 6), oma = c(6, 0, 4, 0), tck = -0.01, mfrow = c(4, 1), cex.lab=1.4, cex.axis=1.4, cex=0.8))
# save as decompose_stl_v2.pdf
acf(WVCtscomp2$time.series[,3], main="", cex.lab=1.2, cex.axis=1.2)
# save as decompose_stl_acfresids_v2.pdf

## using auto.arima and years 1 and 2 of data to start, do one-step forecasting with 
## data updating
## the auto.arima output is the same as for arima, so you can use $arma to get a 	
## compact form of the specification, as a vector giving the number of AR, MA, 
## seasonal AR and seasonal MA coefficients, plus the period and the number of non-seasonal 
## and seasonal differences.

library(forecast)
# save predictive errors on original scale
ferr<-list(aicc=c(rep(99,length=46)),aic=c(rep(99,length=46)),bic=c(rep(99,length=46)))
# save predicted values on differenced scale
fval<-list(aicc=vector("list",46),aic=vector("list",46),bic=vector("list",46))
# save predicted values and 95% CI on original scale
fvalr<-list(aicc=c(rep(99,length=46)),aic=c(rep(99,length=46)),bic=c(rep(99,length=46)))
fvalrl<-list(aicc=c(rep(99,length=46)),aic=c(rep(99,length=46)),bic=c(rep(99,length=46)))
fvalru<-list(aicc=c(rep(99,length=46)),aic=c(rep(99,length=46)),bic=c(rep(99,length=46)))
# save models chosen
fmod<-list(aicc=vector("list",46),aic=vector("list",46),bic=vector("list",46))
WVCtsa2<-ts(data$WVC_FNC,frequency=365, start=c(2011, 1))
WVCtsa2d<-c(WVCtsa2[366:730]-WVCtsa2[1:365],WVCtsa2[731:776]-WVCtsa2[366:411])
WVCtsa2dd<-diff(WVCtsa2d, differences=1)
for(i in 1:46){	
	WVCtst2<-ts(data$WVC_FNC[1:(730+i-1)],frequency=365, start=c(2011, 1))
	if(i>1){
			WVCtst2d<-c(WVCtst2[366:730]-WVCtst2[1:365],WVCtst2[731:(730+i-1)]-WVCtst2[366:(365+i-1)])
		} else {
			WVCtst2d<-c(WVCtst2[366:730]-WVCtst2[1:365])
			}
	WVCtst2dd<-diff(WVCtst2d, differences=1)
	# find best ARMA model
	fmod$aicc[[i]]<-auto.arima(WVCtst2dd,max.order=20, stepwise=FALSE)
	fmod$aic[[i]]<-auto.arima(WVCtst2dd,max.order=20, stepwise=FALSE, ic="aic")
	fmod$bic[[i]]<-auto.arima(WVCtst2dd,max.order=20, stepwise=FALSE, ic="bic")
	# get one-step prediction
	fval$aicc[[i]]<-forecast(fmod$aicc[[i]],1)
	fval$aic[[i]]<-forecast(fmod$aic[[i]],1)
	fval$bic[[i]]<-forecast(fmod$bic[[i]],1)
	## get one-step predictive errors on ORIGINAL scale!
	fvalr$aicc[i]<-WVCtsa2[730+i-1]-WVCtsa2[365+i-1]+WVCtsa2[365+i]+fval$aicc[[i]]$mean
	fvalrl$aicc[i]<-WVCtsa2[730+i-1]-WVCtsa2[365+i-1]+WVCtsa2[365+i]+fval$aicc[[i]]$lower[1,2]
	fvalru$aicc[i]<-WVCtsa2[730+i-1]-WVCtsa2[365+i-1]+WVCtsa2[365+i]+fval$aicc[[i]]$upper[1,2]
	ferr$aicc[i]<-fvalr$aicc[i]-WVCtsa2[365+i]
	fvalr$aic[i]<-WVCtsa2[730+i-1]-WVCtsa2[365+i-1]+WVCtsa2[365+i]+fval$aic[[i]]$mean
	fvalrl$aic[i]<-WVCtsa2[730+i-1]-WVCtsa2[365+i-1]+WVCtsa2[365+i]+fval$aic[[i]]$lower[1,2]
	fvalru$aic[i]<-WVCtsa2[730+i-1]-WVCtsa2[365+i-1]+WVCtsa2[365+i]+fval$aic[[i]]$upper[1,2]
	ferr$aic[i]<-fvalr$aic[i]-WVCtsa2[730+i]
	fvalr$bic[i]<-WVCtsa2[730+i-1]-WVCtsa2[365+i-1]+WVCtsa2[365+i]+fval$bic[[i]]$mean
	fvalrl$bic[i]<-WVCtsa2[730+i-1]-WVCtsa2[365+i-1]+WVCtsa2[365+i]+fval$bic[[i]]$lower[1,2]
	fvalru$bic[i]<-WVCtsa2[730+i-1]-WVCtsa2[365+i-1]+WVCtsa2[365+i]+fval$bic[[i]]$upper[1,2]
	ferr$bic[i]<-fvalr$bic[i]-WVCtsa2[730+i]
}
# aic and aicc results same
par(mfrow=c(2,1))
plot(c(731:776),WVCtsa2[731:776],type='l',col='black')
lines(c(731:776),fvalr$aic,col='red')
lines(c(731:776),fvalr$bic,col='green')
plot(c(731:776),ferr$aic,type='l',col='red')
lines(c(731:776),ferr$bic,col='green')
#
par(mfrow=c(2,1))
plot(c(731:776),WVCtsa2[731:776],type='l',col='black', xlab='time point (days)', ylab='prediction', cex.lab=1.2, cex.axis=1.2)
lines(c(731:776),fvalr$bic,lty=2)
plot(c(731:776),WVCtsa2[731:776],type='l',col='black', xlab='time point (days)', ylab='prediction', cex.lab=1.2, cex.axis=1.2)
lines(c(731:776),fvalrl$bic,lty=2)
lines(c(731:776),fvalru$bic,lty=2)
# save as sarima_freq_results2_v2.pdf

# chosen model
ar<-list(aic=rep(99,46),aicc=rep(99,46),bic=rep(99,46))
ma<-ar
for(i in 1:46){
	ar$aic[i]<-fmod$aic[[i]]$arma[1]
	ar$aicc[i]<-fmod$aicc[[i]]$arma[1]
	ar$bic[i]<-fmod$bic[[i]]$arma[1]
	ma$aic[i]<-fmod$aic[[i]]$arma[2]
	ma$aicc[i]<-fmod$aicc[[i]]$arma[2]
	ma$bic[i]<-fmod$bic[[i]]$arma[2]
}
plot(c(731:776),ar$aic,type='l',xlab="time point (days)", ylab='order', cex.lab=1.2, cex.axis=1.2)
lines(c(731:776),ma$aic,lty=2)
lines(c(731:776),ar$bic,lty=3)
lines(c(731:776),ma$bic,lty=4)
legend('topright',legend=c('aic ar','aic ma','bic ar','bic ma'),lty=c(1,2,3,4), cex=1.2)
# save as sarima_freq_models_v2.pdf

# using ALL of the data find out what ARMA is selected by hand and by auto.arima
mset2<-vector("list", 36)
maic2<-c(rep(0,36))
for(i in 1:6){
	for(j in 1:6){
		mset2[[6*(i-1)+j]]<-arima(WVCtsa2dd, order=c(i,0,j), method="ML", include.mean=FALSE)
		maic2[6*(i-1)+j]<-mset2[[6*(i-1)+j]]$aic
	}
}
# maic2 chooses 36th model or ARMA(6,6)
# shortcut using auto.arima() assume non-seasonal (using AICc)
library(forecast)
mauto2<-auto.arima(WVCtsa2dd,max.order=20, stepwise=FALSE)
# ARIMA(1,0,2) with zero mean  

## assess ARMA(6,6) and ARMA(1,2)
m1a<-arima(WVCts2dd,order=c(6,0,6), include.mean=FALSE)
m2a<-arima(WVCts2dd,order=c(1,0,2), include.mean=FALSE)
Box.test(m1a$residuals, type="Ljung")
	Box-Ljung test
data:  m1a$residuals
X-squared = 0.0011, df = 1, p-value = 0.9741
Box.test(m2a$residuals, type="Ljung")
	Box-Ljung test
data:  m2$residuals
X-squared = 0.0062, df = 1, p-value = 0.9375
par(mfrow=c(2,1))
plotForecastErrors2(m1a$residuals, main="resids from ARMA(6,6)")
plotForecastErrors2(m2a$residuals, main="resids from ARMA(1,2)")
# save as forecast_errors_SARIMA2_v2.pdf

## Bayes ARMA (?) on the same data
library(dlm)
# uncover ARMA in dlm form
tmp<-dlmModARMA(ar=c(0.3,0.4,0.5),ma=c(0.6,0.7))
tmp
$FF
     [,1] [,2] [,3]
[1,]    1    0    0
$V
     [,1]
[1,]    0
$GG
     [,1] [,2] [,3]
[1,]  0.3    1    0
[2,]  0.4    0    1
[3,]  0.5    0    0
$W
     [,1] [,2] [,3]
[1,]  1.0 0.60 0.70
[2,]  0.6 0.36 0.42
[3,]  0.7 0.42 0.49
$m0
[1] 0 0 0
$C0
      [,1]  [,2]  [,3]
[1,] 1e+07 0e+00 0e+00
[2,] 0e+00 1e+07 0e+00
[3,] 0e+00 0e+00 1e+07
# specify same model in dlm form
# recall y_t = FF\theta_t + v_t where v_t~N(0,V); \theta_t = GG\theta_(t-1)+w_t where w_t~N(0,W)
# FF is coefficient of theta, V is observed error variance, GG is innovation coefficients, W is innovation variance, m0 is pre-sample state, C0 is pre-sample variance, 
# length of FF and dim(GG) and dim(W) and dim(C0) and length(m0) is max(p,q) if p bigger and max(p,q)+1 if q bigger
myModel<-dlm(FF=matrix(data=c(1,0,0),nrow=1,ncol=3),V=0,GG=matrix(data=c(0.3,0.4,0.5,1,0,0,0,1,0),nrow=3,ncol=3),W=matrix(data=c(1,0.6,0.7,0.6,0.36,0.42,0.7,0.42,0.49),nrow=3,ncol=3),m0=c(0,0,0),C0=diag(x=1e+07,3,3))

## OK ARMA is not easily done, so let's just use existing facilities in DLM package to fit 
# try Gibbs sampling to fit approximation to ARMA(1,2) with AR coeff fixed at 1
outGibbs<-dlmGibbsDIG(WVCts2dd, dlmModPoly(2), a.y=1000, b.y=1, shape.theta=2, rate.theta=2, n.sample=5500, save.states=FALSE)
  |=====================================================| 100%
burn<-1000
attach(outGibbs)
dV<-dV[-(1:burn)]
dW<-dW[-(1:burn),]
detach()
par(mfrow=c(2,3), mar=c(3.1, 2.1, 2.1, 1.1))
plot(dV, type='l', xlab="", ylab="", main=expression(sigma^2))
plot(dW[,1], type='l', xlab="", ylab="", main=expression(sigma[beta1]^2))
plot(dW[,2], type='l', xlab="", ylab="", main=expression(sigma[beta2]^2))
use<-length(dV)-burn
from<-0.05*use
at<-pretty(c(0,use), n=3); at<-at[at >= from]
plot(ergMean(dV, from), type='l', xaxt='n',xlab="", ylab="")
axis(1, at=at-from, labels=format(at))
plot(ergMean(dW[,1], from), type='l', xaxt='n', xlab="", ylab="")
axis(1, at=at-from, labels=format(at))
plot(ergMean(dW[,2], from), type='l', xaxt='n', xlab="", ylab="")
axis(1, at=at-from, labels=format(at))
mcmcMean(cbind(dV[-(1:burn)], dW[-(1:burn), ]))
# look at theoretical and prediction bounds
n<-20
m<-1; p<-2
new<-100
fore<-dlmForecast(m, nAhead=n, sampleNew=new)
ciTheory <- (outer(sapply(fore$Q, FUN=function(x) sqrt(diag(x))), qnorm(c(0.1,0.9))) +as.vector(t(fore$f)))
ciSample <- t(apply(array(unlist(fore$newObs), dim=c(n,m,new))[,1,], 1,FUN=function(x) quantile(x, c(0.1,0.9))))
plot.ts(cbind(ciTheory,fore$f[,1]),plot.type="s", col=c("red","red","green"),ylab="y")
for (j in 1:2) lines(ciSample[,j], col="blue")
legend(2,-40,legend=c("forecast mean", "theoretical bounds", "Monte Carlo bounds"),col=c("green","red","blue"), lty=1, bty="n")

## implement bayesian ARMA ala Monahan
library(Matrix)
# specify prior values
## assume gamma=zero, gamma^*=zero, tau=\infty,(tau^*)^(-1)=zero, mu=zero (as data is centered)
fprior<-list(alpha=10, beta=10/0.01, gamma=1)
#fprior<-list(alpha=10, beta=10/10e-4, gamma=1)
##fprior<-list(alpha=seq(10,40,10), beta=seq(10,40,10)/1e-04, gamma=seq(-1,2,1))
nmod<-1
##nmod<-64
# save predictive errors on original scale
ferr<-matrix(99,nmod,46)
# save predicted values on differenced scale
fval<-matrix(99,nmod,46)
# save predicted values and 95% CI on original scale
fvalr<-matrix(99,nmod,46)
fvalrl<-matrix(99,nmod,46)
fvalru<-matrix(99,nmod,46)
# save models chosen
fmod<-matrix(99,46,2)
WVCtsa2<-ts(data$WVC_FNC,frequency=365, start=c(2011, 1))
WVCtsa2d<-c(WVCtsa2[366:730]-WVCtsa2[1:365],WVCtsa2[731:776]-WVCtsa2[366:411])
WVCtsa2dd<-diff(WVCtsa2d, differences=1)
WVCtsa2dd<-WVCtsa2dd-mean(WVCtsa2dd)
# function for location/scale parameterized t distribution; thx to Christian Groll
qt_ls <- function(prob, df, mu, a) qt(prob, df)*a + mu
for(i in 1:46){	
	WVCtst2<-ts(data$WVC_FNC[1:(730+i-1)],frequency=365, start=c(2011, 1))
	if(i>1){
			WVCtst2d<-c(WVCtst2[366:730]-WVCtst2[1:365],WVCtst2[731:(730+i-1)]-WVCtst2[366:(365+i-1)])
		} else {
			WVCtst2d<-c(WVCtst2[366:730]-WVCtst2[1:365])
			}
	WVCtst2dd<-diff(WVCtst2d, differences=1)
	WVCtst2dd<-WVCtst2dd-mean(WVCtst2dd)
	nob<-length(WVCtst2dd)
	# calculate A_T, A_l, and A_21 using the acf covariance
	tmp<-acf(WVCtst2dd, lag.max=(nob-1), type="covariance", plot=FALSE)
	diags<-vector("list", length=nob)
	for(j in 1:(nob-1)){
		diags[[j]]<-rep(tmp[[1]][j],nob-j+1)
	}
	diags[[nob]]<-0
	AT<-bandSparse(nob, k=-c(0:(nob-1)), diag=diags, symm=TRUE)
	Al<-var(WVCtst2dd)
	A12<-matrix(rev(tmp[[1]]),nob,1)
	# compute location and scale of predictive distribution
	MVSloc<-0+t(A12)%*%chol2inv(chol(AT))%*%matrix(WVCtst2dd,nob,1)
	MVSloc<-MVSloc[1,1]
	betaST<-fprior$beta+matrix(WVCtst2dd,1,nob)%*%chol2inv(chol(AT))%*%matrix(WVCtst2dd,nob,1)
	betaST<-betaST[1,1]
	ATl<-(Al-t(A12)%*%chol2inv(chol(AT))%*%matrix(WVCtst2dd,nob,1))
	ATl<-ATl[1,1]
	MVSscale<-(2*fprior$alpha + (nob/2)*betaST)/ATl
	MVSdf<-2*fprior$alpha+nob
	fmod[i,]<-c(MVSloc, MVSscale)
	# get one-step prediction
	fval[nmod,i]<-MVSloc
	## get one-step predictive errors on ORIGINAL scale!
	fvalr[nmod,i]<-WVCtsa2[730+i-1]-WVCtsa2[365+i-1]+WVCtsa2[365+i]+(fval[,i]+mean(WVCtst2dd))
	fvalrl[nmod,i]<-WVCtsa2[730+i-1]-WVCtsa2[365+i-1]+WVCtsa2[365+i]+qt_ls(0.05,MVSdf,MVSloc,MVSscale)+mean(WVCtst2dd)
	fvalru[nmod,i]<-WVCtsa2[730+i-1]-WVCtsa2[365+i-1]+WVCtsa2[365+i]+qt_ls(0.95,MVSdf,MVSloc,MVSscale)+mean(WVCtst2dd)
	ferr[nmod,i]<-fvalr[nmod,i]-WVCtsa2[365+i]
}
par(mfrow=c(2,1))
plot(1:46,WVCtsa2[731:776],type='l',ylim=c(1700,3760), xlab='predicted observations (1 to 46)', ylab='response')
lines(1:46,fvalr[nmod,],lty=2)
lines(1:46,fvalrl[nmod,],lty=3)
lines(1:46,fvalru[nmod,],lty=3)
plot(1:46,ferr[nmod,],type='l',lty=4, xlab='predicted observations (1 to 46)', ylab='predictive errors')
# save as bayesARMApriorsm.pdf (for prior beta=10/0.01) or bayesARMApriorbig.pdf (for prior beta=10/1e-04)

