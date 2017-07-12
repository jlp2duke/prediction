# ch5.5 notes for stochastic modeling
#	NOAA climate data
# http://www.ncdc.noaa.gov/cdo-web/

# suppose data is year, # rainfall events per year, total rainfall per year
			1				y_{1}=Y(1)						Z(1)
			2				y_{2}=Y(2)-Y(1)					Z(2)
			T				y_{t}=Y(T)-Y(T-1)				Z(T)
#							poisson process				compound poisson process
# \bar{a} = average rainfall per event = (\sum_1^T(Z)/\sum_1^T(Y(T)-Y(T-1)))
#	Y(T+1)-Y(T) = # rainfall events in year T+1 ~ Poisson(\lambda) 
#	where \hat{\lambda} = sum_1^T(Y)/T 
#		and prediction for Y(T+1) is \hat{Y(T+1)}=\bar{a}[Y(T+1)-Y(T)]~\bar{a}Poisson(\hat{\lamda})
#	This is a predictive distribution. 
#	For evaluation use two different methods
#	(1) -\log P_\hat{\lambda}(y_{T+1}) where y_{T+1} = Y(T+1)-Y(T)
# 	Over k-ahead forecasts -\log P_\hat{\lambda}(y_{T+1}),\ldots,-\log P_\hat{\lambda}(y_{T+k})
#		we use log score = LS = -(1/k)\sum_{j=1}^{k} \log P_\hat{\lambda}(y_{T+j})
#	(2) PITS (probability integral transform score):	
#			Find P_\hat{\lambda}(y_{T+1}),\ldots,P_\hat{\lambda}(y_{T+k})
#			Form a histogram. If predictor is 'good' histogram ~ uniform
#			Assess PITScore = L1 distance of (EDF(histogram) - uniform) 
#	Need k-ahead forecasts because 1-ahead forecasts are unstable (like LOO CV)
#
#	Z ~ Erlang(Y(T)-Y(T-1), \Theta) as Z(T)=X_1,T + \cdots + X_Y(T)-Y(T-1),T where
#	X_i,T is the ith rainfall event in year T
#	Y(T)-Y(T-1)~Poisson(\lambda) and X_i,T~exponential(\Theta)
#	using method of moments to get estimates yields 
#		\hat(\Theta)=s^2_Z/2\bar(Z)
#		and \hat(\lambda)=2(\bar(Z)^2)/s^2_Z
#	where s^2_Z = variance(Z) and \bar(Z) = mean(Z)
#	So our probabilistic prediction comes from Erlang as Z_{T+1}~Erlang(Poisson(\hat{\lambda}),\hat{\Theta})
#	We can randomly draw from Poisson(\hat{\lambda}) 
#	Then log score -\log P_{\hat{\lambda},\hat{\Theta}}(\hat_{Z_{T+1}}),\ldots,-\log P_{\hat{\lambda},\hat{\Theta}}(\hat_{Z_{T+k}})
#	PITS based on P_{\hat{\lambda},\hat{\Theta}}(\hat_{Z_{T+1}}),\ldots,P_{\hat{\lambda},\hat{\Theta}}(\hat_{Z_{T+k}})

dat<-read.csv('NCDC_237581_edited.csv', header=T, row.names=NULL)
dim(dat)
#[1] 696   5
nyear<-dim(dat)[1]/12
#[1] 58
# data is by month so consolidate to year level
datnew<-list(nevent=rep(99,nyear),totr=rep(99,nyear))
for(i in 1:nyear){
	datnew$nevent[i]<-sum(dat[((i-1)*12+1):(i*12),2])+sum(dat[((i-1)*12+1):(i*12),3])+sum(dat[((i-1)*12+1):(i*12),4])
	datnew$totr[i]<-sum(dat[((i-1)*12+1):(i*12),5])
}
summary(datnew$nevent)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  20.00   30.00   34.00   35.48   41.00   63.00 
summary(datnew$totr)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  567.0   970.5  1158.0  1144.0  1292.0  2186.0 
## use first 30 years for modeling and the remaining for prediction
nfit<-30

## number of rainfall events
bara<-sum(datnew$totr[1:nfit])/sum(datnew$nevent[1:nfit])
hatlambda1<-sum(datnew$nevent[1:nfit])/nfit
# compute log score
LS1all<-(-dpois(datnew$nevent[(nfit+1):nyear], hatlambda, log = TRUE))
LS1<-sum(LS1all)/(nyear-nfit)
library(vioplot)
plot(1,5,type="n", xlim=c(0.5,1.5), ylim=c(2.7,6.7), xlab="log score", xaxt="n",ylab="",cex.lab=1.2, cex.axis=1.2)
vioplot(LS1all, col="lightgrey",names='log score', lwd=2, add=T)
# save as logscore_nevents_v2.pdf

# compute PITS
EDF1<-ppois(datnew$nevent[(nfit+1):nyear],hatlambda, lower.tail=TRUE, log.p=FALSE)
UNI1<-punif(datnew$nevent[(nfit+1):nyear],min=min(datnew$nevent[1:nfit]), max=max(datnew$nevent[1:nfit]), lower.tail=TRUE, log.p=FALSE)
PITS1<-sum(abs(EDF1-UNI1))
plot(datnew$nevent[(nfit+1):nyear],EDF1,pch=2, xlab='number of rainfalls (annual)', ylab='CDF', ylim=c(0,1.0), cex.lab=1.2, cex.axis=1.2, lwd=2, cex=1.2)
points(datnew$nevent[(nfit+1):nyear],UNI1,pch=1,lwd=2,cex=1.2)
# save as neventsEDF1_v2.pdf

datapit1<-list(UNI=UNI, EDF=EDF1)
library(beeswarm)
beeswarm(datapit1, col='black', pch=c(1,2), main="",cex.axis=1.2, cex=1.2, lwd=2)
# save as neventsPITscore1_v2.pdf

## amount of rainfall
hatTheta<-var(datnew$totr[1:nfit])/(2*mean(datnew$totr[1:nfit]))
hatlambda2<-(2*(mean(datnew$totr[1:nfit]))^2)/var(datnew$totr[1:nfit])
set.seed(10101)
nreps<-1000
LS2<-rep(99, length=nreps)	
shape<-rpois(nreps,hatlambda2)
for(i in 1:nreps){
	tmp<-(-dgamma(datnew$totr[(nfit+1):nyear],shape=shape[i], scale=hatTheta, log = TRUE))
	LS2[i]<-sum(tmp)/(nyear-nfit)
}
library(vioplot)
plot(1,13,type="n", xlim=c(0.5,1.5), ylim=c(7,17.5), xlab="log score", xaxt="n",ylab="",cex.lab=1.2, cex.axis=1.2)
vioplot(LS2, col="lightgrey",names='log score', lwd=2, add=T)
# save as logscore_totrain.pdf

PITS2<-rep(99, length=nreps)	
UNI2<-punif(datnew$totr[(nfit+1):nyear],min=min(datnew$totr[1:nfit]), max=max(datnew$totr[1:nfit]), lower.tail=TRUE, log.p=FALSE)
EDF2<-matrix(99,nreps,nyear-nfit)
for(i in 1:nreps){
	EDF2[i,]<-pgamma(datnew$totr[(nfit+1):nyear], shape=shape[i], scale=hatTheta, lower.tail=TRUE, log.p=FALSE)
	PITS2[i]<-sum(abs(EDF2[i,]-UNI2))
}
sum(abs(apply(EDF2,2,mean)-UNI2))
#[1] 4.088643

plot(datnew$totr[(nfit+1):nyear],apply(EDF2,2,mean),pch=2, xlab='total rainfall (0.01 inches/year)', ylab='CDF', ylim=c(0,1.0), cex.lab=1.2, cex.axis=1.2, lwd=2, cex=1.2)
points(datnew$totr[(nfit+1):nyear],UNI2,pch=1,lwd=2,cex=1.2)
# save as totrainEDF2_v2.pdf

plot(datnew$totr[(nfit+1):nyear],apply(EDF2,2,mean),pch=2, xlab='total rainfall (0.01 inches/year)', ylab='CDF', ylim=c(0,1.0), cex.lab=1.2, cex.axis=1.2, lwd=2, cex=1.2)
sel<-seq(100,900,200)
for(i in sel){
	tmp<-sort(datnew$totr[(nfit+1):nyear],index.return=T)
	lines(tmp$x,EDF2[i,tmp$ix],lwd=2,cex=1.2)
}
# save as totrainEDF2ex_v2.pdf

datapit2<-list(UNI=UNI2, EDF=apply(EDF2,2,mean))
library(beeswarm)
beeswarm(datapit2, col='black', pch=c(1,2), main="", ylim=c(0,1.0),cex.axis=1.2, cex=1.2, lwd=2)
# save as totrainPITscore2.pdf



