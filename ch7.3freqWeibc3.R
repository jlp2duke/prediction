%		survival data
%	http://www.umass.edu/statdata/statdata/data/uissurv.txt

% do frequentist Weibull regression in pieces to match psbcGroup bayesian data
%	testing psbcGroup prior to implementation
%	cox PH: j permutations of the data order, fit to first 300, then predict ahead in steps of
%	n observations, get predictive errors, then find average and variance of predictive
%	errors over permutations
%	censoring done within permutations
%	for censored data loss function is L1 if pred<cens and 0 if pred>cens
%   three sets of data: original, censor as N(max_obs, var_obs) [draw from censoring for 
%	each obs, if draw<obs censor and no censor otherwise], censor as limit follow-up time

setwd("/Users/jclarke/Documents/predictive/Book/ch7/uissurv")
library(psbcGroup)
library(coda)
library(survival)

datorig<-read.csv('uissurv.csv', header=T, row.names=NULL)
dim(datorig)
#[1] 628  12
sum(datorig$censor)
#[1] 508
datorig$age<-as.numeric(datorig$age)
datorig$agecat<-as.numeric(datorig$age>median(datorig$age))

# implement permute, model, and predict function; use subset of data
nobs<-120
dat<-datorig[sample(1:dim(datorig)[1],nobs),]
nperm<-20
nbase<-40
npred<-nobs-nbase

# variables to hold final results (all runs)
permord<-vector('list', length=nperm)
modPH<-permord
medPH<-permord
perr<-permord
FPE<-permord
for(i in 1:nperm) permord[[i]]<-c(rep(0,nobs))
for(i in 1:nperm) modPH[[i]]<-list("freq"=vector('list', length=npred),"bayes1"=vector('list', length=npred),"bayes2"=vector('list', length=npred))
for(i in 1:nperm) medPH[[i]]<-list("freq"=c(rep(99,npred)),"bayes1"=c(rep(99,npred)),"bayes2"=c(rep(99,npred)))
for(i in 1:nperm) perr[[i]]<-list("freq"=c(rep(99,npred)),"bayes1"=c(rep(99,npred)),"bayes2"=c(rep(99,npred)))
for(i in 1:nperm) FPE[[i]]<-list("freq"=c(rep(99,npred)),"bayes1"=c(rep(99,npred)),"bayes2"=c(rep(99,npred)))
modPHc1<-modPHc2<-modPHc3<-modPHc4<-modPH
medPHc1<-medPHc2<-medPHc3<-medPHc4<-medPH
perrc1<-perrc2<-perrc3<-perrc4<-perr
FPEc1<-FPEc2<-FPEc3<-FPEc4<-FPE

set.seed(31122)
#set.seed(21122)
#set.seed(21112)
#set.seed(21012)
#set.seed(10010)
# censored data 1 and 2 (random)
newcens<-newcens2<-c(rep(0,nobs))
while(sum(newcens[1:nbase])<ceiling(nbase/2)){
	draw<-rnorm(nobs,quantile(dat$time[dat$cens>0],0.75),sqrt(var(dat$time)))
	newcens<-dat$cens*(1-as.numeric(draw<dat$time))
}
while(sum(newcens2[1:nbase])<ceiling(nbase/2)){
	draw2<-rnorm(nobs,quantile(dat$time[dat$cens>0],0.95),sqrt(var(dat$time)))
	newcens2<-dat$cens*(1-as.numeric(draw2<dat$time))
}
# censored data 3 and 4 (limit followup) to give same number censored
limit<-172
lcens<-dat$cens*(1-as.numeric(limit<dat$time))
limit2<-230
lcens2<-dat$cens*(1-as.numeric(limit2<dat$time))
# generate nperm `proper' censoring orders, i.e., KM has a median 
# done to match KM runs
i<-1
while(i<(nperm+1)){
	permord[[i]]<-sample(nobs)
	kk<-1
	jj<-nbase
	while((kk<2)&(jj<(nobs+1))){
		p1<-survfit(Surv(dat$time[permord[[i]]][1:jj],dat$cens[permord[[i]]][1:jj])~1)
		p2<-survfit(Surv(dat$time[permord[[i]]][1:jj],newcens[permord[[i]]][1:jj])~1)
		p3<-survfit(Surv(dat$time[permord[[i]]][1:jj],newcens2[permord[[i]]][1:jj])~1)
		p4<-survfit(Surv(dat$time[permord[[i]]][1:jj],lcens[permord[[i]]][1:jj])~1)
		p5<-survfit(Surv(dat$time[permord[[i]]][1:jj],lcens2[permord[[i]]][1:jj])~1)
		ptmp<-rbind(p1$surv,p2$surv,p3$surv,p4$surv,p5$surv)
		if(sum(apply(ptmp,1,min)<0.5)<5) kk<-2
		jj<-jj+1
	}
	if(jj>nobs) i<-i+1
}

# frequentist PH to match psbcGroup
# NOTE: covariates must be scaled
i<-1
for(i in 1:nperm){
	j<-1
	survDat<-list()
	# original data
	for(j in 1:npred){
		datfit<-dat[permord[[i]][1:(nbase+j-1)],]
		survDat$t<-datfit$time
		survDat$di<-lcens[permord[[i]][1:(nbase+j-1)]]
		x<-cbind(datfit$age,datfit$treat,datfit$site)
		x<-scale(x)
		ctr<-attr(x,"scaled:center")
		scl<-attr(x,"scaled:scale")
		survDat$age<-x[,1]
		survDat$treat<-x[,2]
		survDat$site<-x[,3]
		fitw <- survreg(Surv(t,di) ~ age+treat+site, data=survDat, dist="weibull")
		newobs<-c(dat$age[nbase+j],dat$treat[nbase+j],dat$site[nbase+j])
		newobs<-(newobs-ctr)/scl
		xt<-log(datfit$time)
		ppred2<- exp((xt-(c(1,newobs)%*%fitw$coefficients))*(1/fitw$scale))
		ppred<-exp(-ppred2)
		index<-which(abs(ppred-0.5)==min(abs(ppred-0.5)))	
		my.med<-datfit$time[index][1]
		newc<-lcens[permord[[i]][nbase+j]]
		newt<-dat$time[permord[[i]][nbase+j]]
		if(newc>0){ 
			perr[[i]]$freq[j]<-abs(newt-my.med) 
		} else if (newt >= my.med){
			perr[[i]]$freq[j]<-abs(newt-my.med) 
		} else {
			perr[[i]]$freq[j]<-0
			}
		if(j<2) FPE[[i]]$freq[j]<-perr[[i]]$freq[j]
		if(j>1) FPE[[i]]$freq[j]<-FPE[[i]]$freq[j-1]+perr[[i]]$freq[j]
		modPH[[i]]$freq[[j]]<-fitw
		medPH[[i]]$freq[j]<-my.med
	}
}
save.image('ch7computing_freqWc3_31122.RData')
#save.image('ch7computing_freqWc3_21122.RData')
#save.image('ch7computing_freqWc3_21112.RData')
#save.image('ch7computing_freqWc3_21012.RData')
#save.image('ch7computing_freqWc3_10010.RData')

## integrating results across permutations; load each workspace, save FPE results
# 20 runs
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_freqWc3_10010.RData")
FPE1wc3<-FPE
modPH1wc3<-modPH
medPH1wc3<-medPH
# 20 runs each
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_freqWc3_21012.RData")
FPE2wc3<-FPE
modPH2wc3<-modPH
medPH2wc3<-medPH
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_freqWc3_21112.RData")
FPE3wc3<-FPE
modPH3wc3<-modPH
medPH3wc3<-medPH
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_freqWc3_21122.RData")
FPE4wc3<-FPE
modPH4wc3<-modPH
medPH4wc3<-medPH
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_freqWc3_31122.RData")
FPE5wc3<-FPE
modPH5wc3<-modPH
medPH5wc3<-medPH
# 
nperm<-100
nperm2<-20
npred<-80
#
FPEavgfwc3<-c(rep(0,npred))
for(i in 1:nperm2) FPEavgfwc3<-FPEavgfwc3+FPE1wc3[[i]]$freq
for(i in 1:nperm2) FPEavgfwc3<-FPEavgfwc3+FPE2wc3[[i]]$freq
for(i in 1:nperm2) FPEavgfwc3<-FPEavgfwc3+FPE3wc3[[i]]$freq
for(i in 1:nperm2) FPEavgfwc3<-FPEavgfwc3+FPE4wc3[[i]]$freq
for(i in 1:nperm2) FPEavgfwc3<-FPEavgfwc3+FPE5wc3[[i]]$freq
FPEavgfwc3<-FPEavgfwc3/nperm
FPEsdfwc3<-c(rep(0,npred))
for(i in 1:nperm2) FPEsdfwc3<-FPEsdfwc3+(FPEavgfwc3-FPE1wc3[[i]]$freq)^2
for(i in 1:nperm2) FPEsdfwc3<-FPEsdfwc3+(FPEavgfwc3-FPE2wc3[[i]]$freq)^2
for(i in 1:nperm2) FPEsdfwc3<-FPEsdfwc3+(FPEavgfwc3-FPE3wc3[[i]]$freq)^2
for(i in 1:nperm2) FPEsdfwc3<-FPEsdfwc3+(FPEavgfwc3-FPE4wc3[[i]]$freq)^2
for(i in 1:nperm2) FPEsdfwc3<-FPEsdfwc3+(FPEavgfwc3-FPE5wc3[[i]]$freq)^2
FPEsdfwc3<-sqrt(FPEsdfwc3/(nperm-1))
plot(1:npred,FPE1wc3[[1]]$freq, type="l", xlab='time step', ylab='L1 cumulative predictive error')
lines(1:npred,FPE2wc3[[2]]$freq, type="l", col='red') 
lines(1:npred,FPEavgfwc3, type="l", col='blue')
lines(1:npred,FPEavgfwc3-2*FPEsdfwc3, type="l", lty=2, col='green')
lines(1:npred,FPEavgfwc3+2*FPEsdfwc3, type="l", lty=2, col='green')

FPEavgbc3<-c(rep(0,npred))
for(i in 1:nperm) FPEavgbc3<-FPEavgbc3+FPEc3[[i]]$bayes
FPEavgbc3<-FPEavgbc3/nperm
FPEsdbc3<-c(rep(0,npred))
for(i in 1:nperm) FPEsdbc3<-FPEsdbc3+(FPEavgbc3-FPEc3[[i]]$bayes)^2
FPEsdbc3<-sqrt(FPEsdbc3/(nperm-1))
plot(1:npred,FPEc3[[1]]$bayes, type="l", xlab='time step', ylab='L1 cumulative predictive error')
lines(1:npred,FPEc3[[2]]$bayes, type="l", col='red') 
lines(1:npred,FPEavgbc3, type="l", col='blue')
lines(1:npred,FPEavgbc3-2*FPEsdbc3, type="l", lty=2, col='green')
lines(1:npred,FPEavgbc3+2*FPEsdbc3, type="l", lty=2, col='green')


