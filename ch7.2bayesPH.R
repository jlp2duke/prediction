%		survival data
%	http://www.umass.edu/statdata/statdata/data/uissurv.txt

%	testing psbcGroup prior to implementation
%	cox PH: j permutations of the data order, fit to first 300, then predict ahead in steps of
%	n observations, get predictive errors, then find average and variance of predictive
%	errors over permutations
%	censoring done within permutations
%	for censored data loss function is L1 if pred<cens and 0 if pred>cens
%   three sets of data: original, censor as N(max_obs, var_obs) [draw from censoring for 
%	each obs, if draw<obs censor and no censor otherwise], censor as limit follow-up time

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
set.seed(1100)
nobs<-120
dat<-datorig[sample(1:dim(datorig)[1],nobs),]
nperm<-100
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

set.seed(10010)
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

# bayesian PH by psbcGroup
# NOTE: covariates must be scaled or doesn't work

p<-3
priorPara <- list()
priorPara$eta0 <- 1
priorPara$kappa0 <- 1
priorPara$c0 <- 2
priorPara$r <- 0.5
priorPara$delta <- 0.0001
# no group lasso
priorPara$groupInd <- 1:p
mcmcPara <- list()
mcmcPara$numBeta <- p
mcmcPara$beta.prop.var <- 1
initial <- list()
initial$beta.ini <- rep(0.5, p)
initial$lambdaSq <- 1
initial$sigmaSq <- runif(1, 0.1, 10)
initial$tauSq <- rexp(length(unique(priorPara$groupInd)), rate = initial$lambdaSq/2)
rw <- FALSE
num.reps <- 2000
chain <- 1
thin <- 5
save <- 50

i<-1
for(i in 1:nperm){
	j<-1
	survDat<-list()
	# original data
	for(j in 1:npred){
		datfit<-dat[permord[[i]][1:(nbase+j-1)],]
		survDat$t<-datfit$time
		survDat$di<-datfit$censor
		survDat$x<-cbind(datfit$age,datfit$treat,datfit$site)
		survDat$x<-scale(survDat$x)
		ctr<-attr(survDat$x,"scaled:center")
		scl<-attr(survDat$x,"scaled:scale")
		priorPara$s <- sort(survDat$t[survDat$di == 1])
		priorPara$s <- c(priorPara$s, 2 * max(survDat$t) - max(survDat$t[-which(survDat$t==max(survDat$t))]))
		priorPara$J <- length(priorPara$s)
		initial$h <- rgamma(priorPara$J, 1, 1)
		fitGL <- psbcGL(survDat, priorPara, initial, rw=FALSE, mcmcPara, num.reps, thin, chain, save)
		newobs<-c(dat$age[nbase+j],dat$treat[nbase+j],dat$site[nbase+j])
		newobs<-(newobs-ctr)/scl
		ppred1<- exp(apply(fitGL$beta.p,1,function(x){newobs%*%x}))
		ppred2<- t(apply(fitGL$h.p,1,cumsum))
		# include sum for when j=0
		ppred2<-cbind(0,ppred2[,1:(priorPara$J-1)])
		ppred<-exp(-ppred1*ppred2)*(1-exp(-ppred1*fitGL$h.p))
		ppredavg<-apply(ppred,2,mean)
		if(sum(cumsum(ppredavg)<0.5)==0) ind<-1 else ind<-max(which(cumsum(ppredavg)<0.5))+1
		if(ind==1) my.med<-priorPara$s[ind]/2 else my.med<-(priorPara$s[ind]+priorPara$s[ind-1])/2
		newc<-dat$censor[permord[[i]][nbase+j]]
		newt<-dat$time[permord[[i]][nbase+j]]
		if(newc>0){ 
			perr[[i]]$bayes2[j]<-abs(newt-my.med) 
		} else if (newt >= my.med){
			perr[[i]]$bayes2[j]<-abs(newt-my.med) 
		} else {
			perr[[i]]$bayes2[j]<-0
			}
		if(j<2) FPE[[i]]$bayes2[j]<-perr[[i]]$bayes2[j]
		if(j>1) FPE[[i]]$bayes2[j]<-FPE[[i]]$bayes2[j-1]+perr[[i]]$bayes2[j]
		modPH[[i]]$bayes2[[j]]<-fitGL
		medPH[[i]]$bayes2[j]<-my.med
	}
}
## integrating results across random permutations; load each workspace, save FPE results
# 23 runs
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/uissurvBayPH2_n23.RData")
FPE1b<-FPE
modPH1b<-modPH
medPH1b<-medPH
# 20 runs each
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_survBayes_21012.RData")
FPE2b<-FPE
modPH2b<-modPH
medPH2b<-medPH
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_survBayes_21112.RData")
FPE3b<-FPE
modPH3b<-modPH
medPH3b<-medPH
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_survBayes_21122.RData")
FPE4b<-FPE
modPH4b<-modPH
medPH4b<-medPH
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_survBayes_31122.RData")
FPE5b<-FPE
modPH5b<-modPH
medPH5b<-medPH
# 
nperm2<-20
nperm<-100
#
FPEavgb<-c(rep(0,npred))
for(i in 1:nperm2) FPEavgb<-FPEavgb+FPE1b[[i]]$bayes2
for(i in 1:nperm2) FPEavgb<-FPEavgb+FPE2b[[i]]$bayes2
for(i in 1:nperm2) FPEavgb<-FPEavgb+FPE3b[[i]]$bayes2
for(i in 1:nperm2) FPEavgb<-FPEavgb+FPE4b[[i]]$bayes2
for(i in 1:nperm2) FPEavgb<-FPEavgb+FPE5b[[i]]$bayes2
FPEavgb<-FPEavgb/nperm
FPEsdb<-c(rep(0,npred))
for(i in 1:nperm2) FPEsdb<-FPEsdb+(FPEavgb-FPE1b[[i]]$bayes2)^2
for(i in 1:nperm2) FPEsdb<-FPEsdb+(FPEavgb-FPE2b[[i]]$bayes2)^2
for(i in 1:nperm2) FPEsdb<-FPEsdb+(FPEavgb-FPE3b[[i]]$bayes2)^2
for(i in 1:nperm2) FPEsdb<-FPEsdb+(FPEavgb-FPE4b[[i]]$bayes2)^2
for(i in 1:nperm2) FPEsdb<-FPEsdb+(FPEavgb-FPE5b[[i]]$bayes2)^2
FPEsdb<-sqrt(FPEsdb/(nperm-1))
plot(1:npred,FPE1b[[1]]$bayes2, type="l", xlab='time step', ylab='L1 cumulative predictive error')
lines(1:npred,FPE2b[[2]]$bayes2, type="l", col='red') 
lines(1:npred,FPEavgb, type="l", col='blue')
lines(1:npred,FPEavgb-2*FPEsdb, type="l", lty=2, col='green')
lines(1:npred,FPEavgb+2*FPEsdb, type="l", lty=2, col='green')
# bayCPE_n80_PH.pdf
## redo plot for Esther, bw and better labeling
plot(1:npred,FPE1b[[1]]$bayes2, type="l", xlab='time step', ylab='L1 cumulative predictive error', col='grey75', cex.lab=1.2, cex.axis=1.2, cex=1.2, ylim=c(0,17000))
lines(1:npred,FPE2b[[2]]$bayes2, type="l", col='grey75') 
lines(1:npred,FPEavgb, type="l", col='black')
lines(1:npred,FPEavgb-2*FPEsdb, type="l", lty=2)
lines(1:npred,FPEavgb+2*FPEsdb, type="l", lty=2)
# save as bayCPE_n80_PH_v2.pdf

# remove lots of stuff, then save as ch7computing_survBayes.RData

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

plot(1:npred,FPEavgb, type="l", xlab='time step', ylab='avg L1 cumulative predictive error')
lines(1:npred,FPEavgbc1, col="blue")
lines(1:npred,FPEavgbc3, col="green")
legend('topleft',legend=c('original','rcens0.75','tcens0.75'),col=c('black','blue','green'),lty=1)

## test plots

set.seed(1020)
ind<-sample(1:nperm2,10)
ind2<-sample(1:npred,10)
fitGL<-modPH2[[ind[1]]]$bayes2[[ind2[1]]]
ppred1<- exp(apply(fitGL$beta.p,1,function(x){newobs%*%x}))
ppred2<- t(apply(fitGL$h.p,1,cumsum))
# include sum for when j=0
ppred2<-cbind(0,ppred2[,1:(fitGL$mcmcOutcome$priorPara$J-1)])
ppred<-exp(-ppred1*ppred2)*(1-exp(-ppred1*fitGL$h.p))
ppredavg<-apply(ppred,2,mean)
plot(fitGL$mcmcOutcome$priorPara$s, 1-cumsum(ppredavg), type='l', col=rainbow(10)[1])
lines(c(medPH2[[ind[1]]]$bayes2[ind2[1]],medPH2[[ind[1]]]$bayes2[ind2[1]]),c(0,1),col=rainbow(10)[1])
for(kk in 2:10){
	fitGL<-modPH2[[ind[kk]]]$bayes2[[ind2[kk]]]
	ppred1<- exp(apply(fitGL$beta.p,1,function(x){newobs%*%x}))
	ppred2<- t(apply(fitGL$h.p,1,cumsum))
	# include sum for when j=0
	ppred2<-cbind(0,ppred2[,1:(fitGL$mcmcOutcome$priorPara$J-1)])
	ppred<-exp(-ppred1*ppred2)*(1-exp(-ppred1*fitGL$h.p))
	ppredavg<-apply(ppred,2,mean)
	lines(fitGL$mcmcOutcome$priorPara$s, 1-cumsum(ppredavg), col=rainbow(10)[kk])
	lines(c(medPH2[[ind[kk]]]$bayes2[ind2[kk]],medPH2[[ind[kk]]]$bayes2[ind2[kk]]),c(0,1),col=rainbow(10)[kk])
}
# finalPH_bayes_orig.pdf

## plots for book
## load workspaces for both bayes and freq results, original data
nval<-(nbase+1):nobs
plot(nval,FPEavgf, type='n', ylim=c(0,12000), xlab='sample size', ylab='L1 cumulative predictive error')
polygon(c(rev(nval),nval),c(rev(FPEavgf+FPEsdf),FPEavgf-FPEsdf),col=rgb(89/255,89/255,89/255,0.8),border=NA)
lines(nval,FPEavgf)
lines(nval,FPEavgf+FPEsdf,lty='dashed',col='red')
lines(nval,FPEavgf-FPEsdf,lty='dashed',col='red')
polygon(c(rev(nval),nval),c(rev(FPEavgb+FPEsdb),FPEavgb-FPEsdb),col=rgb(153/255,153/255,153/255,0.8),border=NA)
lines(nval,FPEavgb,col='blue')
lines(nval,FPEavgb+FPEsdb,lty='dashed',col='green')
lines(nval,FPEavgb-FPEsdb,lty='dashed',col='green')

% CPE with +/- SD

plot(nval,FPEavgf, type='n', ylim=c(0,12000), xlab='sample size', ylab='L1 cumulative predictive error')
#polygon(c(rev(nval),nval),c(rev(FPEavgf+FPEsdf),FPEavgf-FPEsdf),col=rgb(89/255,89/255,89/255,0.8),border=NA)
lines(nval,FPEavgf)
lines(nval,FPEavgf+FPEsdf,lty='dashed')
lines(nval,FPEavgf-FPEsdf,lty='dashed')
polygon(c(rev(nval),nval),c(rev(FPEavgb+FPEsdb),FPEavgb-FPEsdb),col=rgb(153/255,153/255,153/255,0.8),border=NA)
lines(nval,FPEavgb,lty='dotted')
#lines(nval,FPEavgb+FPEsdb,lty='dashed',col='green')
#lines(nval,FPEavgb-FPEsdb,lty='dashed',col='green')
legend('topleft',legend=c('freq CPE','freq +/- SD','bayes CPE', 'bayes +/- SD'),col=c('black','black','black',rgb(153/255,153/255,153/255,0.8)),lty=c(1,2,3,1),lwd=c(1,1,1,5))

plot(nval,FPEavgfc1, type='n', ylim=c(0,16000), xlab='sample size', ylab='L1 cumulative predictive error')
lines(nval,FPEavgfc1)
lines(nval,FPEavgfc1+FPEsdfc1,lty='dashed')
lines(nval,FPEavgfc1-FPEsdfc1,lty='dashed')
polygon(c(rev(nval),nval),c(rev(FPEavgbc1+FPEsdbc1),FPEavgbc1-FPEsdbc1),col=rgb(153/255,153/255,153/255,0.8),border=NA)
lines(nval,FPEavgbc1,lty='dotted')
legend('topleft',legend=c('freq CPE','freq +/- SD','bayes CPE', 'bayes +/- SD'),col=c('black','black','black',rgb(153/255,153/255,153/255,0.8)),lty=c(1,2,3,1),lwd=c(1,1,1,5))

plot(nval,FPEavgfc3, type='n', ylim=c(0,16000), xlab='sample size', ylab='L1 cumulative predictive error')
lines(nval,FPEavgfc3)
lines(nval,FPEavgfc3+FPEsdfc3,lty='dashed')
lines(nval,FPEavgfc3-FPEsdfc3,lty='dashed')
polygon(c(rev(nval),nval),c(rev(FPEavgbc3+FPEsdbc3),FPEavgbc3-FPEsdbc3),col=rgb(153/255,153/255,153/255,0.8),border=NA)
lines(nval,FPEavgbc3,lty='dotted')
legend('topleft',legend=c('freq CPE','freq +/- SD','bayes CPE', 'bayes +/- SD'),col=c('black','black','black',rgb(153/255,153/255,153/255,0.8)),lty=c(1,2,3,1),lwd=c(1,1,1,5))

plot(nval,FPEavgf, type='n', ylim=c(0,16000), xlab='sample size', ylab='L1 cumulative predictive error')
lines(nval,FPEavgf)
lines(nval,FPEavgb,lty=2)
lines(nval,FPEavgfc1,lty=3)
lines(nval,FPEavgbc1,lty=4)
lines(nval,FPEavgfc3,lty=5,col=rgb(153/255,153/255,153/255,0.8))
lines(nval,FPEavgbc3,lty=6,col=rgb(153/255,153/255,153/255,0.8))
legend('topleft',legend=c('freq orig','bayes orig','freq c1','bayes c1','freq c3','bayes c3'),lty=c(1:6),col=c(rep('black',4),rep(rgb(153/255,153/255,153/255,0.8),2)))

% KM with 95% CI
plot(modKM[[1]]$freq[[npred]],xlab='time', ylab='probability of survival')
polygon(c(rev(modKM[[1]]$bayes[[npred]][,1]),modKM[[1]]$bayes[[npred]][,1]),c(rev(modKM[[1]]$bayes[[npred]][,2]+1.96*modKM[[1]]$bayes[[npred]][,3]),modKM[[1]]$bayes[[npred]][,2]-1.96*modKM[[1]]$bayes[[npred]][,3]),col=rgb(229/255,229/255,229/255,0.8),border=NA)
lines(modKM[[1]]$bayes[[npred]],lty='dotted')
legend('topright',legend=c('freq KM','freq CI','bayes KM', 'bayes CI'),col=c('black','black','black',rgb(229/255,229/255,229/255,0.8)),lty=c(1,2,3,1))

plot(modKMc1[[1]]$freq[[npred]],xlab='time', ylab='probability of survival')
btime<-modKMc1[[1]]$bayes[[npred]][,1]
bavg<-modKMc1[[1]]$bayes[[npred]][,2]
bsd<-modKMc1[[1]]$bayes[[npred]][,3]
polygon(c(rev(btime),btime),c(rev(bavg+1.96*bsd),bavg-1.96*bsd),col=rgb(229/255,229/255,229/255,0.8),border=NA)
lines(modKMc1[[1]]$bayes[[npred]],lty='dotted')
legend('topright',legend=c('freq KM','freq CI','bayes KM', 'bayes CI'),col=c('black','black','black',rgb(229/255,229/255,229/255,0.8)),lty=c(1,2,3,1))

plot(modKMc2[[1]]$freq[[npred]],xlab='time', ylab='probability of survival')
btime<-modKMc2[[1]]$bayes[[npred]][,1]
bavg<-modKMc2[[1]]$bayes[[npred]][,2]
bsd<-modKMc2[[1]]$bayes[[npred]][,3]
polygon(c(rev(btime),btime),c(rev(bavg+1.96*bsd),bavg-1.96*bsd),col=rgb(229/255,229/255,229/255,0.8),border=NA)
lines(modKMc2[[1]]$bayes[[npred]],lty='dotted')
legend('topright',legend=c('freq KM','freq CI','bayes KM', 'bayes CI'),col=c('black','black','black',rgb(229/255,229/255,229/255,0.8)),lty=c(1,2,3,1))

plot(modKMc3[[1]]$freq[[npred]],xlab='time', ylab='probability of survival')
btime<-modKMc3[[1]]$bayes[[npred]][,1]
bavg<-modKMc3[[1]]$bayes[[npred]][,2]
bsd<-modKMc3[[1]]$bayes[[npred]][,3]
polygon(c(rev(btime),btime),c(rev(bavg+1.96*bsd),bavg-1.96*bsd),col=rgb(229/255,229/255,229/255,0.8),border=NA)
lines(modKMc3[[1]]$bayes[[npred]],lty='dotted')
legend('topright',legend=c('freq KM','freq CI','bayes KM', 'bayes CI'),col=c('black','black','black',rgb(229/255,229/255,229/255,0.8)),lty=c(1,2,3,1))

plot(modKMc4[[1]]$freq[[npred]],xlab='time', ylab='probability of survival')
btime<-modKMc4[[1]]$bayes[[npred]][,1]
bavg<-modKMc4[[1]]$bayes[[npred]][,2]
bsd<-modKMc4[[1]]$bayes[[npred]][,3]
polygon(c(rev(btime),btime),c(rev(bavg+1.96*bsd),bavg-1.96*bsd),col=rgb(229/255,229/255,229/255,0.8),border=NA)
lines(modKMc4[[1]]$bayes[[npred]],lty='dotted')
legend('topright',legend=c('freq KM','freq CI','bayes KM', 'bayes CI'),col=c('black','black','black',rgb(229/255,229/255,229/255,0.8)),lty=c(1,2,3,1))
