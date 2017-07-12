#		survival data
#	http://www.umass.edu/statdata/statdata/data/uissurv.txt

#	testing psbcGroup prior to implementation
#	cox PH: j permutations of the data order, fit to first 300, then predict ahead in steps of
#	n observations, get predictive errors, then find average and variance of predictive
#	errors over permutations
#	censoring done within permutations
#	for censored data loss function is L1 if pred<cens and 0 if pred>cens
#   three sets of data: original, censor as N(max_obs, var_obs) [draw from censoring for 
#	each obs, if draw<obs censor and no censor otherwise], censor as limit follow-up time

seed<-10012
set.seed(seed)

library(LearnBayes)
library(survival)
datorig<-read.csv('uissurv.csv', header=T, row.names=NULL)
datorig$age<-as.numeric(datorig$age)
datorig$agecat<-as.numeric(datorig$age>median(datorig$age))

# implement permute, model, and predict function; use subset of data
nobs<-120
dat<-datorig[sample(1:dim(datorig)[1],nobs),]
nperm<-5
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

# bayesian weibuill regression by LearnBayes
# test first
library(LearnBayes)
# replace version of this function in pkg with the below
weibullregpost<-function(theta,data)
{
s=dim(data); k=s[2]; p=k-2
sp=dim(theta); N=sp[1]
t=data[,1]; c=data[,2]; X=data[,3:k]
sigma=exp(theta[1])
mu=theta[2]
beta=array(theta[3:k],c(N,p))
n=length(t)
o=0*mu
for (i in 1:n)
{
lp=0
for (j in 1:p) lp=lp+beta[j]*X[i,j]
zi=(log(t[i])-mu-lp)/sigma
fi=1/sigma*exp(zi-exp(zi))
Si=exp(-exp(zi))
o=o+c[i]*log(fi)+(1-c[i])*log(Si)
}
return(o)
}

# NOTE: covariates must be scaled or doesn't work

num.reps <- 5000
#chain <- 1
#thin <- 5
#save <- 50

i<-1
for(i in 1:nperm){
	j<-1
	survDat<-list()
	# original data
	for(j in 1:npred){
		datfit<-dat[permord[[i]][1:(nbase+j-1)],]
		survDat$t<-datfit$time
		survDat$di<-newcens[permord[[i]][1:(nbase+j-1)]]
		survDat$x<-cbind(datfit$age,datfit$treat,datfit$site)
		survDat$x<-scale(survDat$x)
		ctr<-attr(survDat$x,"scaled:center")
		scl<-attr(survDat$x,"scaled:scale")
		stmp<-survreg(Surv(t,di)~x,data=survDat,dist="weibull")
		start<-c(0.5,stmp$coefficients)
		d<-cbind(survDat$t,survDat$di,survDat$x)
		fit<-laplace(weibullregpost,mode=start,d)
		proposal<-list(var=fit$var,scale=1.5)
 		bayesfit<-rwmetrop(weibullregpost,proposal,fit$mode,num.reps,d)
		newobs<-c(dat$age[nbase+j],dat$treat[nbase+j],dat$site[nbase+j])
		newobs<-(newobs-ctr)/scl
		sigma=exp(bayesfit$par[,1])
 		mu=bayesfit$par[,2]
 		beta1=bayesfit$par[,3]
 		beta2=bayesfit$par[,4]
		beta3<-bayesfit$par[,5]
		xt<-array(log(datfit$time),dim=c(length(datfit$time),1))
		ppred1<-apply(cbind(beta1,beta2,beta3),1,function(x){newobs%*%x})
		ppred<-exp(-exp(apply(xt,1,function(x){x-mu-ppred1})*(1/sigma)))
		ppredavg<-apply(ppred,2,mean)
		if(sum(cumsum(ppredavg)<0.5)==0) ind<-1 else ind<-max(which(cumsum(ppredavg)<0.5))+1
		my.med<-datfit$time[ind][1]
		newc<-newcens[permord[[i]][nbase+j]]
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
		modW[[i]]$bayes2[[j]]<-bayesfit
		medW[[i]]$bayes2[j]<-my.med
	}
}

# save workspace
#filen<-paste("ch7computing_survBayes_wc",seed,".RData",sep="")
#save.image(file=filen)

## integrating results across permutations; load each workspace, save FPE results
# 23 runs
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_survBayes_wc10010.RData")
FPE1bwc1<-FPE
modPH1bwc1<-modPH
medPH1bwc1<-medPH
# 20 runs each
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_survBayes_wc21012.RData")
FPE2bwc1<-FPE
modPH2bwc1<-modPH
medPH2bwc1<-medPH
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_survBayes_wc21112.RData")
FPE3bwc1<-FPE
modPH3bwc1<-modPH
medPH3bwc1<-medPH
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_survBayes_wc21122.RData")
FPE4bwc1<-FPE
modPH4bwc1<-modPH
medPH4bwc1<-medPH
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_survBayes_wc31122.RData")
FPE5bwc1<-FPE
modPH5bwc1<-modPH
medPH5bwc1<-medPH
# 
nperm2<-20
nperm<-100
#
FPEavgbwc1<-c(rep(0,npred))
for(i in 1:nperm2) FPEavgbwc1<-FPEavgbwc1+FPE1bwc1[[i]]$bayes2
for(i in 1:nperm2) FPEavgbwc1<-FPEavgbwc1+FPE2bwc1[[i]]$bayes2
for(i in 1:nperm2) FPEavgbwc1<-FPEavgbwc1+FPE3bwc1[[i]]$bayes2
for(i in 1:nperm2) FPEavgbwc1<-FPEavgbwc1+FPE4bwc1[[i]]$bayes2
for(i in 1:nperm2) FPEavgbwc1<-FPEavgbwc1+FPE5bwc1[[i]]$bayes2
FPEavgbwc1<-FPEavgbwc1/nperm
FPEsdbwc1<-c(rep(0,npred))
for(i in 1:nperm2) FPEsdbwc1<-FPEsdbwc1+(FPEavgbwc1-FPE1bwc1[[i]]$bayes2)^2
for(i in 1:nperm2) FPEsdbwc1<-FPEsdbwc1+(FPEavgbwc1-FPE2bwc1[[i]]$bayes2)^2
for(i in 1:nperm2) FPEsdbwc1<-FPEsdbwc1+(FPEavgbwc1-FPE3bwc1[[i]]$bayes2)^2
for(i in 1:nperm2) FPEsdbwc1<-FPEsdbwc1+(FPEavgbwc1-FPE4bwc1[[i]]$bayes2)^2
for(i in 1:nperm2) FPEsdbwc1<-FPEsdbwc1+(FPEavgbwc1-FPE5bwc1[[i]]$bayes2)^2
FPEsdbwc1<-sqrt(FPEsdbwc1/(nperm-1))
plot(1:npred,FPE1bwc1[[1]]$bayes2, type="l", xlab='time step', ylab='L1 cumulative predictive error')
lines(1:npred,FPE2bwc1[[2]]$bayes2, type="l", col='red') 
lines(1:npred,FPEavgbwc1, type="l", col='blue')
lines(1:npred,FPEavgbwc1-2*FPEsdbwc1, type="l", lty=2, col='green')
lines(1:npred,FPEavgbwc1+2*FPEsdbwc1, type="l", lty=2, col='green')
# bayCPE_c1_Weib.pdf


