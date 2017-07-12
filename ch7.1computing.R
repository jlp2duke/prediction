%	ch7.1	survival data
%	http://www.umass.edu/statdata/statdata/data/uissurv.txt

%	KM: j permutations of the data order, fit to first 300, then predict ahead in steps of
%	n observations, get predictive errors, then find average and variance of predictive
%	errors over permutations
%	censoring done within permutations
%	for censored data loss function is L1 if pred<cens and 0 if pred>cens
%   three sets of data: original, censor as N(max_obs, var_obs) [draw from censoring for 
%	each obs, if draw<obs censor and no censor otherwise], censor as limit follow-up time

library(survival)
library(KMsurv)
library(nlme)
library(sm)
library(km.ci)

datorig<-read.csv('uissurv.csv', header=T, row.names=NULL)
sm.density.compare(datorig$time,datorig$censor, xlab="survival time")
attach(datorig)
my.surv<-Surv(time, censor)
my.fit<-survfit(my.surv ~ 1)
plot(my.fit)
my.fit
#Call: survfit(formula = my.surv ~ 1)
#records   n.max n.start  events  median 0.95LCL 0.95UCL 
#    628     628     628     508     166     148     184 
Tmp<-tail(capture.output(my.fit),2)
Res<-read.table(z <- textConnection(Tmp), header=T)
close(z)
my.med<-Res$median

# implement permute, model, and predict function; use subset of data
set.seed(1100)
nobs<-80
dat<-datorig[sample(1:dim(datorig)[1],nobs),]
nperm<-100
nbase<-20
npred<-nobs-nbase

# variables to hold final results (all runs)
permord<-vector('list', length=nperm)
modKM<-permord
medKM<-permord
perr<-permord
FPE<-permord
for(i in 1:nperm) permord[[i]]<-c(rep(0,nobs))
for(i in 1:nperm) modKM[[i]]<-list("freq"=vector('list', length=npred),"bayes"=vector('list', length=npred))
for(i in 1:nperm) medKM[[i]]<-list("freq"=c(rep(99,npred)),"bayes"=c(rep(99,npred)))
for(i in 1:nperm) perr[[i]]<-list("freq"=c(rep(99,npred)),"bayes"=c(rep(99,npred)))
for(i in 1:nperm) FPE[[i]]<-list("freq"=c(rep(99,npred)),"bayes"=c(rep(99,npred)))
modKMc1<-modKMc2<-modKMc3<-modKMc4<-modKM
medKMc1<-medKMc2<-medKMc3<-medKMc4<-medKM
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
# frequentist KM
i<-1
for(i in 1:nperm){
	j<-1
	# original data
	for(j in 1:npred){
		time<-dat$time[permord[[i]]][1:(nbase+j-1)]
		cens<-dat$cens[permord[[i]]][1:(nbase+j-1)]
		my.surv<-Surv(time,cens)
		my.fit<-survfit(my.surv ~ 1)
		Tmp<-tail(capture.output(my.fit),2)
		Res<-read.table(z <- textConnection(Tmp), header=T)
		close(z)
		my.med<-Res$median
		ptime<-dat$time[permord[[i]]][nbase+j]
		pcens<-dat$cens[permord[[i]]][nbase+j]
		if(pcens>0){ 
			perr[[i]]$freq[j]<-abs(ptime-my.med) 
		} else if (ptime <= my.med){
			perr[[i]]$freq[j]<-abs(ptime-my.med) 
		} else {
			perr[[i]]$freq[j]<-0
			}
		if(j<2) FPE[[i]]$freq[j]<-perr[[i]]$freq[j]
		if(j>1) FPE[[i]]$freq[j]<-FPE[[i]]$freq[j-1]+perr[[i]]$freq[j]
		modKM[[i]]$freq[[j]]<-my.fit
		medKM[[i]]$freq[j]<-my.med
	}
	j<-1
	# random censor 1
	for(j in 1:npred){
		time<-dat$time[permord[[i]]][1:(nbase+j-1)]
		cens<-newcens[permord[[i]]][1:(nbase+j-1)]
		my.surv<-Surv(time,cens)
		my.fit<-survfit(my.surv ~ 1)
		Tmp<-tail(capture.output(my.fit),2)
		Res<-read.table(z <- textConnection(Tmp), header=T)
		close(z)
		my.med<-Res$median
		ptime<-dat$time[permord[[i]]][nbase+j]
		pcens<-newcens[permord[[i]]][nbase+j]
		if(pcens>0){ 
			perrc1[[i]]$freq[j]<-abs(ptime-my.med) 
		} else if (ptime <= my.med){
			perrc1[[i]]$freq[j]<-abs(ptime-my.med) 
		} else {
			perrc1[[i]]$freq[j]<-0
			}
		if(j<2) FPEc1[[i]]$freq[j]<-perrc1[[i]]$freq[j]
		if(j>1) FPEc1[[i]]$freq[j]<-FPEc1[[i]]$freq[j-1]+perrc1[[i]]$freq[j]
		modKMc1[[i]]$freq[[j]]<-my.fit
		medKMc1[[i]]$freq[j]<-my.med
	}	
	j<-1
	# random censor 2
	for(j in 1:npred){
		time<-dat$time[permord[[i]]][1:(nbase+j-1)]
		cens<-newcens2[permord[[i]]][1:(nbase+j-1)]
		my.surv<-Surv(time,cens)
		my.fit<-survfit(my.surv ~ 1)
		Tmp<-tail(capture.output(my.fit),2)
		Res<-read.table(z <- textConnection(Tmp), header=T)
		close(z)
		my.med<-Res$median
		ptime<-dat$time[permord[[i]]][nbase+j]
		pcens<-newcens2[permord[[i]]][nbase+j]
		if(pcens>0){ 
			perrc2[[i]]$freq[j]<-abs(ptime-my.med) 
		} else if (ptime <= my.med){
			perrc2[[i]]$freq[j]<-abs(ptime-my.med) 
		} else {
			perrc2[[i]]$freq[j]<-0
			}
		if(j<2) FPEc2[[i]]$freq[j]<-perrc2[[i]]$freq[j]
		if(j>1) FPEc2[[i]]$freq[j]<-FPEc2[[i]]$freq[j-1]+perrc2[[i]]$freq[j]
		modKMc2[[i]]$freq[[j]]<-my.fit
		medKMc2[[i]]$freq[j]<-my.med
	}	
	j<-1
	# truncation censor 1
	for(j in 1:npred){
		time<-dat$time[permord[[i]]][1:(nbase+j-1)]
		cens<-lcens[permord[[i]]][1:(nbase+j-1)]
		my.surv<-Surv(time,cens)
		my.fit<-survfit(my.surv ~ 1)
		Tmp<-tail(capture.output(my.fit),2)
		Res<-read.table(z <- textConnection(Tmp), header=T)
		close(z)
		my.med<-Res$median
		ptime<-dat$time[permord[[i]]][nbase+j]
		pcens<-lcens[permord[[i]]][nbase+j]
		if(pcens>0){ 
			perrc3[[i]]$freq[j]<-abs(ptime-my.med) 
		} else if (ptime <= my.med){
			perrc3[[i]]$freq[j]<-abs(ptime-my.med) 
		} else {
			perrc3[[i]]$freq[j]<-0
			}
		if(j<2) FPEc3[[i]]$freq[j]<-perrc3[[i]]$freq[j]
		if(j>1) FPEc3[[i]]$freq[j]<-FPEc3[[i]]$freq[j-1]+perrc3[[i]]$freq[j]
		modKMc3[[i]]$freq[[j]]<-my.fit
		medKMc3[[i]]$freq[j]<-my.med
	}	
	j<-1
	# truncation censor 2
	for(j in 1:npred){
		time<-dat$time[permord[[i]]][1:(nbase+j-1)]
		cens<-lcens2[permord[[i]]][1:(nbase+j-1)]
		my.surv<-Surv(time,cens)
		my.fit<-survfit(my.surv ~ 1)
		Tmp<-tail(capture.output(my.fit),2)
		Res<-read.table(z <- textConnection(Tmp), header=T)
		close(z)
		my.med<-Res$median
		ptime<-dat$time[permord[[i]]][nbase+j]
		pcens<-lcens2[permord[[i]]][nbase+j]
		if(pcens>0){ 
			perrc4[[i]]$freq[j]<-abs(ptime-my.med) 
		} else if (ptime <= my.med){
			perrc4[[i]]$freq[j]<-abs(ptime-my.med) 
		} else {
			perrc4[[i]]$freq[j]<-0
			}
		if(j<2) FPEc4[[i]]$freq[j]<-perrc4[[i]]$freq[j]
		if(j>1) FPEc4[[i]]$freq[j]<-FPEc4[[i]]$freq[j-1]+perrc4[[i]]$freq[j]
		modKMc4[[i]]$freq[[j]]<-my.fit
		medKMc4[[i]]$freq[j]<-my.med
	}	
}
#
FPEavgf<-c(rep(0,npred))
for(i in 1:nperm) FPEavgf<-FPEavgf+FPE[[i]]$freq
FPEavgf<-FPEavgf/nperm
FPEsdf<-c(rep(0,npred))
for(i in 1:nperm) FPEsdf<-FPEsdf+(FPEavgf-FPE[[i]]$freq)^2
FPEsdf<-sqrt(FPEsdf/(nperm-1))
plot(1:npred,FPE[[1]]$freq, type="l", xlab='time step', ylab='L1 cumulative predictive error')
lines(1:npred,FPE[[2]]$freq, type="l", col='red') 
lines(1:npred,FPEavgf, type="l", col='blue')
lines(1:npred,FPEavgf-2*FPEsdf, type="l", lty=2, col='green')
lines(1:npred,FPEavgf+2*FPEsdf, type="l", lty=2, col='green')

FPEavgfc1<-c(rep(0,npred))
for(i in 1:nperm) FPEavgfc1<-FPEavgfc1+FPEc1[[i]]$freq
FPEavgfc1<-FPEavgfc1/nperm
FPEsdfc1<-c(rep(0,npred))
for(i in 1:nperm) FPEsdfc1<-FPEsdfc1+(FPEavgfc1-FPEc1[[i]]$freq)^2
FPEsdfc1<-sqrt(FPEsdfc1/(nperm-1))
plot(1:npred,FPEc1[[1]]$freq, type="l", xlab='time step', ylab='L1 cumulative predictive error')
lines(1:npred,FPEc1[[2]]$freq, type="l", col='red') 
lines(1:npred,FPEavgfc1, type="l", col='blue')
lines(1:npred,FPEavgfc1-2*FPEsdfc1, type="l", lty=2, col='green')
lines(1:npred,FPEavgfc1+2*FPEsdfc1, type="l", lty=2, col='green')

FPEavgfc2<-c(rep(0,npred))
for(i in 1:nperm) FPEavgfc2<-FPEavgfc2+FPEc2[[i]]$freq
FPEavgfc2<-FPEavgfc2/nperm
FPEsdfc2<-c(rep(0,npred))
for(i in 1:nperm) FPEsdfc2<-FPEsdfc2+(FPEavgfc2-FPEc2[[i]]$freq)^2
FPEsdfc2<-sqrt(FPEsdfc2/(nperm-1))
plot(1:npred,FPEc2[[1]]$freq, type="l", xlab='time step', ylab='L1 cumulative predictive error')
lines(1:npred,FPEc2[[2]]$freq, type="l", col='red') 
lines(1:npred,FPEavgfc2, type="l", col='blue')
lines(1:npred,FPEavgfc2-2*FPEsdfc2, type="l", lty=2, col='green')
lines(1:npred,FPEavgfc2+2*FPEsdfc2, type="l", lty=2, col='green')

FPEavgfc3<-c(rep(0,npred))
for(i in 1:nperm) FPEavgfc3<-FPEavgfc3+FPEc3[[i]]$freq
FPEavgfc3<-FPEavgfc3/nperm
FPEsdfc3<-c(rep(0,npred))
for(i in 1:nperm) FPEsdfc3<-FPEsdfc3+(FPEavgfc3-FPEc3[[i]]$freq)^2
FPEsdfc3<-sqrt(FPEsdfc3/(nperm-1))
plot(1:npred,FPEc3[[1]]$freq, type="l", xlab='time step', ylab='L1 cumulative predictive error')
lines(1:npred,FPEc3[[2]]$freq, type="l", col='red') 
lines(1:npred,FPEavgfc3, type="l", col='blue')
lines(1:npred,FPEavgfc3-2*FPEsdfc3, type="l", lty=2, col='green')
lines(1:npred,FPEavgfc3+2*FPEsdfc3, type="l", lty=2, col='green')

FPEavgfc4<-c(rep(0,npred))
for(i in 1:nperm) FPEavgfc4<-FPEavgfc4+FPEc4[[i]]$freq
FPEavgfc4<-FPEavgfc4/nperm
FPEsdfc4<-c(rep(0,npred))
for(i in 1:nperm) FPEsdfc4<-FPEsdfc4+(FPEavgfc4-FPEc4[[i]]$freq)^2
FPEsdfc4<-sqrt(FPEsdfc4/(nperm-1))
plot(1:npred,FPEc4[[1]]$freq, type="l", xlab='time step', ylab='L1 cumulative predictive error')
lines(1:npred,FPEc4[[2]]$freq, type="l", col='red') 
lines(1:npred,FPEavgfc4, type="l", col='blue')
lines(1:npred,FPEavgfc4-2*FPEsdfc4, type="l", lty=2, col='green')
lines(1:npred,FPEavgfc4+2*FPEsdfc4, type="l", lty=2, col='green')

plot(1:npred,FPEavgf, type="l", xlab='time step', ylab='avg L1 cumulative predictive error')
lines(1:npred,FPEavgfc1, col="blue")
lines(1:npred,FPEavgfc2, col="red")
lines(1:npred,FPEavgfc3, col="green")
lines(1:npred,FPEavgfc4, col="orange")
legend('topleft',legend=c('original','rcens0.75','rcens0.95','tcens0.75','tcens0.95'),col=c('black','blue','red','green','orange'),lty=1)

## KM with median plots

set.seed(1020)
ind<-sample(1:nperm,10)
plot(modKMc3[[ind[1]]]$freq[[npred]],conf.int=F,mark.time=T,col=rainbow(10)[1])
lines(c(medKMc3[[ind[1]]]$freq[npred],medKMc3[[ind[1]]]$freq[npred]),c(0,1),col=rainbow(10)[1])
for(kk in 2:10){
	lines(modKMc3[[ind[kk]]]$freq[[npred]],conf.int=F,mark.time=T,col=rainbow(10)[kk])
	lines(c(medKMc3[[ind[kk]]]$freq[npred],medKMc3[[ind[kk]]]$freq[npred]),c(0,1),col=rainbow(10)[kk])
}


# repeat above for Bayes KM ala Susarla and Van Ryzin (1976)
# code from Tim Hanson for Bayes

S0=function(t){exp(-t.hat*t)}
####################################################
# svr outputs E{S(t)|data} and sd{S(t)|data} using #
# formulae in Susarla and Van Ryzin (1976)         #
# t=argument of survival function, y=event times   #
# d=non-censoring indicator, c0 is DP mass         #
####################################################
svr=function(t,y,d,c0){ 
 o=c(1,1) 
 o[1]=(c0*S0(t)+sum(y>t))/(c0+length(y))
 o[2]=o[1]*(c0*S0(t)+sum(y>t)+1)/(c0+length(y)+1)
 c=y[d==0]; u=c(sort(unique(c)),t+1); j=1
 while(u[j]<=t){
  d1=c0*S0(u[j])+sum(y>=u[j])
  d2=d1/(d1-sum(c==u[j]))
  o[1]=o[1]*d2
  o[2]=o[2]*d2*(d1+1)/(d1+1-sum(c==u[j]))
  j=j+1
 }
 o[2]=sqrt(o[2]-o[1]^2)
 o
}

i<-1
for(i in 1:nperm){
	j<-1
	# original data
	for(j in 1:npred){
		y<-dat$time[permord[[i]]][1:(nbase+j-1)]
		d<-dat$cens[permord[[i]]][1:(nbase+j-1)]
		t.hat=sum(d==1)/sum(y)
		t=seq(0,1200,length=nobs*2); s=matrix(0,nrow=nobs*2,ncol=2)
		# with uniform prior, DP all alphas as 1
		for(ii in 1:(nbase+j-1)){s[ii,]=svr(t[ii],y,d,1)}
		index<-which(abs(s[,1]-0.5)==min(abs(s[,1]-0.5)))	
		my.med<-t[index]
		ptime<-dat$time[permord[[i]]][nbase+j]
		pcens<-dat$cens[permord[[i]]][nbase+j]
		if(pcens>0){ 
			perr[[i]]$bayes[j]<-abs(ptime-my.med) 
		} else if (ptime <= my.med){
			perr[[i]]$bayes[j]<-abs(ptime-my.med) 
		} else {
			perr[[i]]$bayes[j]<-0
			}
		if(j<2) FPE[[i]]$bayes[j]<-perr[[i]]$bayes[j]
		if(j>1) FPE[[i]]$bayes[j]<-FPE[[i]]$bayes[j-1]+perr[[i]]$bayes[j]
		modKM[[i]]$bayes[[j]]<-cbind(t,s)
		medKM[[i]]$bayes[j]<-my.med
	}
	j<-1
	# random censor 1
	for(j in 1:npred){
		y<-dat$time[permord[[i]]][1:(nbase+j-1)]
		d<-newcens[permord[[i]]][1:(nbase+j-1)]
		t.hat=sum(d==1)/sum(y)
		t=seq(0,1200,length=nobs*2); s=matrix(0,nrow=nobs*2,ncol=2)
		# with uniform prior, DP all alphas as 1
		for(ii in 1:(nbase+j-1)){s[ii,]=svr(t[ii],y,d,1)}
		index<-which(abs(s[,1]-0.5)==min(abs(s[,1]-0.5)))	
		my.med<-t[index]
		ptime<-dat$time[permord[[i]]][nbase+j]
		pcens<-newcens[permord[[i]]][nbase+j]
		if(pcens>0){ 
			perrc1[[i]]$bayes[j]<-abs(ptime-my.med) 
		} else if (ptime <= my.med){
			perrc1[[i]]$bayes[j]<-abs(ptime-my.med) 
		} else {
			perrc1[[i]]$bayes[j]<-0
			}
		if(j<2) FPEc1[[i]]$bayes[j]<-perrc1[[i]]$bayes[j]
		if(j>1) FPEc1[[i]]$bayes[j]<-FPEc1[[i]]$bayes[j-1]+perrc1[[i]]$bayes[j]
		modKMc1[[i]]$bayes[[j]]<-cbind(t,s)
		medKMc1[[i]]$bayes[j]<-my.med
	}
	j<-1
	# random censor 2
	for(j in 1:npred){
		y<-dat$time[permord[[i]]][1:(nbase+j-1)]
		d<-newcens2[permord[[i]]][1:(nbase+j-1)]
		t.hat=sum(d==1)/sum(y)
		t=seq(0,1200,length=nobs*2); s=matrix(0,nrow=nobs*2,ncol=2)
		# with uniform prior, DP all alphas as 1
		for(ii in 1:(nbase+j-1)){s[ii,]=svr(t[ii],y,d,1)}
		index<-which(abs(s[,1]-0.5)==min(abs(s[,1]-0.5)))	
		my.med<-t[index]
		ptime<-dat$time[permord[[i]]][nbase+j]
		pcens<-newcens2[permord[[i]]][nbase+j]
		if(pcens>0){ 
			perrc2[[i]]$bayes[j]<-abs(ptime-my.med) 
		} else if (ptime <= my.med){
			perrc2[[i]]$bayes[j]<-abs(ptime-my.med) 
		} else {
			perrc2[[i]]$bayes[j]<-0
			}
		if(j<2) FPEc2[[i]]$bayes[j]<-perrc2[[i]]$bayes[j]
		if(j>1) FPEc2[[i]]$bayes[j]<-FPEc2[[i]]$bayes[j-1]+perrc2[[i]]$bayes[j]
		modKMc2[[i]]$bayes[[j]]<-cbind(t,s)
		medKMc2[[i]]$bayes[j]<-my.med
	}
	j<-1
	# truncation censor 1
	for(j in 1:npred){
		y<-dat$time[permord[[i]]][1:(nbase+j-1)]
		d<-lcens[permord[[i]]][1:(nbase+j-1)]
		t.hat=sum(d==1)/sum(y)
		t=seq(0,1200,length=nobs*2); s=matrix(0,nrow=nobs*2,ncol=2)
		# with uniform prior, DP all alphas as 1
		for(ii in 1:(nbase+j-1)){s[ii,]=svr(t[ii],y,d,1)}
		index<-which(abs(s[,1]-0.5)==min(abs(s[,1]-0.5)))	
		my.med<-t[index]
		ptime<-dat$time[permord[[i]]][nbase+j]
		pcens<-lcens[permord[[i]]][nbase+j]
		if(pcens>0){ 
			perrc3[[i]]$bayes[j]<-abs(ptime-my.med) 
		} else if (ptime <= my.med){
			perrc3[[i]]$bayes[j]<-abs(ptime-my.med) 
		} else {
			perrc3[[i]]$bayes[j]<-0
			}
		if(j<2) FPEc3[[i]]$bayes[j]<-perrc3[[i]]$bayes[j]
		if(j>1) FPEc3[[i]]$bayes[j]<-FPEc3[[i]]$bayes[j-1]+perrc3[[i]]$bayes[j]
		modKMc3[[i]]$bayes[[j]]<-cbind(t,s)
		medKMc3[[i]]$bayes[j]<-my.med
	}
	j<-1
	# truncation censor 2
	for(j in 1:npred){
		y<-dat$time[permord[[i]]][1:(nbase+j-1)]
		d<-lcens2[permord[[i]]][1:(nbase+j-1)]
		t.hat=sum(d==1)/sum(y)
		t=seq(0,1200,length=nobs*2); s=matrix(0,nrow=nobs*2,ncol=2)
		# with uniform prior, DP all alphas as 1
		for(ii in 1:(nbase+j-1)){s[ii,]=svr(t[ii],y,d,1)}
		index<-which(abs(s[,1]-0.5)==min(abs(s[,1]-0.5)))	
		my.med<-t[index]
		ptime<-dat$time[permord[[i]]][nbase+j]
		pcens<-lcens2[permord[[i]]][nbase+j]
		if(pcens>0){ 
			perrc4[[i]]$bayes[j]<-abs(ptime-my.med) 
		} else if (ptime <= my.med){
			perrc4[[i]]$bayes[j]<-abs(ptime-my.med) 
		} else {
			perrc4[[i]]$bayes[j]<-0
			}
		if(j<2) FPEc4[[i]]$bayes[j]<-perrc4[[i]]$bayes[j]
		if(j>1) FPEc4[[i]]$bayes[j]<-FPEc4[[i]]$bayes[j-1]+perrc4[[i]]$bayes[j]
		modKMc4[[i]]$bayes[[j]]<-cbind(t,s)
		medKMc4[[i]]$bayes[j]<-my.med
	}
}
#
FPEavgb<-c(rep(0,npred))
for(i in 1:nperm) FPEavgb<-FPEavgb+FPE[[i]]$bayes
FPEavgb<-FPEavgb/nperm
FPEsdb<-c(rep(0,npred))
for(i in 1:nperm) FPEsdb<-FPEsdb+(FPEavgb-FPE[[i]]$bayes)^2
FPEsdb<-sqrt(FPEsdb/(nperm-1))
plot(1:npred,FPE[[1]]$bayes, type="l", xlab='time step', ylab='L1 cumulative predictive error')
lines(1:npred,FPE[[2]]$bayes, type="l", col='red') 
lines(1:npred,FPEavgb, type="l", col='blue')
lines(1:npred,FPEavgb-2*FPEsdb, type="l", lty=2, col='green')
lines(1:npred,FPEavgb+2*FPEsdb, type="l", lty=2, col='green')

FPEavgbc1<-c(rep(0,npred))
for(i in 1:nperm) FPEavgbc1<-FPEavgbc1+FPEc1[[i]]$bayes
FPEavgbc1<-FPEavgbc1/nperm
FPEsdbc1<-c(rep(0,npred))
for(i in 1:nperm) FPEsdbc1<-FPEsdbc1+(FPEavgbc1-FPEc1[[i]]$bayes)^2
FPEsdbc1<-sqrt(FPEsdbc1/(nperm-1))
plot(1:npred,FPEc1[[1]]$bayes, type="l", xlab='time step', ylab='L1 cumulative predictive error')
lines(1:npred,FPEc1[[2]]$bayes, type="l", col='red') 
lines(1:npred,FPEavgbc1, type="l", col='blue')
lines(1:npred,FPEavgbc1-2*FPEsdbc1, type="l", lty=2, col='green')
lines(1:npred,FPEavgbc1+2*FPEsdbc1, type="l", lty=2, col='green')

FPEavgbc2<-c(rep(0,npred))
for(i in 1:nperm) FPEavgbc2<-FPEavgbc2+FPEc2[[i]]$bayes
FPEavgbc2<-FPEavgbc2/nperm
FPEsdbc2<-c(rep(0,npred))
for(i in 1:nperm) FPEsdbc2<-FPEsdbc2+(FPEavgbc2-FPEc2[[i]]$bayes)^2
FPEsdbc2<-sqrt(FPEsdbc2/(nperm-1))
plot(1:npred,FPEc2[[1]]$bayes, type="l", xlab='time step', ylab='L1 cumulative predictive error')
lines(1:npred,FPEc2[[2]]$bayes, type="l", col='red') 
lines(1:npred,FPEavgbc2, type="l", col='blue')
lines(1:npred,FPEavgbc2-2*FPEsdbc2, type="l", lty=2, col='green')
lines(1:npred,FPEavgbc2+2*FPEsdbc2, type="l", lty=2, col='green')

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

FPEavgbc4<-c(rep(0,npred))
for(i in 1:nperm) FPEavgbc4<-FPEavgbc4+FPEc4[[i]]$bayes
FPEavgbc4<-FPEavgbc4/nperm
FPEsdbc4<-c(rep(0,npred))
for(i in 1:nperm) FPEsdbc4<-FPEsdbc4+(FPEavgbc4-FPEc4[[i]]$bayes)^2
FPEsdbc4<-sqrt(FPEsdbc4/(nperm-1))
plot(1:npred,FPEc4[[1]]$bayes, type="l", xlab='time step', ylab='L1 cumulative predictive error')
lines(1:npred,FPEc4[[2]]$bayes, type="l", col='red') 
lines(1:npred,FPEavgbc4, type="l", col='blue')
lines(1:npred,FPEavgbc4-2*FPEsdbc4, type="l", lty=2, col='green')
lines(1:npred,FPEavgbc4+2*FPEsdbc4, type="l", lty=2, col='green')

plot(1:npred,FPEavgb, type="l", xlab='time step', ylab='avg L1 cumulative predictive error')
lines(1:npred,FPEavgbc1, col="blue")
lines(1:npred,FPEavgbc2, col="red")
lines(1:npred,FPEavgbc3, col="green")
lines(1:npred,FPEavgbc4, col="orange")
legend('topleft',legend=c('original','rcens0.75','rcens0.95','tcens0.75','tcens0.95'),col=c('black','blue','red','green','orange'),lty=1)

## KM with median plots

set.seed(1020)
ind<-sample(1:nperm,10)
plot(modKMc3[[ind[1]]]$bayes[[npred]],col=rainbow(10)[1],type="l")
lines(c(medKMc3[[ind[1]]]$bayes[npred],medKMc3[[ind[1]]]$bayes[npred]),c(0,1),col=rainbow(10)[1])
for(kk in 2:10){
	lines(modKMc3[[ind[kk]]]$bayes[[npred]],col=rainbow(10)[kk])
	lines(c(medKMc3[[ind[kk]]]$bayes[npred],medKMc3[[ind[kk]]]$bayes[npred]),c(0,1),col=rainbow(10)[kk])
}

## plots for book
NOTE: difference between freq and bay very small

nval<-(nbase+1):nobs
plot(nval,FPEavgf, type='n', ylim=c(35,4500), xlab='sample size', ylab='L1 cumulative predictive error')
polygon(c(rev(nval),nval),c(rev(FPEavgf+FPEsdf),FPEavgf-FPEsdf),col=rgb(89/255,89/255,89/255,0.8),border=NA)
lines(nval,FPEavgf)
lines(nval,FPEavgf+FPEsdf,lty='dashed',col='red')
lines(nval,FPEavgf-FPEsdf,lty='dashed',col='red')
polygon(c(rev(nval),nval),c(rev(FPEavgb+FPEsdb),FPEavgb-FPEsdb),col=rgb(153/255,153/255,153/255,0.8),border=NA)
lines(nval,FPEavgb,col='blue')
lines(nval,FPEavgb+FPEsdb,lty='dashed',col='green')
lines(nval,FPEavgb-FPEsdb,lty='dashed',col='green')

% CPE with +/- SD

plot(nval,FPEavgf, type='n', ylim=c(35,4500), xlab='sample size', ylab='L1 cumulative predictive error')
#polygon(c(rev(nval),nval),c(rev(FPEavgf+FPEsdf),FPEavgf-FPEsdf),col=rgb(89/255,89/255,89/255,0.8),border=NA)
lines(nval,FPEavgf)
lines(nval,FPEavgf+FPEsdf,lty='dashed')
lines(nval,FPEavgf-FPEsdf,lty='dashed')
polygon(c(rev(nval),nval),c(rev(FPEavgb+FPEsdb),FPEavgb-FPEsdb),col=rgb(153/255,153/255,153/255,0.8),border=NA)
lines(nval,FPEavgb,lty='dotted')
#lines(nval,FPEavgb+FPEsdb,lty='dashed',col='green')
#lines(nval,FPEavgb-FPEsdb,lty='dashed',col='green')
legend('topleft',legend=c('freq CPE','freq +/- SD','bayes CPE', 'bayes +/- SD'),col=c('black','black','black',rgb(153/255,153/255,153/255,0.8)),lty=c(1,2,3,1),lwd=c(1,1,1,5))

plot(nval,FPEavgfc1, type='n', ylim=c(35,4500), xlab='sample size', ylab='L1 cumulative predictive error')
lines(nval,FPEavgfc1)
lines(nval,FPEavgfc1+FPEsdfc1,lty='dashed')
lines(nval,FPEavgfc1-FPEsdfc1,lty='dashed')
polygon(c(rev(nval),nval),c(rev(FPEavgbc1+FPEsdbc1),FPEavgbc1-FPEsdbc1),col=rgb(153/255,153/255,153/255,0.8),border=NA)
lines(nval,FPEavgbc1,lty='dotted')
legend('topleft',legend=c('freq CPE','freq +/- SD','bayes CPE', 'bayes +/- SD'),col=c('black','black','black',rgb(153/255,153/255,153/255,0.8)),lty=c(1,2,3,1),lwd=c(1,1,1,5))

plot(nval,FPEavgfc2, type='n', ylim=c(35,4500), xlab='sample size', ylab='L1 cumulative predictive error')
lines(nval,FPEavgfc2)
lines(nval,FPEavgfc2+FPEsdfc2,lty='dashed')
lines(nval,FPEavgfc2-FPEsdfc2,lty='dashed')
polygon(c(rev(nval),nval),c(rev(FPEavgbc2+FPEsdbc2),FPEavgbc2-FPEsdbc2),col=rgb(153/255,153/255,153/255,0.8),border=NA)
lines(nval,FPEavgbc2,lty='dotted')
legend('topleft',legend=c('freq CPE','freq +/- SD','bayes CPE', 'bayes +/- SD'),col=c('black','black','black',rgb(153/255,153/255,153/255,0.8)),lty=c(1,2,3,1),lwd=c(1,1,1,5))

plot(nval,FPEavgfc3, type='n', ylim=c(35,4500), xlab='sample size', ylab='L1 cumulative predictive error')
lines(nval,FPEavgfc3)
lines(nval,FPEavgfc3+FPEsdfc3,lty='dashed')
lines(nval,FPEavgfc3-FPEsdfc3,lty='dashed')
polygon(c(rev(nval),nval),c(rev(FPEavgbc3+FPEsdbc3),FPEavgbc3-FPEsdbc3),col=rgb(153/255,153/255,153/255,0.8),border=NA)
lines(nval,FPEavgbc3,lty='dotted')
legend('topleft',legend=c('freq CPE','freq +/- SD','bayes CPE', 'bayes +/- SD'),col=c('black','black','black',rgb(153/255,153/255,153/255,0.8)),lty=c(1,2,3,1),lwd=c(1,1,1,5))

plot(nval,FPEavgfc4, type='n', ylim=c(35,4500), xlab='sample size', ylab='L1 cumulative predictive error')
lines(nval,FPEavgfc4)
lines(nval,FPEavgfc4+FPEsdfc4,lty='dashed')
lines(nval,FPEavgfc4-FPEsdfc4,lty='dashed')
polygon(c(rev(nval),nval),c(rev(FPEavgbc4+FPEsdbc4),FPEavgbc4-FPEsdbc4),col=rgb(153/255,153/255,153/255,0.8),border=NA)
lines(nval,FPEavgbc4,lty='dotted')
legend('topleft',legend=c('freq CPE','freq +/- SD','bayes CPE', 'bayes +/- SD'),col=c('black','black','black',rgb(153/255,153/255,153/255,0.8)),lty=c(1,2,3,1),lwd=c(1,1,1,5))

plot(nval,FPEavgf, type='n', ylim=c(35,4500), xlab='sample size', ylab='L1 cumulative predictive error')
lines(nval,FPEavgf)
lines(nval,FPEavgb,lty=2)
lines(nval,FPEavgfc1,lty=3)
lines(nval,FPEavgbc1,lty=4)
lines(nval,FPEavgfc2,lty=5)
lines(nval,FPEavgbc2,lty=6)
lines(nval,FPEavgfc3,lty=1,col=rgb(153/255,153/255,153/255,0.8))
lines(nval,FPEavgbc3,lty=2,col=rgb(153/255,153/255,153/255,0.8))
lines(nval,FPEavgfc4,lty=3,col=rgb(153/255,153/255,153/255,0.8))
lines(nval,FPEavgbc4,lty=4,col=rgb(153/255,153/255,153/255,0.8))
legend('topleft',legend=c('freq orig','bayes orig','freq c1','bayes c1','freq c2','bayes c2','freq c3','bayes c3','freq c4','bayes c4'),lty=c(1:6,1:4),col=c(rep('black',6),rep(rgb(153/255,153/255,153/255,0.8),4))

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


######################### original code from Tim Hanson #################################

library(KMsurv)
library(survival)
data(tongue)

y=tongue$time[tongue$type==1]; d=tongue$delta[tongue$type==1]
t.hat=sum(d==1)/sum(y)
S0=function(t){exp(-t.hat*t)}
# centering dist'n is exponential with MLE rate

####################################################
# svr outputs E{S(t)|data} and sd{S(t)|data} using #
# formulae in Susarla and Van Ryzin (1976)         #
# t=argument of survival function, y=event times   #
# d=non-censoring indicator, c0 is DP mass         #
####################################################
svr=function(t,y,d,c0){ 
 o=c(1,1) 
 o[1]=(c0*S0(t)+sum(y>t))/(c0+length(y))
 o[2]=o[1]*(c0*S0(t)+sum(y>t)+1)/(c0+length(y)+1)
 c=y[d==0]; u=c(sort(unique(c)),t+1); j=1
 while(u[j]<=t){
  d1=c0*S0(u[j])+sum(y>=u[j])
  d2=d1/(d1-sum(c==u[j]))
  o[1]=o[1]*d2
  o[2]=o[2]*d2*(d1+1)/(d1+1-sum(c==u[j]))
  j=j+1
 }
 o[2]=sqrt(o[2]-o[1]^2)
 o
}

plot(y,rep(0,length(y)),ylim=c(0,1),ylab="Weeks",xlab="Aneuploid Survival")
points(y[d==1],rep(0,length(y[d==1])),pch=16)
 
t=seq(0,400,length=100); s=matrix(0,nrow=100,ncol=2)
# with uniform prior, DP all alphas as 1
for(i in 1:100){s[i,]=svr(t[i],y,d,1)}

lines(t,s[,1],lty=1)
lines(t,s[,1]-1.28*s[,2],lty=3) # 80% CI
lines(t,s[,1]+1.28*s[,2],lty=3)
lines(survfit(Surv(y,d)~1)) # K-M estimate

