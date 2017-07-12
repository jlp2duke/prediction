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

seed<-10010
set.seed(seed)

library(psbcGroup)
library(coda)
library(survival)
datorig<-read.csv('uissurv.csv', header=T, row.names=NULL)
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
		survDat$di<-newcens[permord[[i]][1:(nbase+j-1)]
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
		modPH[[i]]$bayes2[[j]]<-fitGL
		medPH[[i]]$bayes2[j]<-my.med
	}
}

# save workspace
filen<-paste("ch7computing_survBayes_c",seed,".RData",sep="")
save.image(file=filen)

## integrating all results across random permutations
# 20 runs
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_survBayes_c10010.RData")
FPE1bc1<-FPE
modPH1bc1<-modPH
medPH1bc1<-medPH
# 20 runs each
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_survBayes_c21012.RData")
FPE2bc1<-FPE
modPH2bc1<-modPH
medPH2bc1<-medPH
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_survBayes_c21112.RData")
FPE3bc1<-FPE
modPH3bc1<-modPH
medPH3bc1<-medPH
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_survBayes_c21122.RData")
FPE4bc1<-FPE
modPH4bc1<-modPH
medPH4bc1<-medPH
load("/Users/jclarke/Documents/predictive/Book/ch7/uissurv/ch7computing_survBayes_c31122.RData")
FPE5bc1<-FPE
modPH5bc1<-modPH
medPH5bc1<-medPH
# 
nperm2<-20
nperm<-100
#
FPEavgbc1<-c(rep(0,npred))
for(i in 1:nperm2) FPEavgbc1<-FPEavgbc1+FPE1bc1[[i]]$bayes2
for(i in 1:nperm2) FPEavgbc1<-FPEavgbc1+FPE2bc1[[i]]$bayes2
for(i in 1:nperm2) FPEavgbc1<-FPEavgbc1+FPE3bc1[[i]]$bayes2
for(i in 1:nperm2) FPEavgbc1<-FPEavgbc1+FPE4bc1[[i]]$bayes2
for(i in 1:nperm2) FPEavgbc1<-FPEavgbc1+FPE5bc1[[i]]$bayes2
FPEavgbc1<-FPEavgbc1/nperm
FPEsdbc1<-c(rep(0,npred))
for(i in 1:nperm2) FPEsdbc1<-FPEsdbc1+(FPEavgbc1-FPE1bc1[[i]]$bayes2)^2
for(i in 1:nperm2) FPEsdbc1<-FPEsdbc1+(FPEavgbc1-FPE2bc1[[i]]$bayes2)^2
for(i in 1:nperm2) FPEsdbc1<-FPEsdbc1+(FPEavgbc1-FPE3bc1[[i]]$bayes2)^2
for(i in 1:nperm2) FPEsdbc1<-FPEsdbc1+(FPEavgbc1-FPE4bc1[[i]]$bayes2)^2
for(i in 1:nperm2) FPEsdbc1<-FPEsdbc1+(FPEavgbc1-FPE5bc1[[i]]$bayes2)^2
FPEsdbc1<-sqrt(FPEsdbc1/(nperm-1))
plot(1:npred,FPE1bc1[[1]]$bayes2, type="l", xlab='time step', ylab='L1 cumulative predictive error')
lines(1:npred,FPE2bc1[[2]]$bayes2, type="l", col='red') 
lines(1:npred,FPEavgbc1, type="l", col='blue')
lines(1:npred,FPEavgbc1-2*FPEsdbc1, type="l", lty=2, col='green')
lines(1:npred,FPEavgbc1+2*FPEsdbc1, type="l", lty=2, col='green')
# bayCPEc1_n80_PH.pdf
## redo plot for Esther, bw and better labels
plot(1:npred,FPE1bc1[[1]]$bayes2, type="l", xlab='time step', ylab='L1 cumulative predictive error', col='grey75', cex.lab=1.2, cex.axis=1.2, cex=1.2, ylim=c(0,17000))
lines(1:npred,FPE2bc1[[2]]$bayes2, type="l", col='grey75') 
lines(1:npred,FPEavgbc1, type="l", col='black')
lines(1:npred,FPEavgbc1-2*FPEsdbc1, type="l", lty=2)
lines(1:npred,FPEavgbc1+2*FPEsdbc1, type="l", lty=2)
# save as bayCPE-c1-PH_v2.pdf
# remove lots of stuff, then save as ch7computing_survBayesc1.RData

