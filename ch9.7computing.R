#		ch9.7 computing

#	use M-H as a method to find BMA, explore posterior on parameters and posterior on
#		model space
#	simulated annealing
#		split data into three parts, one for train - one for MSE calc - one for 
#			indep assessment
#		burn-in to estimate parameters

#BIC:	estimate parameters by something like backprop for each of 20 models, find BIC 
#		for each model, do model selection and model average with BIC as weights
#BMA/MH:	pick the best models under J, average the best models or select by CV
#SA+CV:	jumps on nodes, once a node is selected use SA to est. parameters

# linear models with SA
library(caret)
laheart<-read.csv('laheart.csv', row.names=1, header=T)
nreps<-10
nobs<-150
niter<-300
burnin<-50
npred<-nobs-burnin
predstep<-5
nvar<-16
# items for storing results
predSA<-matrix(0,nreps,npred)
errSA<-matrix(0,nreps,npred)
modSA<-vector(mode="list", length=nreps)
ii<-1
for(ii in 1:nreps) modSA[[ii]]<-vector(mode="list", length=(npred/predstep))

i<-1
set.seed(104050450)
for(i in 1:nreps){
	j<-burnin
	k<-1
	datsel<-sort(sample.int(dim(laheart)[1],size=nobs))
	pdata<-laheart[datsel,]
	sa_ctrl<-safsControl(functions=caretSA, method="boot", number=1, verbose=T, seeds=c(sample.int(1e5,2)))
	while(j<nobs){
		mod<-safs(x=pdata[1:j,-13],y=pdata$CHOL_62[1:j],iters=niter, safsControl=sa_ctrl, method="lm")
		modSA[[i]][[k]]<-mod
		predSA[i,(j-burnin+1:predstep)]<-predict(mod,newdata=pdata[j+1:predstep,])
		errSA[i,(j-burnin+1:predstep)]<-sqrt((predSA[i,(j-burnin+1:predstep)]-pdata$CHOL_62[j+1:predstep])^2)
		j<-j+predstep
		k<-k+1
	}	
}

plot(mod) + theme_bw(base_size=16)+ theme(axis.text=element_text(size=16),legend.position="none", panel.grid.major=element_line(colour="grey50"), panel.grid.minor=element_line(colour="grey50"))+geom_point(size=4,aes(colour="blue",shape=factor(Estimate)))+annotate("text",x=200,y=39,label="internal=triangle",size=8)+annotate("text",x=200,y=38,label="external=circle",size=8)
plot(modSA[[10]][[10]]) + theme_bw(base_size=16)+ theme(axis.text=element_text(size=16),legend.position="none", panel.grid.major=element_line(colour="grey50"), panel.grid.minor=element_line(colour="grey50"))+geom_point(size=4,aes(colour="blue",shape=factor(Estimate)))+annotate("text",x=200,y=59,label="internal=triangle",size=8)+annotate("text",x=200,y=58,label="external=circle",size=8)

# for error plots
errSAm<-apply(errSA,2,mean)
errSAu<-apply(errSA,2,quantile, probs=0.95)
errSAl<-apply(errSA,2,quantile, probs=0.05)
cerrSA<-t(apply(errSA,1,cumsum))
cerrSAm<-apply(cerrSA,2,mean)
cerrSAu<-apply(cerrSA,2,quantile, probs=0.95)
cerrSAl<-apply(cerrSA,2,quantile, probs=0.05)
# number of terms
termSA<-vector(mode="list", length=nreps)
ii<-1
for(ii in 1:nreps) termSA[[ii]]<-vector(mode="list", length=(npred/predstep))
niter<-300
ii<-jj<-1
for(ii in 1:nreps){
	for(jj in 1:(npred/predstep)){
		termSA[[ii]][[jj]]<-modSA[[ii]][[jj]]$sa$subsets[[modSA[[ii]][[jj]]$optIter]]
	}	
}

# test rf with no subsampling
library(caret)
laheart<-read.csv('laheart.csv', row.names=1, header=T)
nreps<-10
nobs<-150
niter<-300
burnin<-50
npred<-nobs-burnin
predstep<-5
nvar<-16
# items for storing results
predfSA<-matrix(0,nreps,npred)
errfSA<-matrix(0,nreps,npred)
modfSA<-vector(mode="list", length=nreps)
ii<-1
for(ii in 1:nreps) modfSA[[ii]]<-vector(mode="list", length=(npred/predstep))

i<-1
set.seed(104050450)
for(i in 1:nreps){
	j<-burnin
	k<-1
	datsel<-sort(sample.int(dim(laheart)[1],size=nobs))
	pdata<-laheart[datsel,]
	sa_ctrl<-safsControl(functions = rfSA, method="boot", number=1, verbose=TRUE, improve=50, seeds=c(sample.int(1e5,2)))
	while(j<nobs){
		rf_sa<-safs(x=pdata[1:j,-13],y=pdata$CHOL_62[1:j], iters=niter, safsControl=sa_ctrl, ntree=50, importance=TRUE)
		modfSA[[i]][[k]]<-rf_sa
		predfSA[i,(j-burnin+1:predstep)]<-predict(rf_sa,newdata=pdata[j+1:predstep,])
		errfSA[i,(j-burnin+1:predstep)]<-sqrt((predfSA[i,(j-burnin+1:predstep)]-pdata$CHOL_62[j+1:predstep])^2)
		j<-j+predstep
		k<-k+1
	}	
}
plot(mod) + theme_bw(base_size=16)+ theme(axis.text=element_text(size=16),legend.position="none", panel.grid.major=element_line(colour="grey50"), panel.grid.minor=element_line(colour="grey50"))+geom_point(size=4,aes(colour="blue",shape=factor(Estimate)))+annotate("text",x=200,y=39,label="internal=triangle",size=8)+annotate("text",x=200,y=38,label="external=circle",size=8)
plot(modfSA[[8]][[4]]) + theme_bw(base_size=16)+ theme(axis.text=element_text(size=16),legend.position="none", panel.grid.major=element_line(colour="grey50"), panel.grid.minor=element_line(colour="grey50"))+geom_point(size=4,aes(colour="blue",shape=factor(Estimate)))+annotate("text",x=50,y=28,label="internal=triangle",size=8)+annotate("text",x=50,y=27,label="external=circle",size=8)

# number of terms/what terms
termRF<-vector(mode="list", length=nreps)
ii<-1
for(ii in 1:nreps) termRF[[ii]]<-vector(mode="list", length=(npred/predstep))
niter<-300
ii<-jj<-1
for(ii in 1:nreps){
	for(jj in 1:(npred/predstep)){
		termRF[[ii]][[jj]]<-modfSA[[ii]][[jj]]$sa$subsets[[modfSA[[ii]][[jj]]$optIter]]
	}	
}

errSAfm<-apply(errfSA,2,mean)
errSAfu<-apply(errfSA,2,quantile, probs=0.95)
errSAfl<-apply(errfSA,2,quantile, probs=0.05)
cerrSAf<-t(apply(errfSA,1,cumsum))
cerrSAfm<-apply(cerrSAf,2,mean)
cerrSAfu<-apply(cerrSAf,2,quantile, probs=0.95)
cerrSAfl<-apply(cerrSAf,2,quantile, probs=0.05)

# test lm with MCMC

library(BMS)
laheart<-read.csv('laheart.csv', row.names=1, header=T)
nreps<-10
nobs<-150
niter<-50000
burnin<-50
nburn<-10000
nmod<-2000
npred<-nobs-burnin
predstep<-5
nvar<-16
# items for storing results
predMC<-matrix(0,nreps,npred)
errMC<-matrix(0,nreps,npred)
modMC<-vector(mode="list", length=nreps)
ii<-1
for(ii in 1:nreps) modMC[[ii]]<-vector(mode="list", length=(npred/predstep))

i<-1
set.seed(104050450)
for(i in 1:nreps){
	j<-burnin
	k<-1
	datsel<-sort(sample.int(dim(laheart)[1],size=nobs))
	pdata<-laheart[datsel,]
	while(j<nobs){
#		lm_mcmc<-bms(as.matrix(cbind(pdata$CHOL_62[1:j],pdata[1:j,-13])),burn=nburn, iter=niter, g="hyper", mprior="random", nmodel=nmod, mcmc="bd", user.int=F)
		lm_mcmc<-bms(as.matrix(cbind(pdata$CHOL_62[1:j],pdata[1:j,-13])),burn=nburn, iter=niter, g="BRIC", mprior="uniform", nmodel=nmod, mcmc="bd", user.int=F)
		modMC[[i]][[k]]<-lm_mcmc
		predMC[i,(j-burnin+1:predstep)]<-predict(lm_mcmc,newdata=pdata[j+1:predstep,], topmodels=1)
#		predMC[i,(j-burnin+1:predstep)]<-predict(lm_mcmc,newdata=pdata[j+1:predstep,])
		errMC[i,(j-burnin+1:predstep)]<-sqrt((predMC[i,(j-burnin+1:predstep)]-pdata$CHOL_62[j+1:predstep])^2)
		j<-j+predstep
		k<-k+1
	}	
}
plotConv(lm_mcmc)
plotConv(lm_mcmc[1:100])
# number of terms/what terms
termMC<-vector(mode="list", length=nreps)
ii<-1
for(ii in 1:nreps) termMC[[ii]]<-vector(mode="list", length=(npred/predstep))
niter<-300
ii<-jj<-1
for(ii in 1:nreps){
	for(jj in 1:(npred/predstep)){
		termMC[[ii]][[jj]]<-c(sum(coef(modMC[[ii]][[jj]])[,1]),coef(modMC[[ii]][[jj]],order.by.pip=F)[,1])
	}	
}

plotConv(lm_mcmc[1:100], lty=c(1,2), col="black")
image(lm_mcmc, col  = gray(c(0.25,0.75)))
plotModelsize(lm_mcmc, col="black")

plotConv(modMC[[1]][[20]])
plotConv(modMC[[1]][[20]][1:100])
errMCm<-apply(errMC,2,mean)
errMCu<-apply(errMC,2,quantile, probs=0.95)
errMCl<-apply(errMC,2,quantile, probs=0.05)
cerrMC<-t(apply(errMC,1,cumsum))
cerrMCm<-apply(cerrMC,2,mean)
cerrMCu<-apply(cerrMC,2,quantile, probs=0.95)
cerrMCl<-apply(cerrMC,2,quantile, probs=0.05)

errMCm2<-apply(errMC,2,mean)
errMCu2<-apply(errMC,2,quantile, probs=0.95)
errMCl2<-apply(errMC,2,quantile, probs=0.05)
cerrMC2<-t(apply(errMC,1,cumsum))
cerrMCm2<-apply(cerrMC2,2,mean)
cerrMCu2<-apply(cerrMC2,2,quantile, probs=0.95)
cerrMCl2<-apply(cerrMC2,2,quantile, probs=0.05)

errMCm3<-apply(errMC,2,mean)
errMCu3<-apply(errMC,2,quantile, probs=0.95)
errMCl3<-apply(errMC,2,quantile, probs=0.05)
cerrMC3<-t(apply(errMC,1,cumsum))
cerrMCm3<-apply(cerrMC3,2,mean)
cerrMCu3<-apply(cerrMC3,2,quantile, probs=0.95)
cerrMCl3<-apply(cerrMC3,2,quantile, probs=0.05)

#pmp.bma(lm_mcmc)[1:5,]
#       PMP (Exact) PMP (MCMC)
#000000  0.10751411    0.09168
#000010  0.04117801    0.03816
#000200  0.03853308    0.03548
#000008  0.03072671    0.02528
#002000  0.01872916    0.01646
# colSums(pmp.bma(lm_mcmc))
#PMP (Exact)  PMP (MCMC) 
#    0.95836     0.95836 

# stepwise absolute error (over steps of size five)
plot(x=1:100, y=errSAm, ylim=c(1,950),xaxt="n",xlab="timestep",ylab="absolute error", type="n", cex=1.4, cex.axis=1.4, cex.lab=1.4)
axis(1, at=c(0,20,40,60,80,100), labels=c("50","70","90","110","130","150"))
lines(1:100,errSAm, lty=1)
lines(1:100,errSAl, lty=1, col="grey75")
lines(1:100,errSAu, lty=1, col="grey75")
lines(1:100,errSAfm, lty=2)
lines(1:100,errSAfl, lty=2, col="grey75")
lines(1:100,errSAfu, lty=2, col="grey75")
lines(1:100,errMCm, lty=3)
lines(1:100,errMCl, lty=3, col="grey75")
lines(1:100,errMCu, lty=3, col="grey75")
lines(1:100,errMCl3, lty=4, col="grey75")
lines(1:100,errMCm3, lty=4)
lines(1:100,errMCu3, lty=4, col="grey75")
legend('topright',legend=c('lmSA','rfSA','lmMCMC-All','lmMCMC-Top'),lty=1:4,cex=1.4)

# cumulative error
plot(x=1:100, y=cerrSAm, ylim=c(5,35900),xaxt="n",xlab="timestep",ylab="cumulative absolute error", type="n", cex=1.4, cex.axis=1.4, cex.lab=1.4)
axis(1, at=c(0,20,40,60,80,100), labels=c("50","70","90","110","130","150"))
lines(1:100,cerrSAm, lty=1)
lines(1:100,cerrSAl, lty=1, col="grey75")
lines(1:100,cerrSAu, lty=1, col="grey75")
lines(1:100,cerrSAfm, lty=2)
lines(1:100,cerrSAfl, lty=2, col="grey75")
lines(1:100,cerrSAfu, lty=2, col="grey75")
lines(1:100,cerrMCm, lty=3)
lines(1:100,cerrMCl, lty=3, col="grey75")
lines(1:100,cerrMCu, lty=3, col="grey75")
lines(1:100,cerrMCl3, lty=4, col="grey75")
lines(1:100,cerrMCm3, lty=4)
lines(1:100,cerrMCu3, lty=4, col="grey75")
legend('topleft',legend=c('lmSA','rfSA','lmMCMC-All','lmMCMC-Top'),lty=1:4, cex=1.4)

# number of terms/what terms in best models
modterm<-matrix(0,3,dim(laheart)[2]-1)
modsize<-c(0,0,0)
ii<-jj<-1
for(ii in 1:nreps){
	for(jj in 1:(npred/predstep)){
		modterm[1,termSA[[ii]][[jj]]]<-modterm[1,termSA[[ii]][[jj]]]+1
		modsize[1]<-modsize[1]+length(termSA[[ii]][[jj]])
		modterm[2,termRF[[ii]][[jj]]]<-modterm[2,termRF[[ii]][[jj]]]+1
		modsize[2]<-modsize[2]+length(termRF[[ii]][[jj]])
		modterm[3,]<-modterm[3,]+termMC[[ii]][[jj]][-1]
		modsize[3]<-modsize[3]+termMC[[ii]][[jj]][1]
	}	
}
layout(matrix(c(1,2),1,2),widths=c(1,3))
par(las=2)
barplot(modsize/200,names.arg=c("lm_sa","rf_sa","lm_mcmc"),main="avg model size",cex.main=0.8,cex.axis=0.8,cex.names=0.8)
plot(x=c(1,nvar), y=c(0,1), ylim=c(0,1),xaxt="n",xlab="model term",ylab="probability of inclusion", main="model term selection",type="n",cex.main=0.8,cex.axis=0.8,cex.lab=0.8)
axis(1, at=c(1:nvar), labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"),cex.axis=0.8)
lines(1:nvar,modterm[1,]/200,lty=1)
lines(1:nvar,modterm[2,]/200,lty=2)
lines(1:nvar,modterm[3,]/200,lty=3)
legend("topright",legend=c('lm_sa','rf_sa','lm_mc'),lty=1:3,cex=0.8)
