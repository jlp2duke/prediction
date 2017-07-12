%			ch8 computing

data<-read.csv('soil_water_content.csv', row.names=NULL, header=T)
# plot data
library(lattice)
shade.col.fun <- trellis.par.get("shade.colors")$palette
shade.colors <- shade.col.fun(0.5, 0.5, seq(0, 1, length = 100))
wireframe(SWC~X*Y, data=data, xlab="X Coordinate (meters)", ylab="Y Coordinate (meters)", main="Soil Water Content", shade = TRUE, colorkey = list(col = shade.colors, at = do.breaks(range(data$SWC), 100)))
# save as trent_SWC.pdf
wireframe(Z~X*Y, data=data, xlab="X Coordinate (meters)", ylab="Y Coordinate (meters)", main="Elevation (feet above sea level)", shade = TRUE, colorkey = list(col = shade.colors, at = do.breaks(range(data$Z), 100)))
# save as trent_elevation.pdf

#orthonormal expansion (various bases); fourier, wavelets, chebyshev
#			code from PWM paper for cheby, sin, trig
#			treelet package; wavelets and wavethresh
#			http://www.di.fc.ul.pt/~jpn/r/fourier/fourier.html
			
names(data)
#[1] "SWC"    "X"      "Y"      "Z"      "WI"    
#[6] "ER_top" "ER_bot"
set.seed(31218)
# number of replicates
N<-50
# number of observations for training set
ntrain<-2^7
# number of observations to predict
npred<-2^7
# size of prediction steps, in terms of number of obs
predstep<-4
# number of explanatory variables
nex<-6
## wavelets
library(wavethresh)
# storage for results
predw<-matrix(0,N,npred)
errw<-matrix(0,N,npred)
anovaw<-vector(mode="list", length=N)
for(ii in 1:N) anovaw[[ii]]<-vector(mode="list", length=(npred/predstep))
# labels for explanatory variables
lab<-c("X","Y","Z","A","B","C")
# unique X values; for use in selecting data subset datan
xuniq<-unique(data$X)
j<-1
for(j in 1:N){
	# choose a reasonable subset of the data with which to work, n=256
	# select 512 observations so can predict from wavelet functions
	datan<-NULL
	for(i in 1:((2^9)/8)){
		tmp<-which(data$X==xuniq[i])
		tmp2<-sample(tmp,8)
		datan<-rbind(datan,data[tmp2,])
	}
	#	train on first 2^7 and predict next 2^7
	# build explanatory variables
	k<-1
	kk<-0
	while(k<npred){
		kk<-kk+1
		st<-k
		end<-ntrain+k-1
		fmla0<-"response ~"
		datat<-NULL
		jj<-1
		for(jj in 1:nex){
			namew<-paste("wpst",lab[jj], sep="")
			a<-paste(namew,'<-makewpstRO(datan[st:end,(jj+1)], datan[st:end,1], family="DaubLeAsymm")',sep="")
			eval(parse(text=a))
			a<-paste('l<-length(names(',namew,'$df))',sep="")
			eval(parse(text=a))
 			a<-paste('names(',namew,'$df)<-c("response",paste(lab[jj], 1:(l-1), sep=""))',sep="")
			eval(parse(text=a))
 			xnam1 <- paste(lab[jj], 1:(l-1), sep="")
			if(jj<nex) fmla0 <-paste(fmla0,paste(xnam1, collapse= "+"),"+")
			if(jj==nex) fmla0<-paste(fmla0,paste(xnam1, collapse= "+"))
			if(jj==1){ 
				a<-paste('datat<-',namew,'$df',sep="")
				eval(parse(text=a))
				}
			if(jj>1){
				a<-paste('datat<-cbind(datat,',namew,'$df[,2:l])',sep="")
				eval(parse(text=a))
			}
		}
		fmla1<-as.formula(fmla0)
		datax.lm1<-lm(fmla1, data=datat)
		# select variables to keep
		keep<-rownames(anova(datax.lm1))[which(anova(datax.lm1)[,5]<0.05)]
		fmla2<-as.formula(paste("response ~",paste(keep, collapse= "+")))
		datax.lm2<-lm(fmla2, data=datat)
		anovaw[[j]][[kk]]<-anova(datax.lm2)
		# predict next responses
		datap<-NULL
		jj<-1
		for(jj in 1:nex){
			namew<-paste("wpstp",lab[jj], sep="")
			a<-paste(namew,'<-makewpstRO(datan[end+1:npred,(jj+1)], datan[end+1:npred,1], family="DaubLeAsymm")',sep="")
			eval(parse(text=a))
			a<-paste('l<-length(names(',namew,'$df))',sep="")
			eval(parse(text=a))
 			a<-paste('names(',namew,'$df)<-c("response",paste(lab[jj], 1:(l-1), sep=""))',sep="")
			eval(parse(text=a))
 			xnam1 <- paste(lab[jj], 1:(l-1), sep="")
			if(jj<nex) fmla0 <-paste(fmla0,paste(xnam1, collapse= "+"),"+")
			if(jj==nex) fmla0<-paste(fmla0,paste(xnam1, collapse= "+"))
			if(jj==1){ 
				a<-paste('datap<-',namew,'$df',sep="")
				eval(parse(text=a))
				}
			if(jj>1){
				a<-paste('datap<-cbind(datap,',namew,'$df[,2:l])',sep="")
				eval(parse(text=a))
			}
		}
		pred<-predict(datax.lm2, newdata=datap)
		predw[j,(k-1)+1:predstep]<-pred[1:predstep]
		errw[j,(k-1)+1:predstep]<-datan$SWC[end+1:predstep]-predw[j,(k-1)+1:predstep]
		k<-k+predstep
	}
}

## chebyshev
library(orthopolynom)
set.seed(31218)
# storage for results
predc<-matrix(0,N,npred)
errc<-matrix(0,N,npred)
anovac<-vector(mode="list", length=N)
for(ii in 1:N) anovac[[ii]]<-vector(mode="list", length=(npred/predstep))
# labels for explanatory variables
lab<-c("X","Y","Z","A","B","C")
# unique X values; for use in selecting data subset datan
xuniq<-unique(data$X)
nlev<-10
j<-1
for(j in 1:N){
	# choose a reasonable subset of the data with which to work, n=256
	# select 512 observations so can predict from wavelet functions
	datan<-NULL
	for(i in 1:((2^9)/8)){
		tmp<-which(data$X==xuniq[i])
		tmp2<-sample(tmp,8)
		datan<-rbind(datan,data[tmp2,])
	}
	#	train on first 2^7 and predict next 2^7
	# build explanatory variables
	k<-1
	kk<-0
	while(k<npred){
		kk<-kk+1
		st<-k
		end<-ntrain+k-1
		fmla0<-"response ~"
		datat<-NULL
		jj<-1
		for(jj in 1:nex){
			namec<-paste("cheb",lab[jj], sep="")
			tmp<-polynomial.values(chebyshev.t.polynomials(nlev),datan[st:end,(jj+1)])
			a<-paste(namec,'<-as.data.frame(matrix(unlist(tmp),length(st:end),nlev+1))',sep="")
			eval(parse(text=a))
 			a<-paste('names(',namec,')<-c(paste(lab[jj], 1:(nlev+1), sep=""))',sep="")
			eval(parse(text=a))
 			xnam1 <- paste(lab[jj], 1:(nlev+1), sep="")
			if(jj<nex) fmla0 <-paste(fmla0,paste(xnam1, collapse= "+"),"+")
			if(jj==nex) fmla0<-paste(fmla0,paste(xnam1, collapse= "+"))
			if(jj==1){ 
				a<-paste('datat<-',namec,sep="")
				eval(parse(text=a))
				}
			if(jj>1){
				a<-paste('datat<-cbind(datat,',namec,')',sep="")
				eval(parse(text=a))
			}
		}
		datat<-cbind(datan[st:end,1],datat)
		names(datat)[1]<-c('response')
		fmla1<-as.formula(fmla0)
		datax.lm1<-lm(fmla1, data=datat)
		# select variables to keep
		keep<-rownames(anova(datax.lm1))[which(anova(datax.lm1)[,5]<0.05)]
		fmla2<-as.formula(paste("response ~",paste(keep, collapse= "+")))
		datax.lm2<-lm(fmla2, data=datat)
		anovac[[j]][[kk]]<-anova(datax.lm2)
		# predict next responses
		datap<-NULL
		jj<-1
		for(jj in 1:nex){
			namec<-paste("chebp",lab[jj], sep="")
			tmp<-polynomial.values(chebyshev.t.polynomials(nlev),datan[end+1:predstep,(jj+1)])
			a<-paste(namec,'<-as.data.frame(matrix(unlist(tmp),length(end+1:predstep),nlev+1))',sep="")
			eval(parse(text=a))
 			a<-paste('names(',namec,')<-c(paste(lab[jj], 1:(nlev+1), sep=""))',sep="")
			eval(parse(text=a))
 			xnam1 <- paste(lab[jj], 1:(nlev+1), sep="")
			if(jj<nex) fmla0 <-paste(fmla0,paste(xnam1, collapse= "+"),"+")
			if(jj==nex) fmla0<-paste(fmla0,paste(xnam1, collapse= "+"))
			if(jj==1){ 
				a<-paste('datap<-',namec,sep="")
				eval(parse(text=a))
				}
			if(jj>1){
				a<-paste('datap<-cbind(datap,',namec,')',sep="")
				eval(parse(text=a))
			}
		}
		pred<-predict(datax.lm2, newdata=datap)
		predc[j,(k-1)+1:predstep]<-pred
		errc[j,(k-1)+1:predstep]<-datan$SWC[end+1:predstep]-predc[j,(k-1)+1:predstep]
		k<-k+predstep
	}
}


## fourier
library(fda)
set.seed(31218)
# storage for results
predf<-matrix(0,N,npred)
errf<-matrix(0,N,npred)
anovaf<-vector(mode="list", length=N)
for(ii in 1:N) anovaf[[ii]]<-vector(mode="list", length=(npred/predstep))
# labels for explanatory variables
lab<-c("X","Y","Z","A","B","C")
# unique X values; for use in selecting data subset datan
xuniq<-unique(data$X)
nlev<-10
j<-1
for(j in 1:N){
	# choose a reasonable subset of the data with which to work, n=256
	# select 512 observations so can predict from wavelet functions
	datan<-NULL
	for(i in 1:((2^9)/8)){
		tmp<-which(data$X==xuniq[i])
		tmp2<-sample(tmp,8)
		datan<-rbind(datan,data[tmp2,])
	}
	#	train on first 2^7 and predict next 2^7
	# build explanatory variables
	k<-1
	kk<-0
	while(k<npred){
		kk<-kk+1
		st<-k
		end<-ntrain+k-1
		fmla0<-"response ~"
		datat<-NULL
		jj<-1
		for(jj in 1:nex){
			namef<-paste("four",lab[jj], sep="")
			tmp<-fourier(datan[st:end,(jj+1)],nbasis=nlev)
			a<-paste(namef,'<-as.data.frame(tmp)',sep="")
			eval(parse(text=a))
 			a<-paste('names(',namef,')<-c(paste(lab[jj], 1:(nlev+1), sep=""))',sep="")
			eval(parse(text=a))
 			xnam1 <- paste(lab[jj], 1:(nlev+1), sep="")
			if(jj<nex) fmla0 <-paste(fmla0,paste(xnam1, collapse= "+"),"+")
			if(jj==nex) fmla0<-paste(fmla0,paste(xnam1, collapse= "+"))
			if(jj==1){ 
				a<-paste('datat<-',namef,sep="")
				eval(parse(text=a))
				}
			if(jj>1){
				a<-paste('datat<-cbind(datat,',namef,')',sep="")
				eval(parse(text=a))
			}
		}
		datat<-cbind(datan[st:end,1],datat)
		names(datat)[1]<-c('response')
		fmla1<-as.formula(fmla0)
		datax.lm1<-lm(fmla1, data=datat)
		# select variables to keep
		keep<-rownames(anova(datax.lm1))[which(anova(datax.lm1)[,5]<0.05)]
		fmla2<-as.formula(paste("response ~",paste(keep, collapse= "+")))
		datax.lm2<-lm(fmla2, data=datat)
		anovaf[[j]][[kk]]<-anova(datax.lm2)
		# predict next responses
		datap<-NULL
		jj<-1
		for(jj in 1:nex){
			namef<-paste("fourp",lab[jj], sep="")
			tmp<-fourier(datan[end+1:npred,(jj+1)],nbasis=nlev)
			a<-paste(namef,'<-as.data.frame(tmp)',sep="")
			eval(parse(text=a))
 			a<-paste('names(',namef,')<-c(paste(lab[jj], 1:(nlev+1), sep=""))',sep="")
			eval(parse(text=a))
  			xnam1 <- paste(lab[jj], 1:(nlev+1), sep="")
			if(jj<nex) fmla0 <-paste(fmla0,paste(xnam1, collapse= "+"),"+")
			if(jj==nex) fmla0<-paste(fmla0,paste(xnam1, collapse= "+"))
			if(jj==1){ 
				a<-paste('datap<-',namef,sep="")
				eval(parse(text=a))
				}
			if(jj>1){
				a<-paste('datap<-cbind(datap,',namef,')',sep="")
				eval(parse(text=a))
			}
		}
		pred<-predict(datax.lm2, newdata=datap)
		predf[j,(k-1)+1:predstep]<-pred[1:predstep]
		errf[j,(k-1)+1:predstep]<-datan$SWC[end+1:predstep]-predf[j,(k-1)+1:predstep]
		k<-k+predstep
	}
}

#kernel methods	density estimation converted to point predictor
#			http://stackoverflow.com/questions/14551610/getting-values-from-kernel-density-estimation-in-r
#			https://stat.ethz.ch/R-manual/R-patched/library/stats/html/density.html
#			R packages ASH, KernSmooth
#			Choice of lambda in Kernel methods:
#					AMISE or CV

library(kedd)
set.seed(31218)
# storage for results
predk1<-matrix(0,N,npred)
errk1<-matrix(0,N,npred)
anovak1<-vector(mode="list", length=N)
for(ii in 1:N) anovak1[[ii]]<-vector(mode="list", length=(npred/predstep))
predk2<-matrix(0,N,npred)
errk2<-matrix(0,N,npred)
anovak2<-vector(mode="list", length=N)
for(ii in 1:N) anovak2[[ii]]<-vector(mode="list", length=(npred/predstep))
# unique X values; for use in selecting data subset datan
xuniq<-unique(data$X)
j<-1
for(j in 1:N){
	# choose a reasonable subset of the data with which to work, n=256
	# select 512 observations so can predict from wavelet functions
	datan<-NULL
	for(i in 1:((2^9)/8)){
		tmp<-which(data$X==xuniq[i])
		tmp2<-sample(tmp,8)
		datan<-rbind(datan,data[tmp2,])
	}
	#	train on first 2^7 and predict next 2^7
	# find kernel bandwidth by AMISE
	k<-1
	kk<-0
	while(k<npred){
		kk<-kk+1
		st<-k
		end<-ntrain+k-1
		k1<-h.amise(datan[st:end,1],kernel="epanechnikov")
		anovak1[[j]][[kk]]<-k1
		# predict next observation
		tmp<-density(datan[st:end,1], bw=k1$h, kernel="epanechnikov")
		predk1[j,(k-1)+1:predstep]<-rep(mean(tmp$x), predstep)
		errk1[j,(k-1)+1:predstep]<-datan$SWC[end+1:predstep]-predk1[j,(k-1)+1:predstep]
		k<-k+predstep
	}
	# find kernel bandwidth by ucv
	k<-1
	kk<-0
	while(k<npred){
		kk<-kk+1
		st<-k
		end<-ntrain+k-1
		k2<-h.ucv(datan[st:end,1],kernel="epanechnikov")
		anovak2[[j]][[kk]]<-k2
		# predict next observation
		tmp<-density(datan[st:end,1], bw=k2$h, kernel="epanechnikov")
		predk2[j,(k-1)+1:predstep]<-rep(mean(tmp$x), predstep)
		errk2[j,(k-1)+1:predstep]<-datan$SWC[end+1:predstep]-predk2[j,(k-1)+1:predstep]
		k<-k+predstep
	}	
}

#k-NN			density estimation converted to point predictor
#			Choice of k for NNs:
#					AMISE or CV
					
# density estimation, GCV and CV?
library(fda.usc)
set.seed(31218)
 # storage for results
predknn1<-matrix(0,N,npred)
errknn1<-matrix(0,N,npred)
anovaknn1<-vector(mode="list", length=N)
for(ii in 1:N) anovaknn1[[ii]]<-vector(mode="list", length=(npred/predstep))
predknn2<-matrix(0,N,npred)
errknn2<-matrix(0,N,npred)
anovaknn2<-vector(mode="list", length=N)
for(ii in 1:N) anovaknn2[[ii]]<-vector(mode="list", length=(npred/predstep))
# unique X values; for use in selecting data subset datan
xuniq<-unique(data$X)
j<-1
h<-seq(5,50,5)
for(j in 1:N){
	# choose a reasonable subset of the data with which to work, n=256
	# select 512 observations so can predict from wavelet functions
	datan<-NULL
	for(i in 1:((2^9)/8)){
		tmp<-which(data$X==xuniq[i])
		tmp2<-sample(tmp,8)
		datan<-rbind(datan,data[tmp2,])
	}
	#	train on first 2^7 and predict next 2^7
	# find kernel bandwidth by gcv
	k<-1
	kk<-0
	while(k<npred){
		kk<-kk+1
		st<-k
		end<-ntrain+k-1
		tmp<-fdata(datan$SWC[st:end])
		out1<-min.np(tmp, type.S=S.KNN, h=h, Ker=Ker.unif)
		fit<-out1$fdata.est$data
		anovaknn1[[j]][[kk]]<-out1
		# predict next observation
		predknn1[j,(k-1)+1:predstep]<-rep(mean(fit), predstep)
		errknn1[j,(k-1)+1:predstep]<-datan$SWC[end+1:predstep]-predknn1[j,(k-1)+1:predstep]
		k<-k+predstep
	}
	# find kernel bandwidth by cv
	k<-1
	kk<-0
	while(k<npred){
		kk<-kk+1
		st<-k
		end<-ntrain+k-1
		tmp<-fdata(datan$SWC[st:end])
		out1<-min.np(tmp, type.CV=CV.S, type.S=S.KNN, h=h, Ker=Ker.unif)
		fit<-out1$fdata.est$data
		anovaknn2[[j]][[kk]]<-out1
		# predict next observation
		predknn2[j,(k-1)+1:predstep]<-rep(mean(fit), predstep)
		errknn2[j,(k-1)+1:predstep]<-datan$SWC[end+1:predstep]-predknn2[j,(k-1)+1:predstep]
		k<-k+predstep
	}	
}

#k-NN		regression
#			knn function in package class
#			kknn package
#			ipred package function ipredknn
#			fnn package function knn.reg
#			Choice of k for NNs:
#					AMISE or CV

## NW estimator
library(np)
set.seed(31218)
 # storage for results
prednw1<-matrix(0,N,npred)
errnw1<-matrix(0,N,npred)
anovanw1<-vector(mode="list", length=N)
for(ii in 1:N) anovanw1[[ii]]<-vector(mode="list", length=(npred/predstep))
xuniq<-unique(data$X)
j<-1
for(j in 1:N){
	# choose a reasonable subset of the data with which to work, n=256
	# select 512 observations so can predict from wavelet functions
	datan<-NULL
	for(i in 1:((2^9)/8)){
		tmp<-which(data$X==xuniq[i])
		tmp2<-sample(tmp,8)
		datan<-rbind(datan,data[tmp2,])
	}
	#	train on first 2^7 and predict next 2^7
	# find NW estimate with cv.aic
	k<-1
	kk<-0
	while(k<npred){
		kk<-kk+1
		st<-k
		end<-ntrain+k-1
		fmla2<-as.formula(paste("SWC ~",paste(names(datan)[2:dim(datan)[2]], collapse= "+")))
		datax.bw<-npregbw(fmla2, data=datan, subset=st:end, regtype="lc", ckertype="epanechnikov", bwmethod="cv.aic")
		datax.np<-npreg(bws=datax.bw)
#		# select variables to keep
#		tmp<-npsigtest(datax.np)
#		keep<-names(datan)[which(tmp$P<0.05)+1]
#		fmla2<-as.formula(paste("SWC ~",paste(keep, collapse= "+")))
#		datax.bw<-npregbw(fmla2, data=datan, subset=st:end, regtype="lc", ckertype="epanechnikov", bwmethod="cv.aic")
#		datax.np<-npreg(bws=datax.bw)
		anovanw1[[j]][[kk]]<-datax.np
		# predict next responses
		pred<-predict(datax.np, newdata=datan[end+1:predstep,])
		prednw1[j,(k-1)+1:predstep]<-pred
		errnw1[j,(k-1)+1:predstep]<-datan$SWC[end+1:predstep]-prednw1[j,(k-1)+1:predstep]
		k<-k+predstep
	}
}

## kknn package
library(kknn)
set.seed(31218)
 # storage for results
predknnr1<-matrix(0,N,npred)
errknnr1<-matrix(0,N,npred)
anovaknnr1<-vector(mode="list", length=N)
for(ii in 1:N) anovaknnr1[[ii]]<-vector(mode="list", length=(npred/predstep))
xuniq<-unique(data$X)
j<-1
for(j in 1:N){
	# choose a reasonable subset of the data with which to work, n=256
	# select 512 observations so can predict from wavelet functions
	datan<-NULL
	for(i in 1:((2^9)/8)){
		tmp<-which(data$X==xuniq[i])
		tmp2<-sample(tmp,8)
		datan<-rbind(datan,data[tmp2,])
	}
	#	train on first 2^7 and predict next 2^7
	# find knn reg estimate with cv by simulation
	k<-1
	kk<-0
	while(k<npred){
		kk<-kk+1
		st<-k
		end<-ntrain+k-1
		fmla2<-as.formula(paste("SWC ~",paste(names(datan)[2:dim(datan)[2]], collapse= "+")))
		tmp<-train.kknn(fmla2, data=datan[st:end,], kmax=50, distance=1, kernel=c("epanechnikov","optimal"))		
		# select k and kernel
		datax.kknn1<-kknn(fmla2, train=datan[st:end,], test=datan[end+1:predstep,], k=tmp$best.parameters$k, kernel=tmp$best.parameters$kernel)
		anovaknnr1[[j]][[kk]]<-datax.kknn1
		# predict next responses
		pred<-datax.kknn1$fitted.values
		predknnr1[j,(k-1)+1:predstep]<-pred
		errknnr1[j,(k-1)+1:predstep]<-datan$SWC[end+1:predstep]-predknnr1[j,(k-1)+1:predstep]
		k<-k+predstep
	}
}
	
#Bayes			polya trees (univariate density)
#			Choice of prior parameters:
#					binary splits of intervals, all betas 1/2??
#			R package DPpackage - PTprior, PTdensity, PTsampler (polya trees)
#				GPPs
#			R package http://becs.aalto.fi/en/research/bayes/gpstuff/
#			http://www.r-bloggers.com/gaussian-process-regression-with-r/
#			gp-regression-demo.r
#			Choice of covariance function for GPPs:
#					AMISE or CV
		
## polya trees
library(DPpackage)
set.seed(31218)
 # storage for results
predpt1<-matrix(0,N,npred)
errpt1<-matrix(0,N,npred)
anovapt1<-vector(mode="list", length=N)
for(ii in 1:N) anovapt1[[ii]]<-vector(mode="list", length=(npred/predstep))
xuniq<-unique(data$X)
j<-1
for(j in 1:N){
	# choose a reasonable subset of the data with which to work, n=256
	# select 512 observations so can predict from wavelet functions
	datan<-NULL
	for(i in 1:((2^9)/8)){
		tmp<-which(data$X==xuniq[i])
		tmp2<-sample(tmp,8)
		datan<-rbind(datan,data[tmp2,])
	}
	#	train on first 2^7 and predict next 2^7
	# find pt density estimate
    # MCMC parameters
    nburn <- 2000
    nsave <- 5000
    nskip <- 49
	ndisplay <- 500
	mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay,tune1=0.03,tune2=0.25,tune3=1.8)
    # Prior information; uses jeffrey's prior 
    prior<-list(a0=1,b0=0.01,M=6,m0=21,S0=100,sigma=20)
	state<-NULL
	k<-1
	kk<-0
	while(k<npred){
		kk<-kk+1
		st<-k
		end<-ntrain+k-1
      	fit1 <- PTdensity(y=datan$SWC[st:end],ngrid=1000,prior=prior,mcmc=mcmc,state=state,status=TRUE)
		# predict next responses
		pred<-fit1$x1[which(fit1$dens==max(fit1$dens))]
		predpt1[j,(k-1)+1:predstep]<-rep(pred,predstep)
		errpt1[j,(k-1)+1:predstep]<-datan$SWC[end+1:predstep]-predpt1[j,(k-1)+1:predstep]
		k<-k+predstep
		state<-fit1$state
	}
}

## GPPs by library(kernlab)
library(kernlab)
set.seed(31218)
 # storage for results
predgp2<-matrix(0,N,npred)
errgp2<-matrix(0,N,npred)
anovagp2<-vector(mode="list", length=N)
for(ii in 1:N) anovagp2[[ii]]<-vector(mode="list", length=(npred/predstep))
xuniq<-unique(data$X)
j<-1
for(j in 1:N){
	# choose a reasonable subset of the data with which to work, n=256
	# select 512 observations so can predict from wavelet functions
	datan<-NULL
	for(i in 1:((2^9)/8)){
		tmp<-which(data$X==xuniq[i])
		tmp2<-sample(tmp,8)
		datan<-rbind(datan,data[tmp2,])
	}
	#	train on first 2^7 and predict next 2^7
	# find gpr estimate
    # prior parameters
    kpar<-"automatic"
	k<-1
	kk<-0
	while(k<npred){
		kk<-kk+1
		st<-k
		end<-ntrain+k-1
      	fmla2<-as.formula(paste("response ~",paste(names(datan)[2:length(datan)], collapse= "+")))
		fit1 <- gausspr(x=datan[st:end,2:dim(datan)[2]], y=datan[st:end,1], type="regression", kernel="rbfdot", kpar=kpar)
		# predict next responses
		predgp2[j,(k-1)+1:predstep]<- predict(fit1,datan[end+1:predstep,2:dim(datan)[2]])		
		errgp2[j,(k-1)+1:predstep]<-datan$SWC[end+1:predstep]-predgp2[j,(k-1)+1:predstep]
		k<-k+predstep
	}
}
# save.image('trentgp.RData')

## compile all results for both plots and ensembling
nval<-c(ntrain+1:npred)
## wavelets ...
cpew<-t(apply(errw,1,function(x){cumsum(abs(x))}))
resw<-matrix(0,2,npred)
resw[1,]<-apply(cpew,2,median)
resw[2,]<-apply(cpew,2,var)/nval
resw[2,]<-sqrt(resw[2,])
plot(1:npred,resw[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="wavelets", cex.lab=1.4, cex.axis=1.4, cex.main=1.4)
lines(1:npred,resw[1,]+2*resw[2,],lty=2)
lines(1:npred,resw[1,]-2*resw[2,],lty=2)
# save as pred_wavelets_v2.pdf

# cheby ...
cpec<-t(apply(errc,1,function(x){cumsum(abs(x))}))
resc<-matrix(0,2,npred)
resc[1,]<-apply(cpec,2,median)
resc[2,]<-apply(cpec,2,var)/nval
resc[2,]<-sqrt(resc[2,])
plot(1:npred,resc[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="chebyshev", cex.lab=1.4, cex.axis=1.4, cex.main=1.4)
lines(1:npred,resc[1,]+2*resc[2,],lty=2)
lines(1:npred,resc[1,]-2*resc[2,],lty=2)
# save as pred_cheby_v2.pdf

# fourier ...
cpef<-t(apply(errf,1,function(x){cumsum(abs(x))}))
resf<-matrix(0,2,npred)
resf[1,]<-apply(cpef,2,median)
resf[2,]<-apply(cpef,2,var)/nval
resf[2,]<-sqrt(resf[2,])
plot(1:npred,resf[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="fourier", cex.lab=1.4, cex.axis=1.4, cex.main=1.4)
lines(1:npred,resf[1,]+2*resf[2,],lty=2)
lines(1:npred,resf[1,]-2*resf[2,],lty=2)
# save as pred_fourier_v2.pdf

# kernel AMISE
cpek1<-t(apply(errk1,1,function(x){cumsum(abs(x))}))
resk1<-matrix(0,2,npred)
resk1[1,]<-apply(cpek1,2,median)
resk1[2,]<-apply(cpek1,2,var)/nval
resk1[2,]<-sqrt(resk1[2,])
plot(1:npred,resk1[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="kernel AMISE", cex.lab=1.4, cex.axis=1.4, cex.main=1.4)
lines(1:npred,resk1[1,]+2*resk1[2,],lty=2)
lines(1:npred,resk1[1,]-2*resk1[2,],lty=2)
# save as pred_kernel_AMISE_v2.pdf

# kernel UCV
cpek2<-t(apply(errk2,1,function(x){cumsum(abs(x))}))
resk2<-matrix(0,2,npred)
resk2[1,]<-apply(cpek2,2,median)
resk2[2,]<-apply(cpek2,2,var)/nval
resk2[2,]<-sqrt(resk2[2,])
plot(1:npred,resk2[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="kernel UCV", cex.lab=1.4, cex.axis=1.4, cex.main=1.4)
lines(1:npred,resk2[1,]+2*resk2[2,],lty=2)
lines(1:npred,resk2[1,]-2*resk2[2,],lty=2)
# save as pred_kernel_UCV_v2.pdf

# k-NN density estimation by GCV
cpeknn1<-t(apply(errknn1,1,function(x){cumsum(abs(x))}))
resknn1<-matrix(0,2,npred)
resknn1[1,]<-apply(cpeknn1,2,median)
resknn1[2,]<-apply(cpeknn1,2,var)/nval
resknn1[2,]<-sqrt(resknn1[2,])
plot(1:npred,resknn1[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="k-NN density est. GCV", cex.lab=1.4, cex.axis=1.4, cex.main=1.4)
lines(1:npred,resknn1[1,]+2*resknn1[2,],lty=2)
lines(1:npred,resknn1[1,]-2*resknn1[2,],lty=2)
# save as pred_kNNdens_GCV_v2.pdf

# k-NN density estimation by CV
cpeknn2<-t(apply(errknn2,1,function(x){cumsum(abs(x))}))
resknn2<-matrix(0,2,npred)
resknn2[1,]<-apply(cpeknn2,2,median)
resknn2[2,]<-apply(cpeknn2,2,var)/nval
resknn2[2,]<-sqrt(resknn2[2,])
plot(1:npred,resknn2[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="k-NN density est. CV", cex.lab=1.4, cex.axis=1.4, cex.main=1.4)
lines(1:npred,resknn2[1,]+2*resknn2[2,],lty=2)
lines(1:npred,resknn2[1,]-2*resknn2[2,],lty=2)
# save as pred_kNNdens_CV_v2.pdf

# Nadaraya-Watson 
cpenw1<-t(apply(errnw1,1,function(x){cumsum(abs(x))}))
resnw1<-matrix(0,2,npred)
resnw1[1,]<-apply(cpenw1,2,median)
resnw1[2,]<-apply(cpenw1,2,var)/nval
resnw1[2,]<-sqrt(resnw1[2,])
plot(1:npred,resnw1[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="Nadaraya-Watson", cex.lab=1.4, cex.axis=1.4, cex.main=1.4)
lines(1:npred,resnw1[1,]+2*resnw1[2,],lty=2)
lines(1:npred,resnw1[1,]-2*resnw1[2,],lty=2)
# save as pred_nadarayawatson_v2.pdf

# k-NN regression
cpeknnr<-t(apply(errknnr1,1,function(x){cumsum(abs(x))}))
resknnr<-matrix(0,2,npred)
resknnr[1,]<-apply(cpeknnr,2,median)
resknnr[2,]<-apply(cpeknnr,2,var)/nval
resknnr[2,]<-sqrt(resknnr[2,])
plot(1:npred,resknnr[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="k-NN regression", cex.lab=1.4, cex.axis=1.4, cex.main=1.4)
lines(1:npred,resknnr[1,]+2*resknnr[2,],lty=2)
lines(1:npred,resknnr[1,]-2*resknnr[2,],lty=2)
# save as pred_kNNreg_v2.pdf

# Bayes polya trees
cpept<-t(apply(errpt1,1,function(x){cumsum(abs(x))}))
respt<-matrix(0,2,npred)
respt[1,]<-apply(cpept,2,median)
respt[2,]<-apply(cpept,2,var)/nval
respt[2,]<-sqrt(respt[2,])
plot(1:npred,respt[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="Polya Trees", cex.lab=1.4, cex.axis=1.4, cex.main=1.4)
lines(1:npred,respt[1,]+2*respt[2,],lty=2)
lines(1:npred,respt[1,]-2*respt[2,],lty=2)
# save as pred_polyatree_v2.pdf

# Gaussian process prior
cpegp2<-t(apply(errgp2,1,function(x){cumsum(abs(x))}))
resgp2<-matrix(0,2,npred)
resgp2[1,]<-apply(cpegp2,2,median)
resgp2[2,]<-apply(cpegp2,2,var)/nval
resgp2[2,]<-sqrt(resgp2[2,])
plot(1:npred,resgp2[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="GPP regression", cex.lab=1.4, cex.axis=1.4, cex.main=1.4)
lines(1:npred,resgp2[1,]+2*resgp2[2,],lty=2)
lines(1:npred,resgp2[1,]-2*resgp2[2,],lty=2)
# save as pred_GPPreg_v2.pdf

## prediction by straight model averaging
cperes<-abind(cpew,cpec,cpef,cpek1,cpek2,cpeknn1,cpeknn2,cpenw1,cpeknnr,cpept,cpegp2,along=3)
modavg<-apply(cperes,c(1,2),mean)
dim(modavg)
#[1]  50 128
resavg<-matrix(0,2,npred)
resavg[1,]<-apply(modavg,2,median)
resavg[2,]<-apply(modavg,2,var)/nval
resavg[2,]<-sqrt(modavg[2,])
plot(1:npred,resavg[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="model average", cex.lab=1.4, cex.axis=1.4, cex.main=1.4)
lines(1:npred,resavg[1,]+2*resavg[2,],lty=2)
lines(1:npred,resavg[1,]-2*resavg[2,],lty=2)
# save as pred_modavg_v2.pdf

# remove wavelets
cperes2<-abind(cpec,cpef,cpek1,cpek2,cpeknn1,cpeknn2,cpenw1,cpeknnr,cpept,cpegp2,along=3)
modavg2<-apply(cperes2,c(1,2),mean)
dim(modavg)
#[1]  50 128
resavg2<-matrix(0,2,npred)
resavg2[1,]<-apply(modavg2,2,median)
resavg2[2,]<-apply(modavg2,2,var)/nval
resavg2[2,]<-sqrt(modavg2[2,])
plot(1:npred,resavg2[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="model average (no wavelets)", cex.lab=1.4, cex.axis=1.4, cex.main=1.4)
lines(1:npred,resavg2[1,]+2*resavg2[2,],lty=2)
lines(1:npred,resavg2[1,]-2*resavg2[2,],lty=2)
# save as pred_modavg2_v2.pdf

## adaptive model selection based on most recent predicted values, starting with ntrain+1

i<-j<-1
burnin<-64
modavg3<-matrix(0,N,npred-1)
modavgmeth<-matrix(0,N,npred-1)
for(i in 1:N){
	for(j in 1:(npred-1)){
		if(j<burnin){ 
		tmp<-apply(cperes[i,1:(j+1),],2,sum)
		}
		if(j>=burnin){ 
		tmp<-apply(cperes[i,(j-burnin+1):(j+1),],2,sum)
		}
		modavg3[i,j]<-cperes[i,j+1,which(tmp==min(tmp))[1]]
		modavgmeth[i,j]<-which(tmp==min(tmp))[1]
	}
}
hist(modavgmeth, breaks=seq(0.5,11.5,1), xlab='selected method', main="frequencies of models selected")
# save as pred_opt_select.pdf
plot(1:127,rep(1,127),type="n",ylim=c(1,50),xlab="predictive step", ylab="iteration",main="method selected")
text(1:127,rep(1,127),as.character(modavgmeth[1,]),cex=0.4)
for(i in 2:50) text(1:127,rep(i,127),as.character(modavgmeth[i,]),cex=0.4)
# save as pred_opt_select_sequence.pdf
resavg3<-matrix(0,2,npred-1)
resavg3[1,]<-apply(modavg3,2,median)
resavg3[2,]<-apply(modavg3,2,var)/nval[2:npred]
resavg3[2,]<-sqrt(modavg3[2,])
plot(1:(npred-1),resavg3[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="optimal predictive selection", cex.lab=1.4, cex.axis=1.4, cex.main=1.4)
lines(1:(npred-1),resavg3[1,]+2*resavg3[2,],lty=2)
lines(1:(npred-1),resavg3[1,]-2*resavg3[2,],lty=2)
# save as pred_opt_select_v2.pdf

