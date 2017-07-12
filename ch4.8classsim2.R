#	computed examples for chapter 4
# classification
# LDA, logistic regression, k-NN
# CPE/average predictive classification error at time (n); sum over n will give 
# cumulative average predictive error (where average is over runs)
# scenario 2
# 2) clusters separated by x2 = x1^2 and 
#	i. only x1 and x2 in model

## 
#seed<-1054
#seed<-1099
#seed<-1215
seed<-2115

set.seed(seed)
library(MASS)
library(aod)
library(ggplot2)
library(class)
library(kknn)
library(car)

# burn-in? last step? last step cannot be greater than N-1
bn<-11
# number of runs/idependent data sets
nruns<-25
# power of error
powerr<-1

# n is multiples of 10, i.e., seq(20,100,10)
nval<-seq(20,200,10)

## for saving results

CPE<-list(lda=matrix(0,nruns,length(nval)),logc=matrix(0,nruns,length(nval)),knn=matrix(0,nruns,length(nval)),kknn=matrix(0,nruns,length(nval)))
# LDA parameters
# logc parameters
# knn parameters
K=3
# kknn parameters
K=3
Dist=2
Kernel="optimal"
# define global storage
err<-vector("list",length(nval))
predv<-err
fval<-err
# create storage
for(i in 1:length(nval)){
	lst<-nval[i]
# storage for fitted values, predicted values, predictive errors, CPE for each method
	err[[i]]<-list(lda=matrix(0,nruns,lst-bn),logc=matrix(0,nruns,lst-bn),knn=matrix(0,nruns,lst-bn),kknn=matrix(0,nruns,lst-bn))
	predv[[i]]<-list(lda=matrix(99,nruns,lst-bn),logc=matrix(0,nruns,lst-bn),knn=matrix(0,nruns,lst-bn),kknn=matrix(0,nruns,lst-bn))
	fval[[i]]<-list(lda=matrix(0,nruns,lst-1),logc=matrix(0,nruns,lst-1),knn=matrix(0,nruns,lst-1),kknn=matrix(0,nruns,lst-1))
}
# start outer loop	
t1<-proc.time()
for(j in 1:nruns){
# for each run ...
# generate dataset
	sigma <- 1.0*diag(2)
	M <- t(chol(sigma))
	Z <- matrix(rnorm(nval[length(nval)]),2,nval[length(nval)]) # 2 rows, nval[length(nval)] columns
   	X <- t(M %*% Z)
   	x12<-X[,1]^2
   	x22<-X[,2]^2
   	X<-cbind(X,x12,x22)
   	y<-rep(2,dim(X)[1])
   	for(i in 1:dim(X)[1]){
		if(X[i,2]>(X[i,3]-1.1)) y[i]<-1		
	}
# incorporate some error
   	for(i in 1:dim(X)[1]){
		if(i%%10 < 1) y[i]<-(y[i]%%2)+1		
	}
	colnames(X)<-c('x1','x2','x12','x22')
	x<-as.data.frame(X)
   	Y <- as.factor(y)
	for(i in 1:length(nval)){
		lst2<-nval[i]
## LDA
# for each run do CPE ...
		kk<-bn
		for(kk in bn:(lst2-1)){
			clda <- lda(Y ~ x1+x2, data=x, prior = c(1,1)/2, subset = c(1:kk))
			predv[[i]]$lda[j,kk-bn+1]<-predict(clda,x[kk+1,],method="predictive")$class
			pr<-as.numeric(predv[[i]]$lda[j,kk-bn+1])
			err[[i]]$lda[j,kk-bn+1]<-(abs(pr-y[kk+1]))^powerr
		}
		fval[[i]]$lda[j,]<-predict(clda, x[1:kk,])$class
		CPE$lda[j,i]<-sum(err[[i]]$lda[j,])

## logistic regression
# for each run do CPE ...
		kk<-bn
		y2<-y-1
		Y2<-as.factor(y2)
		for(kk in bn:(lst2-1)){
			logc<-glm(y2 ~ x1+x2, data=x, family="binomial", subset = c(1:kk))
			# get predicted class from predicted probability
			predv[[i]]$logc[j,kk-bn+1]<-round(predict(logc, x[kk+1,], type="response"))
			err[[i]]$logc[j,kk-bn+1]<-(abs(predv[[i]]$logc[j,kk-bn+1]-y2[kk+1]))^powerr			
		}
		fval[[i]]$logc[j,]<-round(predict(logc, x[1:kk,], type="response"))
		CPE$logc[j,i]<-sum(err[[i]]$logc[j,])

## k-NN
# for each run do CPE ...
		kk<-bn
		for(kk in bn:(lst2-1)){
			cknn<-knn(x[1:kk,1:2],x[1:kk,1:2], cl=Y[1:kk], k=K)
			predv[[i]]$knn[j,kk-bn+1]<-knn(x[1:kk,],x[kk+1,], cl=Y[1:kk], k=K)
			pr<-as.numeric(predv[[i]]$knn[j,kk-bn+1])			
			err[[i]]$knn[j,kk-bn+1]<-(abs(pr-y[kk+1]))^powerr			
		}
		fval[[i]]$knn[j,]<-cknn
		CPE$knn[j,i]<-sum(err[[i]]$knn[j,])

## kknn
# for each run do CPE ...
		kk<-bn
		Xsub<-x[1:kk,]
		Ysub<-Y[1:kk]
		for(kk in bn:(lst2-1)){
			ckknn<-kknn(Ysub ~ x1+x2, train=Xsub, test=x[kk+1,], k=K, distance=Dist, kernel=Kernel)
			predv[[i]]$kknn[j,kk-bn+1]<-fitted.values(ckknn)
			err[[i]]$kknn[j,kk-bn+1]<-(abs(predv[[i]]$kknn[j,kk-bn+1]-y[kk+1]))^powerr			
		}
		fval[[i]]$kknn[j,]<-fitted.values(kknn(Ysub ~ x1+x2, train=Xsub, test=x[1:kk,], k=K, distance=Dist, kernel=Kernel))
		CPE$kknn[j,i]<-sum(err[[i]]$kknn[j,])
	}
}
# for timing
t2<-proc.time()
t2-t1

# save workspace
filen<-paste("ch4computing_cld2m1_",seed,".RData",sep="")
save.image(file=filen)

# load and compile results across saved workspaces
# key components are predv, err, CPE
.libPaths("/Library/Frameworks/R.framework/Versions/Current/library")
load('ch4computing_cld2m1_1054.RData')
predv1<-predv
err1<-err
fval1<-fval
CPE1<-CPE
load('ch4computing_cld2m1_1099.RData')
predv2<-predv
err2<-err
fval2<-fval
CPE2<-CPE
load('ch4computing_cld2m1_1215.RData')
predv3<-predv
err3<-err
fval3<-fval
CPE3<-CPE
load('ch4computing_cld2m1_2115.RData')
predv4<-predv
err4<-err
fval4<-fval
CPE4<-CPE
#
nrtot<-100
nr<-25
nval<-seq(20,200,10)
bn<-11
err0<-vector("list",length(nval))
predv0<-err0
fval0<-err0
CPE0<-list(lda=matrix(0,nrtot,length(nval)),logc=matrix(0,nrtot,length(nval)),knn=matrix(0,nrtot,length(nval)),kknn=matrix(0,nrtot,length(nval)),blrg=matrix(0,nrtot,length(nval)),fqr=matrix(0,nrtot,length(nval)),bqr=matrix(0,nrtot,length(nval)))
for(i in 1:length(nval)){
	lst<-nval[i]
	err0[[i]]<-list(lda=matrix(0,nrtot,lst-bn),logc=matrix(0,nrtot,lst-bn),knn=matrix(0,nrtot,lst-bn),kknn=matrix(0,nrtot,lst-bn))
	predv0[[i]]<-list(lda=matrix(99,nrtot,lst-bn),logc=matrix(0,nrtot,lst-bn),knn=matrix(0,nrtot,lst-bn),kknn=matrix(0,nrtot,lst-bn))
	fval0[[i]]<-list(lda=matrix(0,nrtot,lst-1),logc=matrix(0,nrtot,lst-1),knn=matrix(0,nrtot,lst-1),kknn=matrix(0,nrtot,lst-1))
}
# compile all results for all methods
for(k in 1:4){
	CPE0$lda[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('CPE',k,'$lda[1:',nr,',]',sep="")))
	CPE0$logc[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('CPE',k,'$logc[1:',nr,',]',sep="")))
	CPE0$knn[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('CPE',k,'$knn[1:',nr,',]',sep="")))
	CPE0$kknn[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('CPE',k,'$kknn[1:',nr,',]',sep="")))
}
for(k in 1:4){
	for(i in 1:length(nval)){
		err0[[i]]$lda[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('err',k,'[[',i,']]$lda[1:',nr,',]',sep="")))
		err0[[i]]$logc[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('err',k,'[[',i,']]$logc[1:',nr,',]',sep="")))
		err0[[i]]$knn[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('err',k,'[[',i,']]$knn[1:',nr,',]',sep="")))
		err0[[i]]$kknn[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('err',k,'[[',i,']]$kknn[1:',nr,',]',sep="")))
		predv0[[i]]$lda[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('predv',k,'[[',i,']]$lda[1:',nr,',]',sep="")))
		predv0[[i]]$logc[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('predv',k,'[[',i,']]$logc[1:',nr,',]',sep="")))
		predv0[[i]]$knn[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('predv',k,'[[',i,']]$knn[1:',nr,',]',sep="")))
		predv0[[i]]$kknn[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('predv',k,'[[',i,']]$kknn[1:',nr,',]',sep="")))
		fval0[[i]]$lda[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('fval',k,'[[',i,']]$lda[1:',nr,',]',sep="")))
		fval0[[i]]$logc[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('fval',k,'[[',i,']]$logc[1:',nr,',]',sep="")))
		fval0[[i]]$knn[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('fval',k,'[[',i,']]$knn[1:',nr,',]',sep="")))
		fval0[[i]]$kknn[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('fval',k,'[[',i,']]$kknn[1:',nr,',]',sep="")))
	}
}	
#
# save workspace as ch4class_sim_cld2m1.RData
#
.libPaths("/Library/Frameworks/R.framework/Versions/Current/library")
load('ch4class_sim_cld2m1.RData')
CPEres<-rbind(apply(CPE0$lda,2,mean),apply(CPE0$logc,2,mean),apply(CPE0$knn,2,mean),apply(CPE0$kknn,2,mean))
CPEresVar<-rbind(apply(CPE0$lda,2,var),apply(CPE0$logc,2,var),apply(CPE0$knn,2,var),apply(CPE0$kknn,2,var))
CPEresSE<-sqrt(t(t(CPEresVar)/nval))

# make plots of results
plot(seq(20,200,10),CPEres[1,],type='b',pch="1",cex.lab=1.2, cex.axis=1.2,ylim=c(0,61),xlab="Sample size", ylab="CPE", main="CPE by model type")
lines(seq(20,200,10),CPEres[2,],type='b',pch="2")
lines(seq(20,200,10),CPEres[3,],type='b',pch="3")
lines(seq(20,200,10),CPEres[4,],type='b',pch="4")
legend('topleft',legend=c('lda','logc','knn','kknn'),pch=c("1","2","3","4"),lty=1,cex=1.3)
# save as ch4class_cd2m1_CPE_bw_v2.pdf

CPEres1<-CPEres
for(i in 1:dim(CPEres)[1]){ CPEres1[i,]<-CPEres[i,]-CPEres[i,1]}
plot(seq(20,200,10),CPEres1[1,],type='b',pch="1",cex.lab=1.2, cex.axis=1.2,ylim=c(0,58),xlab="Sample size", ylab="CPE", main="CPE by model type (baseline adjusted)")
lines(seq(20,200,10),CPEres1[2,],type='b',pch="2")
lines(seq(20,200,10),CPEres1[3,],type='b',pch="3")
lines(seq(20,200,10),CPEres1[4,],type='b',pch="4")
legend('topleft',legend=c('lda','logc','knn','kknn'),pch=c("1","2","3","4"),lty=1,cex=1.3)
# save as ch4class_cd2m1_CPEadj_bw_v2.pdf

plot(seq(20,200,10),CPEres[1,],type='n',pch="1",cex.lab=1.2, cex.axis=1.2,ylim=c(0,61),xlab="Sample size", ylab="CPE", main="CPE by model type")
polygon(c(rev(nval),nval),c(rev(CPEres[1,]+CPEresSE[1,]),CPEres[1,]-CPEresSE[1,]),col='grey15',border=NA)
lines(seq(20,200,10),CPEres[1,])
points(seq(20,200,10),CPEres[1,],pch="1")
lines(nval,CPEres[1,]+CPEresSE[1,],lty='dashed',col='red')
lines(nval,CPEres[1,]-CPEresSE[1,],lty='dashed',col='red')
# col2rgb('grey35') is (89,89,89)
polygon(c(rev(nval),nval),c(rev(CPEres[2,]+CPEresSE[2,]),CPEres[2,]-CPEresSE[2,]),col=rgb(89/255,89/255,89/255,0.8),border=NA)
lines(seq(20,200,10),CPEres[2,])
points(seq(20,200,10),CPEres[2,],pch="2")
lines(nval,CPEres[2,]+CPEresSE[2,],lty='dashed',col='red')
lines(nval,CPEres[2,]-CPEresSE[2,],lty='dashed',col='red')
# col2rgb('grey60') is (153,153,153)
polygon(c(rev(nval),nval),c(rev(CPEres[3,]+CPEresSE[3,]),CPEres[3,]-CPEresSE[3,]),col=rgb(153/255,153/255,153/255,0.8),border=NA)
lines(seq(20,200,10),CPEres[3,])
points(seq(20,200,10),CPEres[3,],pch="3")
lines(nval,CPEres[3,]+CPEresSE[3,],lty='dashed',col='red')
lines(nval,CPEres[3,]-CPEresSE[3,],lty='dashed',col='red')
# col2rgb('grey90') is (229,229,229)
polygon(c(rev(nval),nval),c(rev(CPEres[4,]+CPEresSE[4,]),CPEres[4,]-CPEresSE[4,]),col=rgb(229/255,229/255,229/255,0.8),border=NA)
lines(seq(20,200,10),CPEres[4,])
points(seq(20,200,10),CPEres[4,],pch="4")
lines(nval,CPEres[4,]+CPEresSE[4,],lty='dashed',col='red')
lines(nval,CPEres[4,]-CPEresSE[4,],lty='dashed',col='red')
legend('topleft',legend=c('lda','logc','knn','kknn'),pch=c("1","2","3","4"),lty=1,cex=1.3)
# save as ch4class_cd2m1_CPEerr_v2.pdf




