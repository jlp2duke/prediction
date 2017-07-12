#	real examples for chapter 4
# regression; no intercept
# World Bank data
# Frequentist linear models
# bayesian linear models with conjugate priors (R BLR package)
#						with noninformative priors (R BLR package)
#						with spike-and-slab priors (R spikeslab package)
#						with g-priors (R BMS package)
# frequentist quantile regression (median); \alpha/2 and 1-\alpha/2 quantiles (R quantreg)
# bayesian quantile regression (median); \alpha/2 and 1-\alpha/2 quantiles (R bayesQR and BSquare)
# point predictions and predictive intervals
# CPE/average predictive error at time (n); sum over n will give cumulative average
# predictive error (where average is over runs)
## 

set.seed(1054)
library(BLR)
library(spikeslab)
library(BMS)
library(quantreg)
#library(bayesQR)
#library(BSquare)

# year 2010
tmp<-read.csv('WorldBank_AgAndRural_1_Topic_MetaData_en_EXCEL.csv',header=T, row.names=NULL)
# reshape data into matrix format with one column per variable
wide<-reshape(tmp, idvar=names(tmp)[1], timevar=names(tmp)[2], direction="wide")
#dim(wide)
##[1]	246	33
row.names(wide)<-wide[,1]
wide<-wide[,-1]
# how much data do we really have?
obs<-c(rep(99,dim(wide)[2]))
for(i in 1:length(obs)){ 
	obs[i]<-sum(!is.na(wide[,i]))
}
#obs
## [1] 0	103  26 237 237 205 236 236 210 225 238 238  50   0 244   0 197 199 198 
##[20] 244 209 166 177 202  14  30  95 241 235 241 173 170
# select variables with at least 160 obs
sel1<-c(4:12,15,17:24,28:32)
widesel1<-wide[,sel1]
varobs<-c(rep(99,dim(wide)[1]))
for(i in 1:length(varobs)){ varobs[i]<-sum(!is.na(widesel1[i,])) }
# how many complete observations for selected variables?
#length(varobs[varobs>22])
##[1] 90
widesel2<-widesel1[varobs>22,]
#colnames(widesel2)
## [1] "Year.2010.AG.LND.AGRI.K2"    "Year.2010.AG.LND.AGRI.ZS"   
## [3] "Year.2010.AG.LND.ARBL.HA"    "Year.2010.AG.LND.ARBL.HA.PC"
## [5] "Year.2010.AG.LND.ARBL.ZS"    "Year.2010.AG.LND.CREL.HA"   
## [7] "Year.2010.AG.LND.CROP.ZS"    "Year.2010.AG.LND.FRST.K2"   
## [9] "Year.2010.AG.LND.FRST.ZS"    "Year.2010.AG.LND.TOTL.K2"   
##[11] "Year.2010.AG.PRD.CROP.XD"    "Year.2010.AG.PRD.FOOD.XD"   
##[13] "Year.2010.AG.PRD.LVSK.XD"    "Year.2010.AG.SRF.TOTL.K2"   
##[15] "Year.2010.AG.YLD.CREL.KG"    "Year.2010.EA.PRD.AGRI.KD"   
##[17] "Year.2010.NV.AGR.TOTL.ZS"    "Year.2010.SH.H2O.SAFE.RU.ZS"
##[19] "Year.2010.SP.RUR.TOTL"       "Year.2010.SP.RUR.TOTL.ZG"   
##[21] "Year.2010.SP.RUR.TOTL.ZS"    "Year.2010.TM.VAL.AGRI.ZS.UN"
#[23] "Year.2010.TX.VAL.AGRI.ZS.UN"

# first use Year.2010.AG.PRD.FOOD.XD as response, remove CROP and LVSK production
widesel3<-widesel2[,c(1:10,12,14:23)]
#dim(widesel3)
##[1] 90 21
# define y and x
y0<-widesel3[,11]
x0<-data.frame(widesel3[,-11])

# burn-in? last step? last step cannot be greater than N-1
bn<-21
# number of runs/idependent data sets
nruns<-10
# power of error
powerr<-1
# number of explanatory variables
nvar<-20

# n is multiples of 10, i.e., seq(20,100,10)
nval<-seq(30,dim(widesel3)[1],10)

## for saving results

CPE<-list(flm=matrix(0,nruns,length(nval)),blrc=matrix(0,nruns,length(nval)),blrni=matrix(0,nruns,length(nval)),blrss=matrix(0,nruns,length(nval)),blrg=matrix(0,nruns,length(nval)),fqr=matrix(0,nruns,length(nval)),bqr=matrix(0,nruns,length(nval)))
# BLR parameters
nIter<-5500
burnIn<-500
priorN<-list(varE=list(df=3,S=1))
priorC<-list(varE=list(df=3,S=1),varBR=list(df=3,S=1))
# spikeslab parameters
n.iter1SS<-500
n.iter2SS<-5000	
# g prior
G<-"EBL"
# frequentist quantile regression
Tau<-0.5
# bayesian quantile regression
#	Nbase<-4
Beta0<-rep(0,nvar+1)
V0P<-100*diag(nvar+1)
Shape0<-0.01
Scale0<-0.01
# define global storage
err<-vector("list",length(nval))
predv<-err
# if want to save predictive intervals
#	predvh<-err
#	predvl<-err
fval<-err
# create storage
for(i in 1:length(nval)){
	lst<-nval[i]
# storage for fitted values, predicted values, predictive errors, CPE for each method
	err[[i]]<-list(flm=matrix(0,nruns,lst-bn),blrc=matrix(0,nruns,lst-bn),blrni=matrix(0,nruns,lst-bn),blrss=matrix(0,nruns,lst-bn),blrg=matrix(0,nruns,lst-bn),fqr=matrix(0,nruns,lst-bn),bqr=matrix(0,nruns,lst-bn))
	predv[[i]]<-list(flm=matrix(99,nruns,lst-bn),blrc=matrix(0,nruns,lst-bn),blrni=matrix(0,nruns,lst-bn),blrss=matrix(0,nruns,lst-bn),blrg=matrix(0,nruns,lst-bn),fqr=matrix(0,nruns,lst-bn),bqr=matrix(0,nruns,lst-bn))
	fval[[i]]<-list(flm=matrix(0,nruns,lst-1),blrc=matrix(0,nruns,lst-1),blrni=matrix(0,nruns,lst-1),blrss=matrix(0,nruns,lst-1),blrg=matrix(0,nruns,lst-1),fqr=matrix(0,nruns,lst-1),bqr=matrix(0,nruns,lst-1))
}
# start outer loop	
t1<-proc.time()
for(j in 1:nruns){
# for each run ...
# reorder data
	ord<-sample(dim(widesel3)[1])
	y<-y0[ord]
	x<-x0[ord,]
	for(i in 1:length(nval)){
		lst2<-nval[i]
## frequentist linear models
# for each run do CPE ...
		kk<-bn
		for(kk in bn:(lst2-1)){
			flm<-lm(y~Year.2010.AG.LND.AGRI.K2+Year.2010.AG.LND.AGRI.ZS+Year.2010.AG.LND.ARBL.HA+Year.2010.AG.LND.ARBL.HA.PC+Year.2010.AG.LND.ARBL.ZS+Year.2010.AG.LND.CREL.HA+Year.2010.AG.LND.CROP.ZS+Year.2010.AG.LND.FRST.K2+Year.2010.AG.LND.FRST.ZS+Year.2010.AG.LND.TOTL.K2+Year.2010.AG.SRF.TOTL.K2+Year.2010.AG.YLD.CREL.KG+Year.2010.EA.PRD.AGRI.KD+Year.2010.NV.AGR.TOTL.ZS+Year.2010.SH.H2O.SAFE.RU.ZS+Year.2010.SP.RUR.TOTL+Year.2010.SP.RUR.TOTL.ZG+Year.2010.SP.RUR.TOTL.ZS+Year.2010.TM.VAL.AGRI.ZS.UN+Year.2010.TX.VAL.AGRI.ZS.UN, data=x, subset=c(1:kk))
			predv[[i]]$flm[j,kk-bn+1]<-predict(flm,x[kk+1,],interval="prediction",type="response")[,1]
#			predvl$flm[kk-bn+1]<-predict(flm,x[kk+1,],interval="prediction",type="response")[,2]
#			predvh$flm[kk-bn+1]<-predict(flm,x[kk+1,],interval="prediction",type="response")[,3]			
			err[[i]]$flm[j,kk-bn+1]<-(abs(predv[[i]]$flm[j,kk-bn+1]-y[kk+1]))^powerr
		}
		fval[[i]]$flm[j,]<-flm$fitted
		CPE$flm[j,i]<-sum(err[[i]]$flm[j,])

## bayesian linear models with conjugate priors (R BLR package)
# for each run do CPE ...
		x2<-as.matrix(x)
		kk<-bn
		for(kk in bn:(lst2-1)){
			yNa<-y
			whichNa<-(kk+1):lst2
			yNa[whichNa]<-NA
			# suppress output to screen
			capture.output(blrc<-BLR(y=yNa,XR=x2,prior=priorC,nIter=nIter,burnIn=burnIn,saveAt="example_"),file=NULL)
			predv[[i]]$blrc[j,kk-bn+1]<-blrc$yHat[blrc$whichNa[1]]
			err[[i]]$blrc[j,kk-bn+1]<-(abs(predv[[i]]$blrc[j,kk-bn+1]-y[kk+1]))^powerr			
		}
		fval[[i]]$blrc[j,]<-blrc$yHat[1:kk]
		CPE$blrc[j,i]<-sum(err[[i]]$blrc[j,])

## bayesian linear models with noninformative priors (R BLR package)
# for each run do CPE ...
		kk<-bn
		for(kk in bn:(lst2-1)){
			yNa<-y
			whichNa<-(kk+1):lst2
			yNa[whichNa]<-NA
			# suppress output to screen
			capture.output(blrni<-BLR(y=yNa,XF=x,prior=priorN,nIter=nIter,burnIn=burnIn,saveAt="example_"))
			predv[[i]]$blrni[j,kk-bn+1]<-blrni$yHat[blrni$whichNa[1]]
			err[[i]]$blrni[j,kk-bn+1]<-(abs(predv[[i]]$blrni[j,kk-bn+1]-y[kk+1]))^powerr			
		}
		fval[[i]]$blrni[j,]<-blrni$yHat[1:kk]
		CPE$blrni[j,i]<-sum(err[[i]]$blrni[j,])

## bayesian linear models with spike-and-slab priors (R spikeslab package)
# for each run do CPE ...
		kk<-bn
		x3<-data.frame(y,x$Year.2010.AG.LND.AGRI.K2,x$Year.2010.AG.LND.AGRI.ZS,x$Year.2010.AG.LND.ARBL.HA,x$Year.2010.AG.LND.ARBL.HA.PC,x$Year.2010.AG.LND.ARBL.ZS,x$Year.2010.AG.LND.CREL.HA,x$Year.2010.AG.LND.CROP.ZS,x$Year.2010.AG.LND.FRST.K2,x$Year.2010.AG.LND.FRST.ZS,x$Year.2010.AG.LND.TOTL.K2,x$Year.2010.AG.SRF.TOTL.K2,x$Year.2010.AG.YLD.CREL.KG,x$Year.2010.EA.PRD.AGRI.KD,x$Year.2010.NV.AGR.TOTL.ZS,x$Year.2010.SH.H2O.SAFE.RU.ZS,x$Year.2010.SP.RUR.TOTL,x$Year.2010.SP.RUR.TOTL.ZG,x$Year.2010.SP.RUR.TOTL.ZS,x$Year.2010.TM.VAL.AGRI.ZS.UN,x$Year.2010.TX.VAL.AGRI.ZS.UN)
		for(kk in bn:(lst2-1)){
			yNa<-y
			whichNa<-(kk+1):lst2
			blrss<-spikeslab(y~., x3[-whichNa,], n.iter1=n.iter1SS, n.iter2=n.iter2SS)
			predv[[i]]$blrss[j,kk-bn+1]<-predict(blrss, x3[whichNa,])$yhat.bma[1]
			err[[i]]$blrss[j,kk-bn+1]<-(abs(predv[[i]]$blrss[j,kk-bn+1]-y[kk+1]))^powerr			
		}
		fval[[i]]$blrss[j,]<-predict(blrss, x3)$yhat.bma[1:kk]
		CPE$blrss[j,i]<-sum(err[[i]]$blrss[j,])
#
## bayesian linear models with g-priors (R BMS package)
# for each run do CPE ...
		kk<-bn
		for(kk in bn:(lst2-1)){
			blrg<-zlm(y~Year.2010.AG.LND.AGRI.K2+Year.2010.AG.LND.AGRI.ZS+Year.2010.AG.LND.ARBL.HA+Year.2010.AG.LND.ARBL.HA.PC+Year.2010.AG.LND.ARBL.ZS+Year.2010.AG.LND.CREL.HA+Year.2010.AG.LND.CROP.ZS+Year.2010.AG.LND.FRST.K2+Year.2010.AG.LND.FRST.ZS+Year.2010.AG.LND.TOTL.K2+Year.2010.AG.SRF.TOTL.K2+Year.2010.AG.YLD.CREL.KG+Year.2010.EA.PRD.AGRI.KD+Year.2010.NV.AGR.TOTL.ZS+Year.2010.SH.H2O.SAFE.RU.ZS+Year.2010.SP.RUR.TOTL+Year.2010.SP.RUR.TOTL.ZG+Year.2010.SP.RUR.TOTL.ZS+Year.2010.TM.VAL.AGRI.ZS.UN+Year.2010.TX.VAL.AGRI.ZS.UN, data=x, subset=c(1:kk), g=G)
			predv[[i]]$blrg[j,kk-bn+1]<-predict(blrg, x[kk+1,])
			err[[i]]$blrg[j,kk-bn+1]<-(abs(predv[[i]]$blrg[j,kk-bn+1]-y[kk+1]))^powerr			
		}
		fval[[i]]$blrg[j,]<-blrg$fitted.values
		CPE$blrg[j,i]<-sum(err[[i]]$blrg[j,])

# frequentist quantile regression (median) (R quantreg)
# for each run do CPE ...
		kk<-bn
		for(kk in bn:(lst2-1)){
			fqr<-rq(y~Year.2010.AG.LND.AGRI.K2+Year.2010.AG.LND.AGRI.ZS+Year.2010.AG.LND.ARBL.HA+Year.2010.AG.LND.ARBL.HA.PC+Year.2010.AG.LND.ARBL.ZS+Year.2010.AG.LND.CREL.HA+Year.2010.AG.LND.CROP.ZS+Year.2010.AG.LND.FRST.K2+Year.2010.AG.LND.FRST.ZS+Year.2010.AG.LND.TOTL.K2+Year.2010.AG.SRF.TOTL.K2+Year.2010.AG.YLD.CREL.KG+Year.2010.EA.PRD.AGRI.KD+Year.2010.NV.AGR.TOTL.ZS+Year.2010.SH.H2O.SAFE.RU.ZS+Year.2010.SP.RUR.TOTL+Year.2010.SP.RUR.TOTL.ZG+Year.2010.SP.RUR.TOTL.ZS+Year.2010.TM.VAL.AGRI.ZS.UN+Year.2010.TX.VAL.AGRI.ZS.UN, tau=Tau, data=x, subset=c(1:kk))
			predv[[i]]$fqr[j,kk-bn+1]<-predict(fqr, x[kk+1,])
			err[[i]]$fqr[j,kk-bn+1]<-(abs(predv[[i]]$fqr[j,kk-bn+1]-y[kk+1]))^powerr			
		}
		fval[[i]]$fqr[j,]<-predict(fqr, x[1:kk,])
		CPE$fqr[j,i]<-sum(err[[i]]$fqr[j,])

# bayesian quantile regression (median) (R bayesQR and/or BSquare)
# for each run do CPE ...
#		kk<-bn
#		for(kk in bn:(lst2-1)){
#			y2<-y[1:kk]
#			xsub<-x[1:kk,]
#			bQRprior<-prior(y2~Year.2010.AG.LND.AGRI.K2+Year.2010.AG.LND.AGRI.ZS+Year.2010.AG.LND.ARBL.HA+Year.2010.AG.LND.ARBL.HA.PC+Year.2010.AG.LND.ARBL.ZS+Year.2010.AG.LND.CREL.HA+Year.2010.AG.LND.CROP.ZS+Year.2010.AG.LND.FRST.K2+Year.2010.AG.LND.FRST.ZS+Year.2010.AG.LND.TOTL.K2+Year.2010.AG.SRF.TOTL.K2+Year.2010.AG.YLD.CREL.KG+Year.2010.EA.PRD.AGRI.KD+Year.2010.NV.AGR.TOTL.ZS+Year.2010.SH.H2O.SAFE.RU.ZS+Year.2010.SP.RUR.TOTL+Year.2010.SP.RUR.TOTL.ZG+Year.2010.SP.RUR.TOTL.ZS+Year.2010.TM.VAL.AGRI.ZS.UN+Year.2010.TX.VAL.AGRI.ZS.UN, data=xsub, beta0=Beta0, V0=V0P, shape0=Shape0, scale0=Scale0)
#			# suppress output to screen
#			capture.output(bqr<-bayesQR(y2~Year.2010.AG.LND.AGRI.K2+Year.2010.AG.LND.AGRI.ZS+Year.2010.AG.LND.ARBL.HA+Year.2010.AG.LND.ARBL.HA.PC+Year.2010.AG.LND.ARBL.ZS+Year.2010.AG.LND.CREL.HA+Year.2010.AG.LND.CROP.ZS+Year.2010.AG.LND.FRST.K2+Year.2010.AG.LND.FRST.ZS+Year.2010.AG.LND.TOTL.K2+Year.2010.AG.SRF.TOTL.K2+Year.2010.AG.YLD.CREL.KG+Year.2010.EA.PRD.AGRI.KD+Year.2010.NV.AGR.TOTL.ZS+Year.2010.SH.H2O.SAFE.RU.ZS+Year.2010.SP.RUR.TOTL+Year.2010.SP.RUR.TOTL.ZG+Year.2010.SP.RUR.TOTL.ZS+Year.2010.TM.VAL.AGRI.ZS.UN+Year.2010.TX.VAL.AGRI.ZS.UN, data=xsub, quantile=Tau, ndraw=nIter, prior=bQRprior))
#			sum<-summary(bqr,burnin=n.iter1SS)
#			predv[[i]]$bqr[j,kk-bn+1]<-c(1,as.matrix(x[kk+1,]))%*%sum[[1]]$betadraw[,1]
#			err[[i]]$bqr[j,kk-bn+1]<-(abs(predv[[i]]$bqr[j,kk-bn+1]-y[kk+1]))^powerr			
#		}
#		fval[[i]]$bqr[j,]<-as.vector(as.matrix(cbind(1,xsub))%*%sum[[1]]$betadraw[,1])
#		CPE$bqr[j,i]<-sum(err[[i]]$bqr[j,])
	}
}
# for timing
t2<-proc.time()
t2-t1

# save workspace
save.image(file="ch4reg_WorldBankFood.RData")
# repeat the process above for several different random seeds, save each to a different
# R workspace. Then load all results for plotting.

# load and compile results across saved workspaces
# key components are predv, err, CPE
.libPaths("/Library/Frameworks/R.framework/Versions/Current/library")
library(BMS)
load('ch4reg_WorldBank_1.RData')
dim(predv[[1]]$flm)
#[1] 10  9
predv1<-predv
err1<-err
fval1<-fval
CPE1<-CPE
load('ch4reg_WorldBank_2.RData')
predv2<-predv
err2<-err
fval2<-fval
CPE2<-CPE
load('ch4reg_WorldBank_3.RData')
predv3<-predv
err3<-err
fval3<-fval
CPE3<-CPE
load('ch4reg_WorldBank_4.RData')
predv4<-predv
err4<-err
fval4<-fval
CPE4<-CPE
load('ch4reg_WorldBank_5.RData')
predv5<-predv
err5<-err
fval5<-fval
CPE5<-CPE
load('ch4reg_WorldBank_6.RData')
predv6<-predv
err6<-err
fval6<-fval
CPE6<-CPE
load('ch4reg_WorldBank_7.RData')
predv7<-predv
err7<-err
fval7<-fval
CPE7<-CPE
load('ch4reg_WorldBank_8.RData')
predv8<-predv
err8<-err
fval8<-fval
CPE8<-CPE
load('ch4reg_WorldBank_9.RData')
predv9<-predv
err9<-err
fval9<-fval
CPE9<-CPE
load('ch4reg_WorldBank_Food.RData')
predv10<-predv
err10<-err
fval10<-fval
CPE10<-CPE
#
nrtot<-100
nr<-10
nval<-seq(30,90,10)
bn<-21
err0<-vector("list",length(nval))
predv0<-err0
fval0<-err0
CPE0<-list(flm=matrix(0,nrtot,length(nval)),blrc=matrix(0,nrtot,length(nval)),blrni=matrix(0,nrtot,length(nval)),blrss=matrix(0,nrtot,length(nval)),blrg=matrix(0,nrtot,length(nval)),fqr=matrix(0,nrtot,length(nval)),bqr=matrix(0,nrtot,length(nval)))
for(i in 1:length(nval)){
	lst<-nval[i]
	err0[[i]]<-list(flm=matrix(0,nrtot,lst-bn),blrc=matrix(0,nrtot,lst-bn),blrni=matrix(0,nrtot,lst-bn),blrss=matrix(0,nrtot,lst-bn),blrg=matrix(0,nrtot,lst-bn),fqr=matrix(0,nrtot,lst-bn),bqr=matrix(0,nrtot,lst-bn))
	predv0[[i]]<-list(flm=matrix(99,nrtot,lst-bn),blrc=matrix(0,nrtot,lst-bn),blrni=matrix(0,nrtot,lst-bn),blrss=matrix(0,nrtot,lst-bn),blrg=matrix(0,nrtot,lst-bn),fqr=matrix(0,nrtot,lst-bn),bqr=matrix(0,nrtot,lst-bn))
	fval0[[i]]<-list(flm=matrix(0,nrtot,lst-1),blrc=matrix(0,nrtot,lst-1),blrni=matrix(0,nrtot,lst-1),blrss=matrix(0,nrtot,lst-1),blrg=matrix(0,nrtot,lst-1),fqr=matrix(0,nrtot,lst-1),bqr=matrix(0,nrtot,lst-1))
}
# compile all results for all methods except bqr from separate run files
for(k in 1:10){
	CPE0$flm[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('CPE',k,'$flm[1:',nr,',]',sep="")))
	CPE0$blrc[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('CPE',k,'$blrc[1:',nr,',]',sep="")))
	CPE0$blrni[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('CPE',k,'$blrni[1:',nr,',]',sep="")))
	CPE0$blrss[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('CPE',k,'$blrss[1:',nr,',]',sep="")))
	CPE0$blrg[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('CPE',k,'$blrg[1:',nr,',]',sep="")))
	CPE0$fqr[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('CPE',k,'$fqr[1:',nr,',]',sep="")))
}
for(k in 1:10){
	for(i in 1:length(nval)){
		err0[[i]]$flm[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('err',k,'[[',i,']]$flm[1:',nr,',]',sep="")))
		err0[[i]]$blrc[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('err',k,'[[',i,']]$blrc[1:',nr,',]',sep="")))
		err0[[i]]$blrni[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('err',k,'[[',i,']]$blrni[1:',nr,',]',sep="")))
		err0[[i]]$blrss[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('err',k,'[[',i,']]$blrss[1:',nr,',]',sep="")))
		err0[[i]]$blrg[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('err',k,'[[',i,']]$blrg[1:',nr,',]',sep="")))
		err0[[i]]$fqr[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('err',k,'[[',i,']]$fqr[1:',nr,',]',sep="")))
		predv0[[i]]$flm[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('predv',k,'[[',i,']]$flm[1:',nr,',]',sep="")))
		predv0[[i]]$blrc[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('predv',k,'[[',i,']]$blrc[1:',nr,',]',sep="")))
		predv0[[i]]$blrni[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('predv',k,'[[',i,']]$blrni[1:',nr,',]',sep="")))
		predv0[[i]]$blrss[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('predv',k,'[[',i,']]$blrss[1:',nr,',]',sep="")))
		predv0[[i]]$blrg[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('predv',k,'[[',i,']]$blrg[1:',nr,',]',sep="")))
		predv0[[i]]$fqr[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('predv',k,'[[',i,']]$fqr[1:',nr,',]',sep="")))
		fval0[[i]]$flm[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('fval',k,'[[',i,']]$flm[1:',nr,',]',sep="")))
		fval0[[i]]$blrc[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('fval',k,'[[',i,']]$blrc[1:',nr,',]',sep="")))
		fval0[[i]]$blrni[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('fval',k,'[[',i,']]$blrni[1:',nr,',]',sep="")))
		fval0[[i]]$blrss[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('fval',k,'[[',i,']]$blrss[1:',nr,',]',sep="")))
		fval0[[i]]$blrg[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('fval',k,'[[',i,']]$blrg[1:',nr,',]',sep="")))
		fval0[[i]]$fqr[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('fval',k,'[[',i,']]$fqr[1:',nr,',]',sep="")))
	}
}	
# compile all results for bqr 
load('ch4reg_WorldBankFood_bqr_1054.RData')
dim(predv[[1]]$flm)
#[1] 10  9
predvb1<-predv
errb1<-err
fvalb1<-fval
CPEb1<-CPE
load('ch4reg_WorldBankFood_bqr_1080.RData')
predvb2<-predv
errb2<-err
fvalb2<-fval
CPEb2<-CPE
load('ch4reg_WorldBankFood_bqr_5054.RData')
predvb3<-predv
errb3<-err
fvalb3<-fval
CPEb3<-CPE
load('ch4reg_WorldBankFood_bqr_6054.RData')
predvb4<-predv
errb4<-err
fvalb4<-fval
CPEb4<-CPE
load('ch4reg_WorldBankFood_bqr_7054.RData')
predvb5<-predv
errb5<-err
fvalb5<-fval
CPEb5<-CPE
load('ch4reg_WorldBankFood_bqr_8054.RData')
predvb6<-predv
errb6<-err
fvalb6<-fval
CPEb6<-CPE
load('ch4reg_WorldBankFood_bqr_9054.RData')
predvb7<-predv
errb7<-err
fvalb7<-fval
CPEb7<-CPE
load('ch4reg_WorldBankFood_bqr_10054.RData')
predvb8<-predv
errb8<-err
fvalb8<-fval
CPEb8<-CPE
load('ch4reg_WorldBankFood_bqr_11054.RData')
predvb9<-predv
errb9<-err
fvalb9<-fval
CPEb9<-CPE
load('ch4reg_WorldBankFood_bqr_12054.RData')
predvb10<-predv
errb10<-err
fvalb10<-fval
CPEb10<-CPE
#
# compile all results for bqr from separate run files
for(k in 1:10){
	CPE0$bqr[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('CPEb',k,'$bqr[1:',nr,',]',sep="")))
}
for(k in 1:10){
	for(i in 1:length(nval)){
		err0[[i]]$bqr[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('errb',k,'[[',i,']]$bqr[1:',nr,',]',sep="")))
		predv0[[i]]$bqr[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('predvb',k,'[[',i,']]$bqr[1:',nr,',]',sep="")))
		fval0[[i]]$bqr[((nr*(k-1))+1):(nr*k),]<-eval(parse(text=paste('fvalb',k,'[[',i,']]$bqr[1:',nr,',]',sep="")))
	}
}	
# save workspace as ch4reg_real.RData
#
.libPaths("/Library/Frameworks/R.framework/Versions/Current/library")
load('ch4reg_real.RData')
CPEres<-rbind(apply(CPE0$flm,2,mean),apply(CPE0$blrc,2,mean),apply(CPE0$blrni,2,mean),apply(CPE0$blrss,2,mean),apply(CPE0$blrg,2,mean),apply(CPE0$fqr,2,mean),apply(CPE0$bqr,2,mean))
CPEresVar<-rbind(apply(CPE0$flm,2,var),apply(CPE0$blrc,2,var),apply(CPE0$blrni,2,var),apply(CPE0$blrss,2,var),apply(CPE0$blrg,2,var),apply(CPE0$fqr,2,var),apply(CPE0$bqr,2,var))
CPEresSE<-sqrt(t(t(CPEresVar)/nval))

plot(seq(30,90,10),CPEres[1,],type='b',pch="1",cex.lab=1.2, cex.axis=1.2,ylim=c(90,2300),xlab="Sample size", ylab="CPE", main="CPE by model type")
lines(seq(30,90,10),CPEres[2,],type='b',pch="2")
lines(seq(30,90,10),CPEres[3,],type='b',pch="3")
lines(seq(30,90,10),CPEres[4,],type='b',pch="4")
lines(seq(30,90,10),CPEres[5,],type='b',pch="5")
lines(seq(30,90,10),CPEres[6,],type='b',pch="6")
lines(seq(30,90,10),CPEres[7,],type='b',pch="7")
legend('topleft',legend=c('flm','blrc','blrni','blrss','blrg','fqr','bqr'),pch=c("1","2","3","4","5","6","7"),lty=1, cex=1.3)
# ch4reg_WBFood_CPE_bw_v2.pdf

CPEres1<-CPEres
for(i in 1:dim(CPEres)[1]){ CPEres1[i,]<-CPEres[i,]-CPEres[i,1]}
plot(seq(30,90,10),CPEres1[1,],type='b', pch="1", cex.lab=1.2, cex.axis=1.2,ylim=c(0,1500),xlab="Sample size", ylab="CPE", main="CPE by model type (baseline adjusted)")
lines(seq(30,90,10),CPEres1[2,],type='b',pch="2")
lines(seq(30,90,10),CPEres1[3,],type='b',pch="3")
lines(seq(30,90,10),CPEres1[4,],type='b',pch="4")
lines(seq(30,90,10),CPEres1[5,],type='b',pch="5")
lines(seq(30,90,10),CPEres1[6,],type='b',pch="6")
lines(seq(30,90,10),CPEres1[7,],type='b',pch="7")
legend('topleft',legend=c('flm','blrc','blrni','blrss','blrg','fqr','bqr'),pch=c("1","2","3","4","5","6","7"),lty=1, cex=1.3)
# save as ch4reg_WBFood_CPEadj_bw_v2.pdf

plot(seq(30,90,10),CPEres[1,],type='n',cex.lab=1.2, cex.axis=1.2, ylim=c(90,2800),xlab="Sample size", ylab="CPE", main="CPE by model type")
polygon(c(rev(nval),nval),c(rev(CPEres[1,]+CPEresSE[1,]),CPEres[1,]-CPEresSE[1,]),col='grey80',border=NA)
lines(seq(30,90,10),CPEres[1,])
points(seq(30,90,10),CPEres[1,], pch="1")
lines(nval,CPEres[1,]+CPEresSE[1,],lty='dashed',col='red')
lines(nval,CPEres[1,]-CPEresSE[1,],lty='dashed',col='red')
polygon(c(rev(nval),nval),c(rev(CPEres[2,]+CPEresSE[2,]),CPEres[2,]-CPEresSE[2,]),col='grey20',border=NA)
lines(seq(30,90,10),CPEres[2,],col='grey20')
points(seq(30,90,10),CPEres[2,],pch="2", col='grey20')
lines(nval,CPEres[2,]+CPEresSE[2,],lty='dashed',col='red')
lines(nval,CPEres[2,]-CPEresSE[2,],lty='dashed',col='red')
polygon(c(rev(nval),nval),c(rev(CPEres[3,]+CPEresSE[3,]),CPEres[3,]-CPEresSE[3,]),col='grey40',border=NA)
lines(seq(30,90,10),CPEres[3,],col='grey40')
points(seq(30,90,10),CPEres[3,],pch="3", col='grey40')
lines(nval,CPEres[3,]+CPEresSE[3,],lty='dashed',col='red')
lines(nval,CPEres[3,]-CPEresSE[3,],lty='dashed',col='red')
polygon(c(rev(nval),nval),c(rev(CPEres[4,]+CPEresSE[4,]),CPEres[4,]-CPEresSE[4,]),col='grey60',border=NA)
lines(seq(30,90,10),CPEres[4,],col='grey60')
points(seq(30,90,10),CPEres[4,],pch="4",col='grey60')
lines(nval,CPEres[4,]+CPEresSE[4,],lty='dashed',col='red')
lines(nval,CPEres[4,]-CPEresSE[4,],lty='dashed',col='red')
polygon(c(rev(nval),nval),c(rev(CPEres[5,]+CPEresSE[5,]),CPEres[5,]-CPEresSE[5,]),col='grey10',border=NA)
lines(seq(30,90,10),CPEres[5,],col='grey10')
points(seq(30,90,10),CPEres[5,],pch="5",col='grey10')
lines(nval,CPEres[5,]+CPEresSE[5,],lty='dashed',col='red')
lines(nval,CPEres[5,]-CPEresSE[5,],lty='dashed',col='red')
# grey90 is (229,229,229)
polygon(c(rev(nval),nval),c(rev(CPEres[6,]+CPEresSE[6,]),CPEres[6,]-CPEresSE[6,]),col=rgb(229/255,229/255,229/255,0.8),border=NA)
lines(seq(30,90,10),CPEres[6,])
points(seq(30,90,10),CPEres[6,],pch="6")
lines(nval,CPEres[6,]-CPEresSE[6,],lty='dashed',col='red')
lines(nval,CPEres[6,]+CPEresSE[6,],lty='dashed',col='red')
polygon(c(rev(nval),nval),c(rev(CPEres[7,]+CPEresSE[7,]),CPEres[7,]-CPEresSE[7,]),col='grey95',border=NA)
lines(seq(30,90,10),CPEres[7,])
points(seq(30,90,10),CPEres[7,],pch="7")
lines(nval,CPEres[7,]+CPEresSE[7,],lty='dashed',col='red')
lines(nval,CPEres[7,]-CPEresSE[7,],lty='dashed',col='red')
polygon(c(rev(nval),nval),c(rev(CPEres[1,]+CPEresSE[1,]),CPEres[1,]-CPEresSE[1,]),col='grey80',border=NA)
lines(seq(30,90,10),CPEres[1,])
points(seq(30,90,10),CPEres[1,], pch="1")
lines(nval,CPEres[1,]+CPEresSE[1,],lty='dashed',col='red')
lines(nval,CPEres[1,]-CPEresSE[1,],lty='dashed',col='red')
polygon(c(rev(nval),nval),c(rev(CPEres[6,]+CPEresSE[6,]),CPEres[6,]-CPEresSE[6,]),col=rgb(229/255,229/255,229/255,0.8),border=NA)
lines(seq(30,90,10),CPEres[6,])
points(seq(30,90,10),CPEres[6,],pch="6")
lines(nval,CPEres[6,]-CPEresSE[6,],lty='dashed',col='red')
lines(nval,CPEres[6,]+CPEresSE[6,],lty='dashed',col='red')
legend('topleft',legend=c('flm','blrc','blrni','blrss','blrg','fqr','bqr'),pch=c("1","2","3","4","5","6","7"),lty=1, cex=1.3)
# save as ch4reg_WBFood_CPEsd_v2.pdf
