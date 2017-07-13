# ch10.6 computing

#	trees and NNs, stacking

library(rpart)
library(rpart.plot)
library(ipred)
library(nnet)
library(caret)
library(caretEnsemble)
library(e1071)
library(party)
library(mboost)
library(plyr)
library(BayesTree)
library(bartMachine)
library(kpodclustr)
library(neuralnet)

load('veg.RData')
## stack data so variables that don't change by location vary
## variables 2-9 do not vary with time so only enter once and repeat
varsel<-c(1:9)
step<-37
datn2<-datn[,varsel]
for(j in 2:ny) datn2<-rbind(datn2,datn[,varsel])
start<-9
tmp<-datn[,start+seq(1,step,1)]
for(j in 2:ny) tmp<-rbind(tmp,datn[,start+step*(j-1)+seq(1,step,1)])
datn2<-cbind(datn2,tmp)
tmp<-as.numeric(c(rep(Year[1:ny],each=nloc)))
datn2<-cbind(datn2,tmp)
colnames(datn2)<-c(colnames(datn2)[1:start],colnames(dat)[c(start+1+1:step,10)])
#write.csv(datn2, file="datn2.csv")

myweightmedian2=function(weight, mass)
{
  z = cbind(weight, mass)
  u = rep(0, length(weight))
  v=u
  u[rank(z[,2], ties.method="random")]=z[,1]
  v[rank(z[,2], ties.method="random")]=z[,2]
  ## add next 2 lines to get order of models rank sorted
  t = rank(z[,2], ties.method="random")
  tsort = sort(t,index.return=TRUE)
  z = cbind(cumsum(u),v)
  weightm = z[z[,1]>=0.5,2]
  ## add next two lines and adjust return
  mod = tsort$ix[z[,1]>=0.5]

  return(c(weightm[1],mod[1]))
}

## burn in (years)
nyb<-5
## minimum cluster size
minc<-20

## number of bootstraps
## proportion of data in each bootstrap run
nboot<-100
bootfrac<-0.6
## number of clusters
nclust<-6
cluster<-vector("list",length=ny-nyb)

set.seed(9050)

## do clustering of data FIRST, then use in all subsequent methods
for(k in 1:(ny-nyb)){
	## cluster locations using first five years of data as initial burnin
	##note oceanographic variables don't matter
	nyb2<-nyb+k-1
	# clustering by k-means that allows missing data
	# go with 6 clusters
	clust1<-kpod(datn[,1:(9+step*nyb2)],k=5)
	niter<-length(clust1$fit_list)
	cluster[[k]]<-clust1$cluster_list[[niter]]
}

## storage will be redone for each method, so the below only applies to one
# storage for model information

tmod<-vector("list",ny-nyb)
nnmod<-vector("list",ny-nyb)
for(j in 1:(ny-nyb)){
	tmod[[j]]<-vector("list",nclust)
	nnmod[[j]]<-vector("list",nclust)
}	
tpred<-vector("list",length=ny-nyb)
for(j in 1:(ny-nyb)) tpred[[j]]<-vector("list",length=nclust)
nnpred<-vector("list",length=ny-nyb)
for(j in 1:(ny-nyb)) nnpred[[j]]<-vector("list",length=nclust)
 
# set of possible model parameter settings for neural networks single layer
mygrid<-expand.grid(.decay=c(0.5, 0.1, 5e-4), .size=c(2,4,6))
mygrid2<-expand.grid(.decay=c(0.5, 0.1, 5e-4), .size=c(2,4,6), .bag=c(4,6,8))
 
set.seed(9051)

ptm <- proc.time()
## with ensembling: stacking
for(k in 1:(ny-nyb)){
	nyb2<-nyb+k-1
	## select by cluster, then predict next year
	cmem<-c(rep(cluster[[k]],nyb2+1))
	# fit model within cluster, use to predict
	for(kk in 1:nclust){
		c1<-datn2[1:(nloc*(nyb2+1)),c(1:43,46,47)]
		c1<-c1[which(cmem==kk),]
		c1<-as.data.frame(c1)
		# neuralnet cannot handle any missing data
 		t<-apply(c1,1,function(x){sum(is.na(x))})
		c1<-c1[which(t<1),]
		c1<-c1[which(!is.na(c1$SSG0_6)),]
		if(dim(c1)[1]>0){		
			tmp<-apply(c1,2,var,use="pairwise.complete.obs")
			c1<-c1[,which(tmp>0)]
		#dim(c1)
		#[1] 1200  47
			if(dim(c1)[1]>minc){
			## with ensembling: stacking - trees
				nobs<-dim(c1[c1$Year<as.numeric(Year[nyb2+1]),])[1]
				tmp<-NULL
				mod1<-rpart(SSG0_6~., data=c1[c1$Year<as.numeric(Year[nyb2+1]),], method="anova") 
				tmp<-cbind(tmp,predict(mod1,c1[c1$Year<=as.numeric(Year[nyb2+1]),]))
				mod2<-rpart(SSG0_6~., data=c1[c1$Year<as.numeric(Year[nyb2+1]),], method="anova", control=rpart.control(minsplit=10, minbucket=3, cp=0.05))
				if(!identical(predict(mod1,c1[c1$Year<as.numeric(Year[nyb2+1]),]),predict(mod2,c1[c1$Year<as.numeric(Year[nyb2+1]),]))) tmp<-cbind(tmp,predict(mod2,c1[c1$Year<=as.numeric(Year[nyb2+1]),]))
				mod3<-rpart(SSG0_6~., data=c1[c1$Year<as.numeric(Year[nyb2+1]),], method="anova", control=rpart.control(minsplit=5, minbucket=3))
				if(sum(as.numeric(apply(tmp,2,function(x){identical(predict(mod3,c1[c1$Year<as.numeric(Year[nyb2+1]),]),x)})))<nobs) tmp<-cbind(tmp,predict(mod3,c1[c1$Year<=as.numeric(Year[nyb2+1]),]))
				mod4<-ctree(SSG0_6~., data=c1[c1$Year<as.numeric(Year[nyb2+1]),])
				if(sum(as.numeric(apply(tmp,2,function(x){identical(predict(mod4,c1[c1$Year<as.numeric(Year[nyb2+1]),]),x)})))<nobs) tmp<-cbind(tmp,predict(mod4,c1[c1$Year<=as.numeric(Year[nyb2+1]),]))
				mod5<-ctree(SSG0_6~., data=c1[c1$Year<as.numeric(Year[nyb2+1]),], controls=ctree_control(stump=TRUE))
				if(sum(as.numeric(apply(tmp,2,function(x){identical(predict(mod5,c1[c1$Year<as.numeric(Year[nyb2+1]),]),x)})))<nobs) tmp<-cbind(tmp,predict(mod5,c1[c1$Year<=as.numeric(Year[nyb2+1]),]))
			#
				datt2<-as.data.frame(cbind(c1[c1$Year<=as.numeric(Year[nyb2+1]),44],tmp))
				names(datt2)<-c(paste("V",1:dim(datt2)[2],sep=""))
				tree4t<-lm(V1 ~ ., data=datt2,subset=1:nobs) 
				tmod[[k]][[kk]]<-tree4t
				tpred[[k]][[kk]]<-predict(tree4t,datt2[(nobs+1):(dim(datt2)[1]),])
			## with ensembling: stacking - neural nets
				tmp<-NULL
				for(j in 1:6){
#				for(j in 1:dim(mygrid)[1]){
					mod<-train(SSG0_6~., data=c1[c1$Year<as.numeric(Year[nyb2+1]),], "nnet", tuneGrid=mygrid[j,], linout=TRUE, maxit=100, trace=F)
					if(j<2) tmp<-cbind(tmp,predict(mod, c1[c1$Year<=as.numeric(Year[nyb2+1]),]))
					if(j>1) if(sum(as.numeric(apply(tmp,2,function(x){identical(predict(mod,c1[c1$Year<as.numeric(Year[nyb2+1]),]),x)})))<nobs) tmp<-cbind(tmp,predict(mod,c1[c1$Year<=as.numeric(Year[nyb2+1]),]))
				}
				datt2<-as.data.frame(cbind(c1[c1$Year<=as.numeric(Year[nyb2+1]),44],tmp))
				names(datt2)<-c(paste("V",1:dim(datt2)[2],sep=""))
				nn4t<-lm(V1 ~ ., data=datt2, subset=1:nobs) 
				nnmod[[k]][[kk]]<-nn4t
				nnpred[[k]][[kk]]<-predict(nn4t,datt2[(nobs+1):(dim(datt2)[1]),])
			}
		}
	cat(kk, "of ", nclust, "\r") 
	flush.console()
	}
	cat(k, "of ", (ny-nyb), "\r") 
	flush.console()
}
