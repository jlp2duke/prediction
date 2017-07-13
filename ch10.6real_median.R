# ch10.6 computing

#	trees and NNs, median model

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
## with ensembling - median	tree
for(k in 1:(ny-nyb)){
	nyb2<-nyb+k-1
	## select by cluster, then predict next year
	cmem<-c(rep(cluster[[k]],nyb2+1))
	# fit model within cluster, use to predict
	for(kk in 1:nclust){
		c1<-datn2[1:(nloc*(nyb2+1)),c(1:43,46,47)]
		c1<-c1[which(cmem==kk),]
		c1<-as.data.frame(c1)
		c1<-c1[which(!is.na(c1$SSG0_6)),]
		if(dim(c1)[1]>0){		
			tmp<-apply(c1,2,var,use="pairwise.complete.obs")
			c1<-c1[,which(tmp>0)]
		#dim(c1)
		#[1] 1200  47
			if(dim(c1)[1]>minc){
			## with ensembling median - tree
			## approx BIC with N/sig*[mse + logN*(d/N)*sig]
			## use e^{-0.5BIC}/sum[e^{-0.5BIC}] as weights
				nobs<-dim(c1[c1$Year<as.numeric(Year[nyb2+1]),])[1]
				nobs2<-dim(c1[c1$Year==as.numeric(Year[nyb2+1]),])[1]
				rp.ctrl<-expand.grid(minsplit=c(5,10,15), cp=c(0.01, 0.05), xval=c(3,5,10), maxdepth=c(2,6))
				post<-pred<-NULL
				for(j in 1:dim(rp.ctrl)[1]){
					tree6t<-rpart(SSG0_6~., data=c1[c1$Year<as.numeric(Year[nyb2+1]),], method="anova", control=rpart.control(minsplit=rp.ctrl[j,1], cp=rp.ctrl[j,2], xval=rp.ctrl[j,3], maxdepth=rp.ctrl[j,4]))
					mse<-sum(tree6t$frame[(tree6t$frame$var=="<leaf>"),]$dev/tree6t$frame[(tree6t$frame$var=="<leaf>"),]$wt)
					post<-c(post,rep(exp(-0.5*(sqrt(sqrt(nobs*mse))+log(nobs)*(dim(tree6t$frame)[1]/nobs))),nobs2))
#					post<-c(post,rep(exp(-0.5*((nobs*mse)+log(nobs)*(dim(tree6t$frame)[1]/nobs))),nobs2))
					pred<-c(pred,predict(tree6t, newdata=c1[c1$Year==as.numeric(Year[nyb2+1]),]))
				}
				tmod[[k]][[kk]]<-cbind(post,pred)
				tmp<-matrix(post,nobs2,dim(rp.ctrl)[1])
				tmp1<-t(apply(tmp,1,function(x){x/sum(x)}))
				tmp2<-matrix(pred,nobs2,dim(rp.ctrl)[1])
				## select models
				j<-1
				tmp4<-NULL
				for(j in 1:nobs2){
					tmp3 = myweightmedian2(tmp1[1,],tmp2[j,])
					wmBICnewy = tmp3[1]
					## compute median prediction
					tmp4<-c(tmp4,wmBICnewy)
				}
				tpred[[k]][[kk]]<-cbind(c1[c1$Year==as.numeric(Year[nyb2+1]),44],tmp4)
			}
		}
	cat(kk, "of ", nclust, "\r") 
	flush.console()
	}
	cat(k, "of ", (ny-nyb), "\r") 
	flush.console()
}

# plot results
tmp3<-vector("list",length=(ny-nyb))
for(j in 1:length(tmp3)){
tmp4<-NULL
m<-1
for(m in 1:nclust){if((dim(tpred[[j]][[m]])>0)&&(sum(is.na(tpred[[j]][[m]][,1])<1)))	 tmp4<-rbind(tmp4,cbind(tpred[[j]][[m]],rep(m,length(tpred[[j]][[m]][,1]))))}
tmp3[[j]]<-tmp4
}
par(mfrow=c(4,3))
j<-2
for(j in 2:13){
plot(1:dim(tmp3[[j]])[1],tmp3[[j]][,1],type="l",main=paste("Year",Year[j+nyb],sep=" ",xlab="",ylab=""),xlab="",ylab="")
lines(1:dim(tmp3[[j]])[1],tmp3[[j]][,2],col="red")
}
# 
tmp<-sort(tpred[[1]][[1]][,1],index.return=T)
plot(1:340,tpred[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)")
points(1:340,tpred[[1]][[1]][tmp$ix,2],col="red")

library(beeswarm)
tmp<-NULL
for(i in 1:13){
	if(!is.null(tpred[[i]][[1]])){
		tmp2<-rep(i,length(tpred[[i]][[1]][,1]))
		tmp2<-cbind(tmp2,(tpred[[i]][[1]][,2]-tpred[[i]][[1]][,1]))
		tmp<-rbind(tmp,tmp2)
	}
}
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),pch=16,xlab='',ylab='error',cex=0.5) 


ptm <- proc.time()
## with ensembling - median	neural net
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
			## with ensembling median - neural net
				nobs<-dim(c1[c1$Year<as.numeric(Year[nyb2+1]),])[1]
				nobs2<-dim(c1[c1$Year==as.numeric(Year[nyb2+1]),])[1]
				post<-pred<-NULL
				for(j in 1:dim(mygrid)[1]){
					fitControl<-trainControl(method="boot", number=6)
					nn6t<-train(SSG0_6~., data=c1[c1$Year<as.numeric(Year[nyb2+1]),], "nnet", tuneGrid=mygrid[j,], trControl=fitControl, linout=TRUE, maxit=1000, trace=F)
					post<-c(post,rep(exp(-0.5*(sqrt(unlist(nn6t$results[3]))+log(nobs)*log(2*pi)+(1+mygrid[j,2])*log(nobs))),nobs2))
#					post<-c(post,rep(exp(-0.5*(unlist(nn6t$results[3])+nobs*log(2*pi)+(1+mygrid[j,2])*log(nobs))),nobs2))
					pred<-c(pred,predict(nn6t$finalModel, newdata=c1[c1$Year==as.numeric(Year[nyb2+1]),]))
				}
				nnmod[[k]][[kk]]<-cbind(post,pred)
				tmp<-matrix(post,nobs2,dim(mygrid)[1])
				tmp1<-t(apply(tmp,1,function(x){x/sum(x)}))
				tmp2<-matrix(pred,nobs2,dim(mygrid)[1])
				## select models
				j<-1
				tmp4<-NULL
				for(j in 1:nobs2){
					tmp3 = myweightmedian2(tmp1[1,],tmp2[j,])
					wmBICnewy = tmp3[1]
					## compute median prediction
					tmp4<-c(tmp4,wmBICnewy)
				}
				nnpred[[k]][[kk]]<-cbind(c1[c1$Year==as.numeric(Year[nyb2+1]),44],tmp4)
				}
			}
		}
	cat(kk, "of ", nclust, "\r") 
	flush.console()
	}
	cat(k, "of ", (ny-nyb), "\r") 
	flush.console()
}

# plot results
tmp3<-vector("list",length=(ny-nyb))
for(j in 1:length(tmp3)){
tmp4<-NULL
m<-1
for(m in 1:nclust){if((dim(nnpred[[j]][[m]])>0)&&(sum(is.na(nnpred[[j]][[m]][,1])<1)))	 tmp4<-rbind(tmp4,cbind(nnpred[[j]][[m]],rep(m,length(nnpred[[j]][[m]][,1]))))}
tmp3[[j]]<-tmp4
}
par(mfrow=c(4,3))
j<-2
for(j in 2:13){
plot(1:dim(tmp3[[j]])[1],tmp3[[j]][,1],type="l",main=paste("Year",Year[j+nyb],sep=" ",xlab="",ylab=""),xlab="",ylab="")
lines(1:dim(tmp3[[j]])[1],tmp3[[j]][,2],col="red")
}
# 
tmp<-sort(nnpred[[1]][[1]][,1],index.return=T)
plot(1:297,nnpred[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)")
points(1:297,nnpred[[1]][[1]][tmp$ix,2],col="red")

library(beeswarm)
tmp<-NULL
for(i in 1:13){
	if(!is.null(nnpred[[i]][[1]])){
		tmp2<-rep(i,length(nnpred[[i]][[1]][,1]))
		tmp2<-cbind(tmp2,(nnpred[[i]][[1]][,2]-nnpred[[i]][[1]][,1]))
		tmp<-rbind(tmp,tmp2)
	}
}
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),pch=16,xlab='',ylab='error',cex=0.5) 
