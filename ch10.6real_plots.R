# ch10.6 computing

#	plot combined results for vegout data across all methods

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

# load data
load('veg.RData')

## load workspaces with results of clustering and all analysis methods
# clustering results
load('vegout_clustering.RData')
# no ensembling
load(`vegout_nonens.RData')
tmod1<-tmod
tpred1<-tpred
load(`vegout_nonens_nnet_all.RData')
nnmod1<-nnmod
nnpred1<-nnpred
# bma
load('vegout_bma_tree.RData')
tmod2<-tmod
tpred2<-tpred
load("vegout_bma_nnet.RData")
nnmod2<-nnmod
nnpred2<-nnpred2
# bagging
load("vegout_bagging_trees.RData")
tmod3<-tmod
tpred3<-tpred
load("vegout_bagging_nnet.RData")
nnmod3<-nnmod
nnpred3<-nnpred
# stacking
load("vegout_stacking.RData")
tmod4<-tmod
tpred4<-tpred
nnmod4<-nnmod
nnpred4<-nnpred
# boosting
load("vegout_boosting_trees.RData")
tmod5<-tmod
tpred5<-tpred
load("vegout_boosting_nnet.RData")
nnmod5<-nnmod
nnpred5<-nnpred
# median method
load("vegout_median_trees.RData")
tmod6<-tmod
tpred6<-tpred
load("vegout_median_nnet.RData")
nnmod6<-nnmod
nnpred6<-nnpred

# beeswarm plots of predictive errors
par(mfrow=c(3,2), oma=c(0,0,0,0), mar=c(3,2,2,1))
library(beeswarm)
tmp<-NULL
for(i in 1:13){
 	if(!is.null(tpred1[[i]][[1]])){
 		tmp2<-rep(i,length(tpred1[[i]][[1]][,1]))
  		tmp2<-cbind(tmp2,(tpred1[[i]][[1]][,1]-tpred1[[i]][[1]][,2]))
 		tmp<-rbind(tmp,tmp2)
 	}
 }
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300), pch=16,xlab='year',ylab='error',cex=0.2, cex.main=0.7, main="tree predictive errors")

tmp<-NULL
for(i in 1:13){
	if(!is.null(nnpred1[[i]][[1]])){
		tmp2<-rep(i,length(nnpred1[[i]][[1]][,1]))
		tmp2<-cbind(tmp2,(nnpred1[[i]][[1]][,1]-nnpred1[[i]][[1]][,2]))
		tmp<-rbind(tmp,tmp2)
	}
	if(is.null(nnpred1[[i]][[1]])){
		tmp2<-rep(i,1)
		tmp2<-cbind(tmp2,0)
		tmp<-rbind(tmp,tmp2)
	}

}
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300),pch=16,xlab='year',ylab='error',cex=0.2, cex.main=0.7,main="neural net predictive errors") 

tmp<-NULL
for(i in 1:13){
 	if(!is.null(tpred2[[i]][[1]])){
 		tmp2<-rep(i,length(tpred2[[i]][[1]][,1]))
  		tmp2<-cbind(tmp2,(tpred2[[i]][[1]][,1]-tpred2[[i]][[1]][,2]))
 		tmp<-rbind(tmp,tmp2)
 	}
 }
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300), pch=16,xlab='year',ylab='error',cex=0.2, cex.main=0.7,main="tree BMA predictive errors")

tmp<-NULL
for(i in 1:13){
	if(!is.null(nnpred2[[i]][[1]])){
		tmp2<-rep(i,length(nnpred2[[i]][[1]]))
		tmp2<-cbind(tmp2,(tpred2[[i]][[1]][,1]-nnpred2[[i]][[1]]))
		tmp<-rbind(tmp,tmp2)
	}
	if(is.null(nnpred2[[i]][[1]])){
		tmp2<-rep(i,1)
		tmp2<-cbind(tmp2,0)
		tmp<-rbind(tmp,tmp2)
	}

}
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300),pch=16,xlab='year',ylab='error',cex=0.2, cex.main=0.7,main="nnet BMA predictive errors") 

tmp<-NULL
for(i in 1:13){
 	if(!is.null(tpred3[[i]][[1]])){
 		tmp2<-rep(i,length(tpred3[[i]][[1]][,1]))
  		tmp2<-cbind(tmp2,(tpred3[[i]][[1]][,1]-tpred3[[i]][[1]][,2]))
 		tmp<-rbind(tmp,tmp2)
 	}
 }
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300), pch=16,xlab='year',ylab='error',cex=0.2, cex.main=0.7,main="tree bag predictive errors")

tmp<-NULL
for(i in 1:13){
	if(!is.null(nnpred3[[i]][[1]])){
		tmp2<-rep(i,length(nnpred3[[i]][[1]]))
		tmp2<-cbind(tmp2,(tpred3[[i]][[1]][,1]-as.vector(nnpred3[[i]][[1]])))
		tmp<-rbind(tmp,tmp2)
	}
	if(is.null(nnpred3[[i]][[1]])){
		tmp2<-rep(i,1)
		tmp2<-cbind(tmp2,0)
		tmp<-rbind(tmp,tmp2)
	}

}
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300),pch=16,xlab='year',ylab='error',cex=0.2, cex.main=0.7,main="nnet bag predictive errors") 

par(mfrow=c(3,2), oma=c(0,0,0,0), mar=c(3,2,2,1))
tmp<-NULL
for(i in 1:13){
 	if(!is.null(tpred4[[i]][[1]])){
 		tmp2<-rep(i,length(tpred3[[i]][[1]][,1]))
  		tmp2<-cbind(tmp2,(tpred3[[i]][[1]][,1]-tpred4[[i]][[1]]))
 		tmp<-rbind(tmp,tmp2)
 	}
 }
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300), pch=16,xlab='year',ylab='error',cex=0.2, cex.main=0.7,main="tree stack predictive errors")

tmp<-NULL
for(i in 1:13){
	if(!is.null(nnpred4[[i]][[1]])){
		tmp2<-rep(i,length(nnpred4[[i]][[1]]))
		tmp2<-cbind(tmp2,(tpred3[[i]][[1]][,1]-nnpred4[[i]][[1]]))
		tmp<-rbind(tmp,tmp2)
	}
	if(is.null(nnpred4[[i]][[1]])){
		tmp2<-rep(i,1)
		tmp2<-cbind(tmp2,0)
		tmp<-rbind(tmp,tmp2)
	}

}
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300),pch=16,xlab='year',ylab='error',cex=0.2, cex.main=0.7,main="nnet stack predictive errors") 

tmp<-NULL
for(i in 1:13){
 	if(!is.null(tpred5[[i]][[1]])){
 		tmp2<-rep(i,length(tpred5[[i]][[1]][,1]))
  		tmp2<-cbind(tmp2,(tpred5[[i]][[1]][,1]-tpred5[[i]][[1]][,2]))
 		tmp<-rbind(tmp,tmp2)
 	}
 }
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300), pch=16,xlab='year',ylab='error',cex=0.2, cex.main=0.7,main="tree boost predictive errors")

tmp<-NULL
for(i in 1:13){
	if(!is.null(nnpred5[[i]][[1]])){
		tmp2<-rep(i,length(nnpred5[[i]][[1]]))
		tmp2<-cbind(tmp2,(tpred5[[i]][[1]][,1]-nnpred5[[i]][[1]]))
		tmp<-rbind(tmp,tmp2)
	}
	if(is.null(nnpred5[[i]][[1]])){
		tmp2<-rep(i,1)
		tmp2<-cbind(tmp2,0)
		tmp<-rbind(tmp,tmp2)
	}

}
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300),pch=16,xlab='year',ylab='error',cex=0.2, cex.main=0.7, main="nnet boost predictive errors") 

tmp<-NULL
for(i in 1:13){
 	if(!is.null(tpred6[[i]][[1]])){
 		tmp2<-rep(i,length(tpred6[[i]][[1]][,1]))
  		tmp2<-cbind(tmp2,(tpred6[[i]][[1]][,1]-tpred6[[i]][[1]][,2]))
 		tmp<-rbind(tmp,tmp2)
 	}
 }
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300), pch=16,xlab='year',ylab='error',cex=0.2, cex.main=0.7,main="tree med predictive errors")

tmp<-NULL
for(i in 1:13){
	if(!is.null(nnpred6[[i]][[1]])){
		tmp2<-rep(i,length(nnpred6[[i]][[1]][,1]))
		tmp2<-cbind(tmp2,(nnpred6[[i]][[1]][,1]-nnpred6[[i]][[1]][,2]))
		tmp<-rbind(tmp,tmp2)
	}
	if(is.null(nnpred6[[i]][[1]])){
		tmp2<-rep(i,1)
		tmp2<-cbind(tmp2,0)
		tmp<-rbind(tmp,tmp2)
	}

}
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300),pch=16,xlab='year',ylab='error',cex=0.2, cex.main=0.7,main="nnet med predictive errors") 

# observed verses predicted plots for data in cluster 1, with responses ordered

par(mfrow=c(3,2), oma=c(0,0,0,0), mar=c(3,2,2,1))
tmp<-sort(tpred1[[1]][[1]][,1],index.return=T)
plot(1:340,tpred1[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="tree predictive errors")
points(1:340,tpred1[[1]][[1]][tmp$ix,2],col="red")
tmp<-sort(nnpred1[[1]][[1]][,1],index.return=T)
plot(1:297,nnpred1[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="neural net predictive errors")
points(1:297,nnpred1[[1]][[1]][tmp$ix,2],col="red")

tmp<-sort(tpred2[[1]][[1]][,1],index.return=T)
plot(1:340,tpred2[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="tree BMA predictive errors")
points(1:340,tpred2[[1]][[1]][tmp$ix,2],col="red")
tmp<-sort(tpred2[[1]][[1]][,1],index.return=T)
plot(1:340,tpred2[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="nnet BMA predictive errors")
points(1:340,nnpred2[[1]][[1]][tmp$ix],col="red")

tmp<-sort(tpred3[[1]][[1]][,1],index.return=T)
plot(1:297,tpred3[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="tree bag predictive errors")
points(1:297,tpred3[[1]][[1]][tmp$ix,2],col="red")
tmp<-sort(tpred3[[1]][[1]][,1],index.return=T)
plot(1:297,tpred3[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="nnet bag predictive errors")
points(1:297,as.vector(nnpred3[[1]][[1]][tmp$ix]),col="red")
# actual_pred_y1c1_1
par(mfrow=c(3,2), oma=c(0,0,0,0), mar=c(3,2,2,1))
tmp<-sort(tpred3[[1]][[1]][,1],index.return=T)
plot(1:297,tpred3[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="tree stack predictive errors")
points(1:297,tpred4[[1]][[1]][tmp$ix],col="red")
tmp<-sort(tpred3[[1]][[1]][,1],index.return=T)
plot(1:297,tpred3[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="neural stack predictive errors")
points(1:297,nnpred4[[1]][[1]][tmp$ix],col="red")

tmp<-sort(tpred5[[1]][[1]][,1],index.return=T)
plot(1:297,tpred5[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="tree boost predictive errors")
points(1:297,tpred5[[1]][[1]][tmp$ix,2],col="red")
tmp<-sort(tpred5[[1]][[1]][,1],index.return=T)
plot(1:297,tpred5[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="nnet boost predictive errors")
points(1:297,nnpred5[[1]][[1]][tmp$ix],col="red")

tmp<-sort(tpred6[[1]][[1]][,1],index.return=T)
plot(1:340,tpred6[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="tree med predictive errors")
points(1:340,tpred6[[1]][[1]][tmp$ix,2],col="red")
tmp<-sort(tpred6[[1]][[1]][,1],index.return=T)
plot(1:340,tpred6[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="nnet med predictive errors")
points(1:340,nnpred6[[1]][[1]][tmp$ix],col="red")
# actual_pred_y1c1_2

par(mfrow=c(3,2), oma=c(0,0,0,0), mar=c(3,2,2,1))
tmp<-sort(tpred1[[12]][[1]][,1],index.return=T)
plot(1:275,tpred1[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="tree predictive errors")
points(1:275,tpred1[[12]][[1]][tmp$ix,2],col="red")
tmp<-sort(tpred1[[12]][[1]][,1],index.return=T)
plot(1:275,tpred1[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="neural net predictive errors")
points(1:275,nnpred1[[12]][[1]][tmp$ix],col="red")

tmp<-sort(tpred2[[12]][[1]][,1],index.return=T)
plot(1:733,tpred2[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="tree BMA predictive errors")
points(1:733,tpred2[[12]][[1]][tmp$ix,2],col="red")
tmp<-sort(tpred2[[12]][[1]][,1],index.return=T)
plot(1:733,tpred2[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="nnet BMA predictive errors")
points(1:733,nnpred2[[12]][[1]][tmp$ix],col="red")

tmp<-sort(tpred3[[12]][[1]][,1],index.return=T)
plot(1:587,tpred3[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="tree bag predictive errors")
points(1:587,tpred3[[12]][[1]][tmp$ix,2],col="red")
tmp<-sort(tpred3[[12]][[1]][,1],index.return=T)
plot(1:587,tpred3[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="nnet bag predictive errors")
points(1:587,as.vector(nnpred3[[12]][[1]][tmp$ix]),col="red")
# actual_pred_y12c1_1
par(mfrow=c(3,2), oma=c(0,0,0,0), mar=c(3,2,2,1))
tmp<-sort(tpred3[[12]][[1]][,1],index.return=T)
plot(1:587,tpred3[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="tree stack predictive errors")
points(1:587,tpred4[[12]][[1]][tmp$ix],col="red")
tmp<-sort(tpred3[[12]][[1]][,1],index.return=T)
plot(1:587,tpred3[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="neural stack predictive errors")
points(1:587,nnpred4[[12]][[1]][tmp$ix],col="red")

tmp<-sort(tpred5[[12]][[1]][,1],index.return=T)
plot(1:587,tpred5[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="tree boost predictive errors")
points(1:587,tpred5[[12]][[1]][tmp$ix,2],col="red")
tmp<-sort(tpred5[[12]][[1]][,1],index.return=T)
plot(1:587,tpred5[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="nnet boost predictive errors")
points(1:587,nnpred5[[12]][[1]][tmp$ix],col="red")

tmp<-sort(tpred6[[12]][[1]][,1],index.return=T)
plot(1:733,tpred6[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="tree med predictive errors")
points(1:733,tpred6[[12]][[1]][tmp$ix,2],col="red")
tmp<-sort(tpred6[[12]][[1]][,1],index.return=T)
plot(1:733,tpred6[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", ylim=c(-200,200), cex=0.2, cex.main=0.7, main="nnet med predictive errors")
points(1:733,nnpred6[[12]][[1]][tmp$ix],col="red")

## boxplots of errors by cluster by year for nnet median (nnpred1) and tree boost (tpred5)
tmp3<-vector("list",length=(ny-nyb))
for(j in 1:length(tmp3)){
	tmp4<-NULL
	m<-1
	for(m in 1:nclust){
		if(!is.null(nnpred6[[j]][[m]])){
			tmp4<-rbind(tmp4,cbind(nnpred6[[j]][[m]],rep(m,length(nnpred6[[j]][[m]][,1]))))
			}
		else{
			tmp4<-rbind(tmp4,cbind(rep(NA,2),rep(NA,2),rep(m,2)))
			}
	}	
	tmp3[[j]]<-tmp4
}
tmp3nnMed<-tmp3
tmp4<-NULL
for(i in 1:(ny-nyb)) tmp4<-rbind(tmp4,cbind(tmp3nnMed[[i]],rep(as.numeric(Year[i+nyb]),dim(tmp3nnMed[[i]])[1])))
boxplot(abs(tmp4[,1]-tmp4[,2])~tmp4[,3]+tmp4[,4],col=gray(1:6/6),varwidth=T, notch=T, cex=0.4)
# no data for 2006
tmp<-boxplot(abs(tmp4[,1]-tmp4[,2])~as.factor(tmp4[,3])+as.factor(tmp4[,4]),at=rep(8*(0:12),each=6)+rep(1:6,13),col=gray(1:6/6),varwidth=T, notch=T, cex=0.4, boxwex=0.7, axes=FALSE)
#tmp<-boxplot(abs(tmp4[,1]-tmp4[,2])~as.factor(tmp4[,3])+as.factor(tmp4[,4]),at=rep(8*(0:11),each=6)+rep(1:6,12),col=gray(1:6/6),varwidth=T, notch=T, cex=0.4, boxwex=0.7, axes=FALSE)
axis(1, at=seq(4,92,by=8), labels=Year[(nyb+1):(ny-1)], cex.axis=0.6)
axis(2)

tmp3<-vector("list",length=(ny-nyb))
for(j in 1:length(tmp3)){
	tmp4<-NULL
	m<-1
	for(m in 1:nclust){
		if(!is.null(tpred5[[j]][[m]])){
			tmp4<-rbind(tmp4,cbind(tpred5[[j]][[m]],rep(m,length(tpred5[[j]][[m]][,1]))))
			}
		else{
			tmp4<-rbind(tmp4,cbind(rep(NA,2),rep(NA,2),rep(m,2)))
			}
	}	
	tmp3[[j]]<-tmp4
}
tmp3tBoo<-tmp3
tmp4<-NULL
for(i in 1:(ny-nyb)) tmp4<-rbind(tmp4,cbind(tmp3tBoo[[i]],rep(as.numeric(Year[i+nyb]),dim(tmp3tBoo[[i]])[1])))
boxplot(abs(tmp4[,1]-tmp4[,2])~tmp4[,3]+tmp4[,4],col=gray(1:6/6),varwidth=T, notch=T, cex=0.4)
# no data for 2006
tmp<-boxplot(abs(tmp4[,1]-tmp4[,2])~as.factor(tmp4[,3])+as.factor(tmp4[,4]),at=rep(8*(0:12),each=6)+rep(1:6,13),col=gray(1:6/6),varwidth=T, notch=T, cex=0.4, boxwex=0.7, axes=FALSE)
#tmp<-boxplot(abs(tmp4[,1]-tmp4[,2])~as.factor(tmp4[,3])+as.factor(tmp4[,4]),at=rep(8*(0:11),each=6)+rep(1:6,12),col=gray(1:6/6),varwidth=T, notch=T, cex=0.4, boxwex=0.7, axes=FALSE)
axis(1, at=seq(4,92,by=8), labels=Year[(nyb+1):(ny-1)], cex.axis=0.6)
axis(2)

