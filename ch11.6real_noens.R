# ch11.6 computing
# trees and neural nets, no ensembling

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
## without ensembling
## neuralnet really needs stepmax of 500k to ensure convergence, otherwise crash ...
for(k in 1:(ny-nyb)){
	nyb2<-nyb+k-1
	## select by cluster, then predict next year
	cmem<-c(rep(cluster[[k]],nyb2+1))
	# fit model within cluster, use to predict
	for(kk in 1:nclust){
		c1<-datn2[1:(nloc*(nyb2+1)),c(1:43,46,47)]
		c1<-c1[which(cmem==kk),]
		## if doing neural net ...
		# neuralnet cannot handle any missing data
 		t<-apply(c1,1,function(x){sum(is.na(x))})
		c1<-c1[which(t<1),]
		c1<-as.data.frame(c1)
		c1<-c1[which(!is.na(c1$SSG0_6)),]
		if(dim(c1)[1]>0){		
			tmp<-apply(c1,2,var,use="pairwise.complete.obs")
			c1<-c1[,which(tmp>0)]
		#dim(c1)
		#[1] 1200  45
			if(dim(c1)[1]>0){
			## without ensembling - trees
				tree1t<-rpart(SSG0_6~., data=c1[c1$Year<as.numeric(Year[nyb2+1]),], method="anova") 
				opt <- which.min(tree1t$cptable[,"xerror"])
				cp <- tree1t$cptable[opt, "CP"]
				tree1t_prune<-prune(tree1t, cp = cp)
				tmod[[k]][[kk]]<-tree1t_prune
				ptree1<-predict(tree1t_prune, newdata=c1[c1$Year==as.numeric(Year[nyb2+1]),])
				tpred[[k]][[kk]]<-cbind(c1$SSG0_6[c1$Year==as.numeric(Year[nyb2+1])],ptree1)
			## without ensembling - neural net
				nn1t<-nnet(SSG0_6~., data=c1[c1$Year<as.numeric(Year[nyb2+1]),], size=4, decay=5e-4, maxit=1000, linout=T, trace=F)
				nnmod[[k]][[kk]]<-nn1t
				nnpred[[k]][[kk]]<-predict(nn1t, newdata=c1[c1$Year==as.numeric(Year[nyb2+1]),], type = "raw")

				c0<-c1[c1$Year<as.numeric(Year[nyb2+1]),]
				# scale data
				maxs <- apply(c0, 2, max, na.rm=T) 
				mins <- apply(c0, 2, min, na.rm=T)
				sc0 <- as.data.frame(scale(c0, center = mins, scale = maxs - mins))
				n <- names(sc0)
				f <- as.formula(paste("SSG0_6 ~ ", paste(n[!n %in% "SSG0_6"], collapse = " + ")))
				nnet1t <- neuralnet(f,data=sc0, hidden=4,linear.output=T)
 				nnmod[[k]][[kk]]<-nnet1t
 				c2<-c1[c1$Year==as.numeric(Year[nyb2+1]),]
				sc2 <- as.data.frame(scale(c2, center = mins, scale = maxs - mins))
				pr.nn <- compute(nnet1t,sc2[,-(dim(sc2)[2]-1)])
				pr.nn2 <- pr.nn$net.result*(max(c0$SSG0_6, na.rm=T)-min(c0$SSG0_6, na.rm=T))+min(c0$SSG0_6, na.rm=T)
				test.r <- (sc2$SSG0_6)*(max(c0$SSG0_6, na.rm=T)-min(c0$SSG0_6, na.rm=T))+min(c0$SSG0_6, na.rm=T)
				nnpred[[k]][[kk]]<-cbind(test.r,pr.nn2)
			}
		}
	}
	cat(k, "of ", (ny-nyb), "\r") 
	flush.console()
}

# plot results
tmp3<-vector("list",length=(ny-nyb))
for(j in 1:length(tmp3)){
tmp4<-NULL
m<-1
for(m in 1:nclust){if((length(tpred[[j]][[m]])>0)&&(sum(is.na(tpred[[j]][[m]])<1)))	 tmp4<-rbind(tmp4,cbind(tpred[[j]][[m]],rep(m,dim(tpred[[j]][[m]])[1])))}
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
		tmp2<-cbind(tmp2,(nnpred[[i]][[1]][,1]-nnpred[[i]][[1]][,2]))
		tmp<-rbind(tmp,tmp2)
	}
}
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),pch=16,xlab='',ylab='error',cex=0.5) 
