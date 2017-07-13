# 	ch10.6 computing to examine stability

#	two clusters and three years
#	SD of response is ~90
#	add N(0,sigma) error to the outcome and examine impact on predictive ability
#	see computing_vegout_plots.txt for previous results shown in book

library(rpart)
library(rpart.plot)
library(caret)
library(caretEnsemble)
library(neuralnet)

load('veg.RData')
load('vegout_clustering.RData')
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

## number of clusters and years to predict
nclust<-2
wclust<-c(1,3)
nynew<-8
## levels of noise to add assuming response has SD=90
sd<-90
sdlev<-c(0,0.2,0.6,0.8)
## number of repetitions
nsim<-100

##	treed and doppler, with trees and NNs, no ensemble plus all six ensemble methods
## storage will be redone for each method, so the below only applies to one
# storage for model information
# for one cluster and several levels of added noise

tmod<-vector("list",nynew-nyb)
nnmod<-vector("list",nynew-nyb)
for(j in 1:(nynew-nyb)){
	tmod[[j]]<-vector("list",length(sdlev))
	nnmod[[j]]<-vector("list",length(sdlev))
}	
tpred<-vector("list",length=nynew-nyb)
for(j in 1:(nynew-nyb)) tpred[[j]]<-vector("list",length=length(sdlev))
nnpred<-vector("list",length=nynew-nyb)
for(j in 1:(nynew-nyb)) nnpred[[j]]<-vector("list",length=length(sdlev))
 
set.seed(9051)

ptm <- proc.time()
## with ensembling - boosting trees	- cluster 1 only
j<-1
for(j in 1:nsim){
	assign(paste("tmod", j, sep=""),tmod)
	assign(paste("tpred", j, sep=""),tpred)
	for(k in 1:(nynew-nyb)){
		nyb2<-nyb+k-1
		## select by cluster, then predict next year
		cmem<-c(rep(cluster[[k]],nyb2+1))
	# fit model within cluster, use to predict
		c1<-datn2[1:(nloc*(nyb2+1)),c(1:43,46,47)]
		c1<-c1[which(cmem==wclust[1]),]
		c1<-as.data.frame(c1)
#		# neuralnet cannot handle any missing data
 		t<-apply(c1,1,function(x){sum(is.na(x))})
		c1<-c1[which(t<1),]
		c1<-c1[which(!is.na(c1$SSG0_6)),]
		tmp<-apply(c1,2,var,use="pairwise.complete.obs")
		c1<-c1[,which(tmp>0)]
		## add random error
		for(kk in 1:length(sdlev)){
			c1$SSG0_6<-c1$SSG0_6+rnorm(dim(c1)[1],mean=0,sd=sd*sdlev[kk])
		#dim(c1)
		#[1] 1200  47
			## with ensembling: boosting - trees
			trControl<-trainControl(method="none")
			tree5t<-train(SSG0_6~., data=c1[c1$Year<as.numeric(Year[nyb2+1]),], method="blackboost", metric="RMSE", tuneLength=1, trControl=trControl, tuneGrid=expand.grid(mstop=50, maxdepth=5))
			assign(paste("tmod",j,"[[",k,"]][[",kk,"]]",sep=""),tree5t)
			assign(paste("tpred",j,"[[",k,"]][[",kk,"]]",sep=""),cbind(c1[c1$Year==as.numeric(Year[nyb2+1]),44],predict(tree5t,c1[c1$Year==as.numeric(Year[nyb2+1]),])))
			cat(kk, "of ", length(sdlev), "\r") 
			flush.console()
			}
		cat(k, "of ", (nynew-nyb), "\r") 
		flush.console()
	}
	cat(j, "of ", nsim, "\r") 
	flush.console()
}
save.image('vegout_boosting_trees_sim_c1.RData')

# plot results - run 5
par(mfrow=c(nynew-nyb,length(sdlev)))
j<-1
for(j in 1:(nynew-nyb)){
		for(m in 1:length(sdlev)){
			plot(1:dim(tpred5[[j]][[m]])[1],tpred5[[j]][[m]][,1],type="l",main=paste("Year",Year[j+nyb]," and SD ",sdlev[m],sep=" "),xlab="",ylab="", ylim=c(-350,350), cex.main=0.8)
			lines(1:dim(tpred5[[j]][[m]])[1],tpred5[[j]][[m]][,2],col="red")
		}
}
# 
library(beeswarm)
# plot results - run 1
par(mfrow=c(nynew-nyb,1),oma=c(0,0,0,0), mar=c(3,2,2,1))
for(i in 1:(nynew-nyb)){
tmp<-NULL
	for(m in 1:length(sdlev)){
		tmp2<-rep(sdlev[m],length(tpred[[i]][[m]][,1]))
		tmp2<-cbind(tmp2,(tpred[[i]][[m]][,1]-tpred[[i]][[m]][,2]))
		tmp<-rbind(tmp,tmp2)
	}
	tmp<-as.data.frame(tmp)
	colnames(tmp)<-c('sdlev','error')
	beeswarm(tmp$error~tmp$sdlev, method=c("swarm"),pch=16,xlab='sdlev',ylab='error',main=paste("Year ",Year[nyb+i],sep=""),ylim=c(-300,400), cex=0.5) 
}

# make a beeswarm of average errors over runs
library(beeswarm)
par(mfrow=c(nynew-nyb,1),oma=c(0,0,0,0), mar=c(3,2,2,1))
for(i in 1:(nynew-nyb)){
	tmp<-NULL
	for(j in 1:100){
		eval(parse(text=paste("tpred0<-tpred",j,sep="")))
		for(m in 1:length(sdlev)){
			tmp2<-sdlev[m]
			tmp2<-cbind(tmp2,mean(tpred0[[i]][[m]][,1]-tpred0[[i]][[m]][,2]))
			tmp<-rbind(tmp,tmp2)
		}
	}
	tmp<-as.data.frame(tmp)
	colnames(tmp)<-c('sdlev','error')
#	boxplot(tmp$error~tmp$sdlev, method=c("swarm"),pch=16,xlab='sdlev',ylab='error',main=paste("Year ",Year[nyb+i],sep=""),ylim=c(-40,30), cex=0.5) 
	beeswarm(tmp$error~tmp$sdlev, method=c("swarm"),pch=16,xlab='sdlev',ylab='error',main=paste("Year ",Year[nyb+i],sep=""),ylim=c(-50,50), cex=0.5) 
}

ptm <- proc.time()
## with ensembling - median	neural net - cluster 1 only
set.seed(9051)

mygrid<-expand.grid(.decay=c(0.5, 0.1, 5e-4), .size=c(2,4,6))
mygrid3<-expand.grid(.decay=0.1, .size=c(1:5))

j<-1
for(j in 1:nsim){
	nnmod0<-nnmod
	nnpred0<-nnpred
	assign(paste("nnmod", j, sep=""),nnmod)
	assign(paste("nnpred", j, sep=""),nnpred)
	for(k in 1:(nynew-nyb)){
		nyb2<-nyb+k-1
	## select by cluster, then predict next year
		cmem<-c(rep(cluster[[k]],nyb2+1))
	# fit model within cluster, use to predict
		c1<-datn2[1:(nloc*(nyb2+1)),c(1:43,46,47)]
		c1<-c1[which(cmem==wclust[1]),]
		c1<-as.data.frame(c1)
		# neuralnet cannot handle any missing data
 		t<-apply(c1,1,function(x){sum(is.na(x))})
		c1<-c1[which(t<1),]
		c1<-c1[which(!is.na(c1$SSG0_6)),]
		tmp<-apply(c1,2,var,use="pairwise.complete.obs")
		c1<-c1[,which(tmp>0)]
		## add random error
		for(kk in 1:length(sdlev)){
			c1$SSG0_6<-c1$SSG0_6+rnorm(dim(c1)[1],mean=0,sd=sd*sdlev[kk])
			## with ensembling median - neural net
			nobs<-dim(c1[c1$Year<as.numeric(Year[nyb2+1]),])[1]
			nobs2<-dim(c1[c1$Year==as.numeric(Year[nyb2+1]),])[1]	
				c0<-c1[c1$Year<as.numeric(Year[nyb2+1]),]
				# scale data
				maxs <- apply(c0, 2, max, na.rm=T) 
				mins <- apply(c0, 2, min, na.rm=T)
				sc0 <- as.data.frame(scale(c0, center = mins, scale = maxs - mins))
				n <- names(sc0)
				f <- as.formula(paste("SSG0_6 ~ ", paste(n[!n %in% "SSG0_6"], collapse = " + ")))
			post<-pred<-NULL
			for(jj in 1:dim(mygrid3)[1]){
				nn2t <- neuralnet(f, data=sc0, hidden=mygrid3[jj,2], threshold=0.02, linear.output=T, stepmax=800000)
 					c2<-c1[c1$Year==as.numeric(Year[nyb2+1]),]
					sc2 <- as.data.frame(scale(c2, center = mins, scale = maxs - mins))
					pr.nn <- compute(nn2t,sc2[,-(dim(sc2)[2]-1)])
					pr.nn2 <- pr.nn$net.result*(max(c0$SSG0_6, na.rm=T)-min(c0$SSG0_6, na.rm=T))+min(c0$SSG0_6, na.rm=T)
					test.r <- (sc2$SSG0_6)*(max(c0$SSG0_6, na.rm=T)-min(c0$SSG0_6, na.rm=T))+min(c0$SSG0_6, na.rm=T)
					pred<-c(pred,pr.nn2)
					post<-c(post,rep(exp(-0.5*(sqrt(sum(abs(test.r-pr.nn2))))+log(nobs)*log(2*pi)+(1+mygrid3[j,2])*log(nobs)),nobs2))
			}
			nnmod0[[k]][[kk]]<-cbind(post,pred)
			tmp<-matrix(post,nobs2,dim(mygrid3)[1])
			tmp1<-t(apply(tmp,1,function(x){x/sum(x)}))
			tmp2<-matrix(pred,nobs2,dim(mygrid3)[1])
			## select models
			jj<-1
			tmp4<-NULL
			for(jj in 1:nobs2){
				tmp3 = myweightmedian2(tmp1[1,],tmp2[jj,])
				wmBICnewy = tmp3[1]
				## compute median prediction
				tmp4<-c(tmp4,wmBICnewy)
			}
			nnpred0[[k]][[kk]]<-cbind(c1[c1$Year==as.numeric(Year[nyb2+1]),44],tmp4)
		}
		cat(k, "of ", (nynew-nyb), "\r\r") 
		flush.console()
	}
	assign(paste("nnmod",j,sep=""),nnmod0)
	assign(paste("nnpred",j,sep=""),nnpred0)
	cat(j, "of ", nsim, "\r") 
	flush.console()
}
save.image('vegout_med_nnet_sim_1-50.RData')

# plot results - run 5
par(mfrow=c(nynew-nyb,length(sdlev)))
j<-1
for(j in 1:(nynew-nyb)){
		for(m in 1:length(sdlev)){
			plot(1:dim(nnpred5[[j]][[m]])[1],nnpred5[[j]][[m]][,1],type="l",main=paste("Year",Year[j+nyb]," and SD ",sdlev[m],sep=" "),xlab="",ylab="", ylim=c(-350,350), cex.main=0.8)
			lines(1:dim(nnpred5[[j]][[m]])[1],nnpred5[[j]][[m]][,2],col="red")
		}
}
# 
library(beeswarm)
# plot results - run 1
par(mfrow=c(nynew-nyb,1),oma=c(0,0,0,0), mar=c(3,2,2,1))
for(i in 1:(nynew-nyb)){
tmp<-NULL
	for(m in 1:length(sdlev)){
		tmp2<-rep(sdlev[m],length(nnpred5[[i]][[m]][,1]))
		tmp2<-cbind(tmp2,(nnpred5[[i]][[m]][,1]-nnpred5[[i]][[m]][,2]))
		tmp<-rbind(tmp,tmp2)
	}
	tmp<-as.data.frame(tmp)
	colnames(tmp)<-c('sdlev','error')
	beeswarm(tmp$error~tmp$sdlev, method=c("swarm"),pch=16,xlab='sdlev',ylab='error',main=paste("Year ",Year[nyb+i],sep=""),ylim=c(-300,400), cex=0.5) 
}
# make a beeswarm of average errors over runs
library(beeswarm)
par(mfrow=c(nynew-nyb,1),oma=c(0,0,0,0), mar=c(3,2,2,1))
for(i in 1:(nynew-nyb)){
	tmp<-NULL
	for(j in 1:50){
		eval(parse(text=paste("nnpred0<-nnpred",j,sep="")))
		for(m in 1:length(sdlev)){
			tmp2<-sdlev[m]
			tmp2<-cbind(tmp2,mean(nnpred0[[i]][[m]][,1]-nnpred0[[i]][[m]][,2]))
			tmp<-rbind(tmp,tmp2)
		}
	}
	tmp<-as.data.frame(tmp)
	colnames(tmp)<-c('sdlev','error')
	beeswarm(tmp$error~tmp$sdlev, method=c("swarm"),pch=16,xlab='sdlev',ylab='error',main=paste("Year ",Year[nyb+i],sep=""),ylim=c(-150,100), cex=0.5) 
}

