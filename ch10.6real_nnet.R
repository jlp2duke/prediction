#	VegOut data and computing for Chapter 10.6 - neural net

# For background information about the data see ch10.6real_tree.R

# prediction methods: 
#	ind'l trees
#	NN/single hidden layer+two and three hidden layers

dat<-read.csv('vegout1345test.csv', header=T)
dat2<-as.matrix(dat)
Year<-as.character(seq(1989,2006,1))
ny<-length(Year)
nloc<-dim(dat)[1]/length(Year)
step<-37
datn<-matrix(999,nloc,9+(step*ny))
colnames(datn)<-colnames(datn, do.NULL = FALSE, prefix = "col")
colnames(datn)[1:9]<-names(dat)[1:9]
colnames(datn)[10:dim(datn)[2]]<-as.vector(outer(names(dat)[11:dim(dat)[2]], Year, paste, sep="."))
datn[,1:9]<-dat2[seq(1,(ny*(nloc-1)+1),ny),1:9]

for(i in 1:nloc){
		datn[i,10:dim(datn)[2]]<-as.vector(t(dat2[(i-1)*ny+1:ny,11:dim(dat2)[2]]))
}
#save.image('veg.RData')

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
write.csv(datn2, file="datn2.csv")
	
library(kpodclustr)
library(neuralnet)
# see http://www.r-bloggers.com/fitting-a-neural-network-in-r-neuralnet-package/
# neural net with two hidden layers of sizes 5 and 3 nodes
## make constructs to save results
k<-1
j<-1
nyb<-5
nclust<-6
cluster<-vector("list",length=ny-nyb)
nnets<-vector("list",length=ny-nyb)
for(j in 1:(ny-nyb)) nnets[[j]]<-vector("list",length=nclust)
nnetpred<-vector("list",length=ny-nyb-1)
for(j in 1:(ny-nyb)) nnetpred[[j]]<-vector("list",length=nclust)
for(k in 1:(ny-nyb)){
	## cluster locations using first five years of data as initial burnin
	##note oceanographic variables don't matter
	set.seed(10945)
	nyb2<-nyb+k-1
	# clustering by k-means that allows missing data
	# go with 6 clusters
	clust1<-kpod(datn[,1:(9+step*nyb2)],k=6)
	niter<-length(clust1$fit_list)
	cluster[[k]]<-clust1$cluster_list[[niter]]
	## select by cluster, then predict next year
	cmem<-c(rep(cluster[[k]],nyb2+1))
	# fit model within cluster, use to predict
	# use CV to find architecture if n>30 in cluster
	hidden<-matrix(c(3,3,4,4,2,3,2,3),4,2)
	for(kk in 1:nclust){
		c1<-datn2[1:(nloc*(nyb2+1)),c(1:43,46,47)]
		c1<-c1[which(cmem==kk),]
		c1<-as.data.frame(c1)
		# neuralnet cannot handle any missing data
 		t<-apply(c1,1,function(x){sum(is.na(x))})
		c1<-c1[which(t<1),]
		if(dim(c1)[1]>0){		
			tmp<-apply(c1,2,var,use="pairwise.complete.obs")
			c1<-c1[,which(tmp>0)]
		#dim(c1)
		#[1] 1200  47
			if((dim(c1)[1]>0)&(dim(c1)[1]<30)){
				c0<-c1[c1$Year<as.numeric(Year[nyb2+1]),]
				# scale data
				maxs <- apply(c0, 2, max) 
				mins <- apply(c0, 2, min)
				sc0 <- as.data.frame(scale(c0, center = mins, scale = maxs - mins))
				n <- names(sc0)
				f <- as.formula(paste("SSG0_6 ~", paste(n[!n %in% "SSG0_6"], collapse = " + ")))
				nnet1 <- neuralnet(f,data=sc0, hidden=hidden[1,],linear.output=T)
 				nnets[[k]][[kk]]<-nnet1
 				c2<-c1[c1$Year==as.numeric(Year[nyb2+1]),]
				sc2 <- as.data.frame(scale(c2, center = mins, scale = maxs - mins))
				pr.nn <- compute(nnet1,sc2[,-(dim(sc2)[2]-1)])
				pr.nn2 <- pr.nn$net.result*(max(c0$SSG0_6)-min(c0$SSG0_6))+min(c0$SSG0_6)
				test.r <- (sc2$SSG0_6)*(max(c0$SSG0_6)-min(c0$SSG0_6))+min(c0$SSG0_6)
				nnetpred[[k]][[kk]]<-cbind(test.r,pr.nn2)
			}
			if(dim(c1)[1]>30){
				c0<-c1[c1$Year<as.numeric(Year[nyb2+1]),]
				# scale data
				maxs <- apply(c0, 2, max) 
				mins <- apply(c0, 2, min)
				sc0 <- as.data.frame(scale(c0, center = mins, scale = maxs - mins))
				n <- names(sc0)
				f <- as.formula(paste("SSG0_6 ~", paste(n[!n %in% "SSG0_6"], collapse = " + ")))
				# run cv loop
				kfold <- 10
				nmod<-nrow(hidden)
				cv.error<-rep(99999999,nmod)
				library(plyr) 
				for(ii in 1:nmod){
					pbar <- create_progress_bar('text')
					pbar$init(kfold)
					for(i in 1:kfold){
    					index <- sample(1:nrow(c0),round(0.9*nrow(c0)))
    					train.cv <- sc0[index,]
    					test.cv <- sc0[-index,]
    					nn <- neuralnet(f,data=train.cv,hidden=hidden[ii,],linear.output=T)
    					if(!is.null(nn$net.result)){
    						pr.nn <- compute(nn,test.cv[,-(dim(sc0)[2]-1)])
    						pr.nn <- pr.nn$net.result*(max(c0$SSG0_6)-min(c0$SSG0_6))+min(c0$SSG0_6)
							test.cv.r <- (test.cv$SSG0_6)*(max(c0$SSG0_6)-min(c0$SSG0_6))+min(c0$SSG0_6)
    						cv.error[ii] <- cv.error[ii]+sum((test.cv.r - pr.nn)^2)/nrow(test.cv)
    						pbar$step()
    					}
					}
				}
				# select best architecture and run
				bhidden<-hidden[which(cv.error==min(cv.error))[1],]
				nnet1 <- neuralnet(f,data=sc0, hidden=bhidden,linear.output=T)
 				nnets[[k]][[kk]]<-nnet1
 				c2<-c1[c1$Year==as.numeric(Year[nyb2+1]),]
				sc2 <- as.data.frame(scale(c2, center = mins, scale = maxs - mins))
				if(!is.null(nnet1$net.result)){
					pr.nn <- compute(nnet1,sc2[,-(dim(sc2)[2]-1)])
					pr.nn2 <- pr.nn$net.result*(max(c0$SSG0_6)-min(c0$SSG0_6))+min(c0$SSG0_6)
					test.r <- (sc2$SSG0_6)*(max(c0$SSG0_6)-min(c0$SSG0_6))+min(c0$SSG0_6)
					nnetpred[[k]][[kk]]<-cbind(test.r,pr.nn2)
				}
			}
		}
	}
}

# plot results
tmp3<-vector("list",length=(ny-nyb))
j<-1
for(j in 1:length(tmp3)){
	tmp4<-NULL
	m<-1
	for(m in 1:nclust){if((length(nnetpred[[j]][[m]])>0)&&(sum(is.na(nnetpred[[j]][[m]])<1))) tmp4<-rbind(tmp4,cbind(nnetpred[[j]][[m]],rep(m,dim(nnetpred[[j]][[m]])[1])))}
	if(length(tmp4)>0) tmp3[[j]]<-tmp4
}
par(mfrow=c(4,3))
j<-1
for(j in 1:length(tmp3)){
	if(length(tmp3[[j]]>0)){
		plot(1:dim(tmp3[[j]])[1],tmp3[[j]][,1],type="l",main=paste("Year",Year[j+nyb],sep=" ",xlab="",ylab=""),xlab="",ylab="")
		lines(1:dim(tmp3[[j]])[1],tmp3[[j]][,2],col="red")
	}
}

tmp<-sort(nnetpred[[1]][[1]][,1],index.return=T)
plot(1:dim(nnetpred[[1]][[1]])[1],nnetpred[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)")
points(1:dim(nnetpred[[1]][[1]])[1],nnetpred[[1]][[1]][tmp$ix,2],col="red")

# beeswarm plot
library(beeswarm)
tmp<-NULL
for(i in 1:13){
	if(!is.null(nnetpred[[i]][[1]])){
		tmp2<-rep(i,length(nnetpred[[i]][[1]][,1]))
		tmp2<-cbind(tmp2,(nnetpred[[i]][[1]][,1]-nnetpred[[i]][[1]][,2]))
		tmp<-rbind(tmp,tmp2)
	}
	if(is.null(nnetpred[[i]][[1]])){
		tmp2<-rep(i,1)
		tmp2<-cbind(tmp2,0)
		tmp<-rbind(tmp,tmp2)
	}

}
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),pch=16,xlab='',ylab='error',cex=0.3) 
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),add=T,pch=2,xlab='',ylab='error',cex=0.3)
legend('topright',legend=c('tree','neuralnet'),pch=c(16,2))

nnetsize<-vector('list',length=2)
for(i in 1:2) nnetsize[[i]]<-matrix(NA,nclust,(ny-nyb))
for(jj in 1:nclust){
	for(j in 1:(ny-nyb)){
			if(!is.null(nnets[[j]][[jj]]$weights)){ 
				nnetsize[[1]][jj,j]<-dim(nnets[[j]][[jj]]$weights[[1]][[1]])[2]
				nnetsize[[2]][jj,j]<-dim(nnets[[j]][[jj]]$weights[[1]][[2]])[2]
		}
	}
}
cell<-matrix(0,dim(nnetsize[[1]])[1],dim(nnetsize[[1]])[2])
for(jj in 1:nclust){
	for(j in 1:(ny-nyb)){
		if(!is.null(nnets[[j]][[jj]]$weights)){
			 cell[jj,j]<-paste(nnetsize[[1]][jj,j],nnetsize[[2]][jj,j],sep=",")
		}
		if(is.null(nnets[[j]][[jj]]$weights)){
		 	cell[jj,j]<-"NA"
		}
	}
}
heatmap.2(nnetsize[[1]]+nnetsize[[2]],Rowv=FALSE, Colv=FALSE, dendrogram='none',cellnote=cell,col=gray.colors(30),srtCol=360,adjCol=c(0,1),notecol="white",tracecol="white")

tmp4<-NULL
for(i in 1:(ny-nyb)) tmp4<-rbind(tmp4,cbind(tmp3nnet[[i]],rep(as.numeric(Year[i+nyb]),dim(tmp3nnet[[i]])[1])))
boxplot(abs(tmp4[,1]-tmp4[,2])~tmp4[,3]+tmp4[,4],col=gray(1:6/6),varwidth=T, notch=T, cex=0.4)
# no data for year 2006
tmp<-boxplot(abs(tmp4[,1]-tmp4[,2])~as.factor(tmp4[,3])+as.factor(tmp4[,4]),at=rep(8*(0:11),each=6)+rep(1:6,12),col=gray(1:6/6),varwidth=T, notch=T, cex=0.4, boxwex=0.7, axes=FALSE)
axis(1, at=seq(4,92,by=8), labels=Year[(nyb+1):(ny-1)], cex.axis=0.6)
axis(2)
