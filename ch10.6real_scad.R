#	VegOut data and computing for Chapter 10.6 - scad

# For information about the data see ch10.6real_tree.R

# prediction methods: 
#	linear regression with penalized/SCAD+aLASSO

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
library(ncvreg)
# see https://cran.r-project.org/web/packages/ncvreg/index.html
## make constructs to save results
k<-1
j<-1
nyb<-5
nclust<-6
cluster<-vector("list",length=ny-nyb)
ncvregs<-vector("list",length=ny-nyb)
for(j in 1:(ny-nyb)) ncvregs[[j]]<-vector("list",length=nclust)
ncvregpred<-vector("list",length=ny-nyb-1)
for(j in 1:(ny-nyb)) ncvregpred[[j]]<-vector("list",length=nclust)
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
	for(kk in 1:nclust){
		c1<-datn2[1:(nloc*(nyb2+1)),c(1:43,46,47)]
		c1<-c1[which(cmem==kk),]
		c1<-as.data.frame(c1)
		# ncvreg cannot handle any missing data
 		t<-apply(c1,1,function(x){sum(is.na(x))})
		c1<-c1[which(t<1),]
		if(dim(c1)[1]>0){		
			tmp<-apply(c1,2,var,use="pairwise.complete.obs")
			c1<-c1[,which(tmp>0)]
		#dim(c1)
		#[1] 1200  47
			if(dim(c1)[1]>0){
				c0<-as.matrix(c1[c1$Year<as.numeric(Year[nyb2+1]),])
 				ncvreg1<-ncvreg(X=c0[,-(dim(c1)[2]-1)], y=c0[,(dim(c1)[2]-1)], family="gaussian", penalty="SCAD")
 				ncvregs[[k]][[kk]]<-ncvreg1
 				lambda<-ncvreg1$lambda[which(ncvreg1$loss==min(ncvreg1$loss))[1]]
 				c2<-as.matrix(c1[c1$Year==as.numeric(Year[nyb2+1]),])
				prvm1<-predict(ncvreg1, X=c2[,-(dim(c1)[2]-1)], type="response", lambda=lambda)
				ncvregpred[[k]][[kk]]<-cbind(c2[,(dim(c1)[2]-1)],prvm1)
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
	for(m in 1:nclust){if((length(ncvregpred[[j]][[m]])>0)&&(sum(is.na(ncvregpred[[j]][[m]])<1))) tmp4<-rbind(tmp4,cbind(ncvregpred[[j]][[m]],rep(m,dim(ncvregpred[[j]][[m]])[1])))}
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

tmp<-sort(ncvregpred[[1]][[1]][,1],index.return=T)
plot(1:dim(ncvregpred[[1]][[1]])[1],ncvregpred[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", main="SCAD predictions")
points(1:dim(ncvregpred[[1]][[1]])[1],ncvregpred[[1]][[1]][tmp$ix,2],col="red")

tmp<-NULL
for(i in 1:13){
	if(!is.null(ncvregpred[[i]][[1]])){
		tmp2<-rep(i,length(ncvregpred[[i]][[1]][,1]))
		tmp2<-cbind(tmp2,(ncvregpred[[i]][[1]][,1]-ncvregpred[[i]][[1]][,2]))
		tmp<-rbind(tmp,tmp2)
	}
	if(is.null(svmpred[[i]][[1]])){
		tmp2<-rep(i,1)
		tmp2<-cbind(tmp2,0)
		tmp<-rbind(tmp,tmp2)
	}

}
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),add=T,pch=16,xlab='',ylab='error',cex=0.3)

legend('topright',legend=c('rvm','svm','SCAD'),pch=c(3,4,6),col=c('black','red','green'))
tmp<-NULL
for(i in 1:13){
 	if(!is.null(ncvregpred[[i]][[1]])){
 		tmp2<-rep(i,length(ncvregpred[[i]][[1]][,1]))
 		tmp2<-cbind(tmp2,(ncvregpred[[i]][[1]][,1]-ncvregpred[[i]][[1]][,2]))
 		tmp<-rbind(tmp,tmp2)
 	}
 	if(is.null(ncvregpred[[i]][[1]])){
 		tmp2<-rep(i,1)
 		tmp2<-cbind(tmp2,0)
 		tmp<-rbind(tmp,tmp2)
 	}
 
 }
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),pch=16,xlab='year',ylab='error',cex=0.2, main="scad predictive errors")

scadsize<-matrix(0,nclust,(ny-nyb))
for(jj in 1:nclust){
	for(j in 1:(ny-nyb)){
			if(!is.null(ncvregs[[j]][[jj]])){ 
				lambda<-ncvregs[[j]][[jj]]$lambda[which(ncvregs[[j]][[jj]]$loss==min(ncvregs[[j]][[jj]]$loss))[1]]
				coef<-ncvregs[[j]][[jj]]$beta[,which(ncvregs[[j]][[jj]]$lambda==lambda)]
				scadsize[jj,j]<-length(coef[which(coef>0)])
		}
	}
}
cell<-matrix(0,dim(scadsize)[1],dim(scadsize)[2])
for(jj in 1:nclust){
	for(j in 1:(ny-nyb)){
		if(!is.null(svms[[j]][[jj]])){
			 cell[jj,j]<-scadsize[jj,j]
		}
		if(is.null(trees[[j]][[jj]])){
		 	cell[jj,j]<-"NA"
		}
	}
}
heatmap.2(scadsize,Rowv=FALSE, Colv=FALSE, dendrogram='none',cellnote=cell,col=gray.colors(30),srtCol=360,adjCol=c(0,1),notecol="white",tracecol="white")


tmp4<-NULL
for(i in 1:(ny-nyb)) tmp4<-rbind(tmp4,cbind(tmp3scad[[i]],rep(as.numeric(Year[i+nyb]),dim(tmp3scad[[i]])[1])))
boxplot(abs(tmp4[,1]-tmp4[,2])~tmp4[,3]+tmp4[,4],col=gray(1:6/6),varwidth=T, notch=T, cex=0.4)
# no data for 2006
tmp<-boxplot(abs(tmp4[,1]-tmp4[,2])~as.factor(tmp4[,3])+as.factor(tmp4[,4]),at=rep(8*(0:11),each=6)+rep(1:6,12),col=gray(1:6/6),varwidth=T, notch=T, cex=0.4, boxwex=0.7, axes=FALSE)
axis(1, at=seq(4,92,by=8), labels=Year[(nyb+1):(ny-1)], cex.axis=0.6)
axis(2)
