#	VegOut data and computing for Chapter 10.6: tree

# Data structure is by location and year, outcome variable is SSG
# See Harms et al. 2009 for details
# Each excel file has different SSG0 depending on the lag time/prediction step; steps are
# 2-week, 4-week, 6-week
# Explanatory variable types: precipitation/drought, satellite, oceanic, biophysical (soil)
# 1402 locations. Biophysical variable do not vary with time (constant within location); 
# oceanic variables do not vary with location (constant within year)
# 18 years of data, 1989-2006

# in Harms et al. predict within location, across years; effective n=18, p=mid-30s

# unsupervised clustering (treating previous responses as explanatory) into 
# 6 clusters based initially on first 5 years of data
# sequential clustering for 1-ith year to predict (i+1)st year
# same prediction function for all locations in a given cluster, but different predictive
# values depending on specific values of explanatory variables
# prediction methods: 
#	ind'l trees

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
# save.image('veg.RData')

## stack data so variables that don't change by location vary
## variables 2-9 do not vary with time so only enter once and repeat
#tmp2<-kpod(datn[,c(2:31,47:68,84:105,121:142,158:179],k=8)
#tmp2$obj_vals
#tmp2$fit_list
varsel<-c(1:9)
#nvc<-22
step<-37
#nyt<-nyb+1
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
	
library(kpodclustr)
library(party)
## make constructs to save results
k<-1
j<-1
nyb<-5
nclust<-6
cluster<-vector("list",length=ny-nyb)
trees<-vector("list",length=ny-nyb)
for(j in 1:(ny-nyb)) trees[[j]]<-vector("list",length=nclust)
tpred<-vector("list",length=ny-nyb-1)
for(j in 1:(ny-nyb)) tpred[[j]]<-vector("list",length=nclust)
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
		c1<-c1[which(!is.na(c1$SSG0_6)),]
		if(dim(c1)[1]>0){		
			tmp<-apply(c1,2,var,use="pairwise.complete.obs")
			c1<-c1[,which(tmp>0)]
			if(dim(c1)[1]>0){
				tree1<-ctree(SSG0_6~., data=c1[c1$Year<as.numeric(Year[nyb2+1]),])
				trees[[k]][[kk]]<-tree1
				ptree1<-predict(tree1, newdata=c1[c1$Year==as.numeric(Year[nyb2+1]),])
				tpred[[k]][[kk]]<-cbind(c1$SSG0_6[c1$Year==as.numeric(Year[nyb2+1])],ptree1)
		#plot(1:length(ptree1), c1$SSG0_6[c1$Year==as.numeric(Year[nyb+1])],type='l')
		#lines(1:length(ptree1),ptree1,col="red")
			}
		}
	}
}
clus<-NULL
for(i in 1:13) clus<-cbind(clus,cluster[[k]])
#write.csv(clus, file="ctree_cluster.csv")

# plot results
tmp3<-vector("list",length=13)
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
tmp<-sort(tpred[[1]][[1]][,1],index.return=T)
plot(1:200,tpred[[1]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)")
points(1:200,tpred[[1]][[1]][tmp$ix,2],col="red")

 
ctreesize<-matrix(0,nclust,(ny-nyb))
for(jj in 1:nclust){
	for(j in 1:(ny-nyb)){
			if(!is.null(trees[[j]][[jj]])){ 
				ctreesize[jj,j]<-length(nodes(trees[[j]][[jj]], unique(where(trees[[j]][[jj]]))))
		}
	}
}
cell<-matrix(0,dim(ctreesize)[1],dim(ctreesize)[2])
for(jj in 1:nclust){
	for(j in 1:(ny-nyb)){
		if(!is.null(trees[[j]][[jj]])){
			 cell[jj,j]<-ctreesize[jj,j]
		}
		if(is.null(trees[[j]][[jj]])){
		 	cell[jj,j]<-"NA"
		}
	}
}
heatmap.2(ctreesize,Rowv=FALSE, Colv=FALSE, dendrogram='none',cellnote=cell,col=gray.colors(30),srtCol=360,adjCol=c(0,1),notecol="white",tracecol="white")

 
