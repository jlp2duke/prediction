#	VegOut data and plotting for Chapter 10

# plot results for all five methods: tree, nnet, rvm, svm, ncvreg
# assume all methods have been run and results saved in the appropriate constructs
 
par(mfrow=c(3,2))
tmp<-sort(tpred[[12]][[1]][,1],index.return=T)
plot(1:dim(tpred[[12]][[1]])[1],tpred[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", main="tree predictions")
points(1:dim(tpred[[12]][[1]])[1],tpred[[12]][[1]][tmp$ix,2],col="red")

tmp<-sort(nnetpred[[12]][[1]][,1],index.return=T)
plot(1:dim(nnetpred[[12]][[1]])[1],nnetpred[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", main="neural net predictions")
points(1:dim(nnetpred[[12]][[1]])[1],nnetpred[[12]][[1]][tmp$ix,2],col="red")

tmp<-sort(rvmpred[[12]][[1]][,1],index.return=T)
plot(1:dim(rvmpred[[12]][[1]])[1],rvmpred[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", main="rvm predictions")
points(1:dim(rvmpred[[12]][[1]])[1],rvmpred[[12]][[1]][tmp$ix,2],col="red")

tmp<-sort(svmpred[[12]][[1]][,1],index.return=T)
plot(1:dim(svmpred[[12]][[1]])[1],svmpred[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", main="svm predictions")
points(1:dim(svmpred[[12]][[1]])[1],svmpred[[12]][[1]][tmp$ix,2],col="red")

tmp<-sort(ncvregpred[[12]][[1]][,1],index.return=T)
plot(1:dim(ncvregpred[[12]][[1]])[1],ncvregpred[[12]][[1]][tmp$ix,1],xlab="cluster 1 obs",ylab="response (ordered)", main="SCAD predictions")
points(1:dim(ncvregpred[[12]][[1]])[1],ncvregpred[[12]][[1]][tmp$ix,2],col="red")

# beeswarm plot
library(beeswarm)
tmp<-NULL
for(i in 1:13){
	if(!is.null(tpred[[i]][[1]])){
		tmp2<-rep(i,length(tpred[[i]][[1]][,1]))
		tmp2<-cbind(tmp2,(tpred[[i]][[1]][,1]-tpred[[i]][[1]][,2]))
		tmp<-rbind(tmp,tmp2)
	}
}
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),pch=16,xlab='',ylab='error',cex=0.5) 

par(mfrow=c(2,3))
library(beeswarm)
tmp<-NULL
for(i in 1:13){
 	if(!is.null(tpred[[i]][[1]])){
 		tmp2<-rep(i,length(tpred[[i]][[1]][,1]))
  		tmp2<-cbind(tmp2,(tpred[[i]][[1]][,1]-tpred[[i]][[1]][,2]))
 		tmp<-rbind(tmp,tmp2)
 	}
 }
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300), pch=16,xlab='year',ylab='error',cex=0.2, main="tree predictive errors")

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
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300),pch=16,xlab='year',ylab='error',cex=0.2, main="neural net predictive errors") 

tmp<-NULL
for(i in 1:13){
 	if(!is.null(svmpred[[i]][[1]])){
 		tmp2<-rep(i,length(svmpred[[i]][[1]][,1]))
 		tmp2<-cbind(tmp2,(svmpred[[i]][[1]][,1]-svmpred[[i]][[1]][,2]))
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
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300), pch=16,xlab='year',ylab='error',cex=0.2, main="svm predictive errors")
tmp<-NULL
for(i in 1:13){
 	if(!is.null(rvmpred[[i]][[1]])){
 		tmp2<-rep(i,length(rvmpred[[i]][[1]][,1]))
 		tmp2<-cbind(tmp2,(rvmpred[[i]][[1]][,1]-rvmpred[[i]][[1]][,2]))
 		tmp<-rbind(tmp,tmp2)
 	}
 	if(is.null(rvmpred[[i]][[1]])){
 		tmp2<-rep(i,1)
 		tmp2<-cbind(tmp2,0)
 		tmp<-rbind(tmp,tmp2)
 	}
 
 }
tmp<-as.data.frame(tmp)
colnames(tmp)<-c('year','error')
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300),pch=16,xlab='year',ylab='error',cex=0.2, main="rvm predictive errors")

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
beeswarm(tmp$error~tmp$year, method=c("swarm"),ylim=c(-250,300),pch=16,xlab='year',ylab='error',cex=0.2, main="scad predictive errors")

## boxplots of errors
#require(RCurl)
#source(textConnection(getURL("https://gist.github.com/mages/5339689/raw/576263b8f0550125b61f4ddba127f5aa00fa2014/add.alpha.R")))
#http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

myColours = c("lightblue", "lightgreen")
myColoursAlpha<-add.alpha(myColours, alpha=0.5)
boxplot(abs(tmp3ctree[[1]][,1]-tmp3ctree[[1]][,2])~tmp3ctree[[1]][,3],col=myColoursAlpha[1])
for(i in 2:(ny-nyb-1)) boxplot(abs(tmp3ctree[[i]][,1]-tmp3ctree[[i]][,2])~tmp3ctree[[i]][,3],add=T,axes=F)
boxplot(abs(tmp3ctree[[ny-nyb]][,1]-tmp3ctree[[ny-nyb]][,2])~tmp3ctree[[ny-nyb]][,3],col=myColoursAlpha[2],add=T,axes=F)

tmp4<-NULL
for(i in 1:(ny-nyb)) tmp4<-rbind(tmp4,cbind(tmp3ctree[[i]],rep(as.numeric(Year[i+nyb]),dim(tmp3ctree[[i]])[1])))
boxplot(abs(tmp4[,1]-tmp4[,2])~tmp4[,3]+tmp4[,4],col=gray(1:6/6),varwidth=T, notch=T, cex=0.4)

tmp<-boxplot(abs(tmp4[,1]-tmp4[,2])~as.factor(tmp4[,3])+as.factor(tmp4[,4]),col=rainbow(6), varwidth=T, notch=T, cex=0.4, boxwex=0.7, las=2, cex.axis=0.5)
tmp<-boxplot(abs(tmp4[,1]-tmp4[,2])~as.factor(tmp4[,3])+as.factor(tmp4[,4]),col=gray(1:6/6), names=rep(c('1','2','3','4','5','6'),13),varwidth=T, notch=T, cex=0.4, boxwex=0.7, las=2, cex.axis=0.5)
tmp<-boxplot(abs(tmp4[,1]-tmp4[,2])~as.factor(tmp4[,3])+as.factor(tmp4[,4]),col=gray(1:6/6),varwidth=T, notch=T, cex=0.4, boxwex=0.7, las=2, cex.axis=0.5)
tmp<-boxplot(abs(tmp4[,1]-tmp4[,2])~as.factor(tmp4[,3])+as.factor(tmp4[,4]),at=rep(8*(0:12),each=6)+rep(1:6,13),col=gray(1:6/6),varwidth=T, notch=T, cex=0.4, boxwex=0.7, cex.axis=0.5, names=NULL)
tmp<-boxplot(abs(tmp4[,1]-tmp4[,2])~as.factor(tmp4[,3])+as.factor(tmp4[,4]),at=rep(8*(0:12),each=6)+rep(1:6,13),col=gray(1:6/6),varwidth=T, notch=T, cex=0.4, boxwex=0.7, axes=FALSE)
axis(1, at=seq(4,100,by=8), labels=Year[(nyb+1):ny], cex.axis=0.6)
axis(2)

tmp4ctree<-tmp4

tmp5<-NULL
tmp5<-rbind(tmp5,cbind(tmp4ctree,rep(1,dim(tmp4ctree)[1])))
tmp5<-rbind(tmp5,cbind(tmp4nnet,rep(2,dim(tmp4nnet)[1])))
tmp5<-rbind(tmp5,cbind(tmp4rvm,rep(3,dim(tmp4rvm)[1])))
tmp5<-rbind(tmp5,cbind(tmp4scad,rep(4,dim(tmp4scad)[1])))
tmp5<-rbind(tmp5,cbind(tmp4svm,rep(5,dim(tmp4svm)[1])))
tmp<-boxplot(abs(tmp5[,1]-tmp5[,2])~as.factor(tmp5[,5])+as.factor(tmp5[,4]),at=rep(8*(0:12),each=5)+rep(1:5,13),col=gray(1:5/5),varwidth=T, notch=T, cex=0.4, boxwex=0.7, axes=FALSE)
axis(1, at=seq(2,98,by=8), labels=Year[(nyb+1):ny], cex.axis=0.6)
axis(2)
legend('topright',legend=c('ctree','neuralnet','rvm','scad','svm'),col=gray(1:5/5), lty=1, lwd=4)
