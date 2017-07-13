#	ch 11 computing
#	ensembling: BMA, bagging, stacking, boosting, median methods (2) 
#	components from same or different model classes, mixing strategies
#	
#	simulation with standard functions (e.g., heavisine, blocks, doppler) 
#	components include poly, trees, NN, SVM, RVM, etc. 
#	use individually, ensemble, or ensemble-ensembles

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

## design one-dim tree with breaks at ti and poly of degree d in leaves
treef = function(t,ti,coef)
{
# ti = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)

 y = rep(0, length(t))
 for (j in 1:length(t)){
 	tmp=0
 	var<-c(1,t[j],t[j]^2)
 	k<-1
 	while(tmp==0){
 		if(t[j]<ti[k]){
 			y[j]<-sum(var*coef[k,])
 			tmp=1
 		}
 		k<-k+1
 	}
 }
 return(y) 
}

dopplerf = function(t,e,k)
{
# e = 0.05
 y = rep(0, length(t))
 for (j in 1:length(t))
 {
   a = sqrt(t[j]*(1-t[j]))
  # b = sin(2*pi*(1+e)/(t[j]+e))
   b = sin(k*pi*(1+e)/(t[j]+e))
   y[j] = a*b
 }
 return(y)
}

myweightmedian=function(weight, mass)
{
  z = cbind(weight, mass)
  u = rep(0, length(weight))
  v=u
  u[rank(z[,2])]=z[,1]
  v[rank(z[,2])]=z[,2]
  ## add next 2 lines to get order of models rank sorted
  t = rank(z[,2])
  tsort = sort(t,index.return=TRUE)
  z = cbind(cumsum(u),v)
  weightm = z[z[,1]>=0.5,2]
  ## add next two lines and adjust return
  mod = tsort$ix[z[,1]>=0.5]

  return(c(weightm[1],mod[1]))
}

## sample size: 200 
N<-200
## burn in
bn<-50
## Use more replications, say 50.
R<-30
## parameters of random functions (treed)
tiall<-coefall<-dall<-vector("list", R)
## parameters of random functions (doppler)
eall<-kall<-vector("list", R)
## number of bootstraps
## proportion of data in each bootstrap run
nboot<-100
bootfrac<-0.6

##	treed and doppler, with trees and NNs, no ensemble plus all six ensemble methods
# storage for fitted and predicted values
trest<-nnrest<-tresd<-nnresd<-vector("list",R)
for(i in 1:R){
	trest[[i]]<-tresd[[i]]<-matrix(0,N,8)
	nnrest[[i]]<-nnresd[[i]]<-matrix(0,N,8)
}
# storage for model information
tmodt<-nnmodt<-tmodd<-nnmodd<-vector("list",R)
for(i in 1:R){
	tmodt[[i]]<-tmodd[[i]]<-vector("list",N-bn)
	nnmodt[[i]]<-nnmodd[[i]]<-vector("list",N-bn)
	for(j in 1:(N-bn)){
		tmodt[[i]][[j]]<-tmodd[[i]][[j]]<-vector("list",6)
		nnmodt[[i]][[j]]<-nnmodd[[i]][[j]]<-vector("list",6)
	}	
}
# set of possible model parameter settings for neural networks single layer
mygrid<-expand.grid(.decay=c(0.5, 0.1, 5e-4), .size=c(2,4,6))
mygrid2<-expand.grid(.decay=c(0.5, 0.1, 5e-4), .size=c(2,4,6), .bag=c(4,6,8))

set.seed(9050)

ptm <- proc.time()
for (l in 1:R)
{
	## choose x
	x = runif(N, min=0, max=1)
	sde=1
	error = rnorm(N, sd = sde)

	# randomly generate ti, coef, and d
	# length of ti from [2,10]; values of d from [1,3], values of coef from [-5,5]
	til<-sample(seq(2,10,1),1)
	ti<-runif(til)
	ti<-c(ti,1)
	d<-sample(seq(1,3,1),length(ti),replace=TRUE)
	coef<-matrix(0,length(ti),3)
	for(i in 1:dim(coef)[1]){
		coef[i,1:d[i]]<-runif(d[i],min=-5,max=5)
	}
	yt = treef(x,ti,coef) + error
	tiall[[l]]<-ti
	dall[[l]]<-d
	coefall[[l]]<-coef
## use the commented code below to generate a plot of an example tree function
##	tmp<-sort(x, index.return=T)
##	plot(tmp$x,yt[tmp$ix],type='b', xlab='x', ylab='treed + error')
##	lines(tmp$x,treef(tmp$x,ti,coef),col='red')

	# randomly generate e and k
	# e from [0,0.5]; k from [-10,10]
	e<-sample(seq(0,0.5,0.01),1)
	k1<-runif(1, min=-10, max=10)
	yd = dopplerf(x,e,k1) + error
	eall[l]<-e
	kall[l]<-k1
##	use the commented code below to generate a plot of an example doppler function
##	tmp<-sort(x, index.return=T)
##	plot(tmp$x,yd[tmp$ix],type='b', xlab='x', ylab='doppler + error')
##	lines(tmp$x,dopplerf(tmp$x,e,k1),col='red')
##	library(zoo)
##	rm<-rollmean(dopplerf(tmp$x,e,k1),k=5,na.pad=T)
##	points(tmp$x,rm,pch=12,cex=0.7)

	trest[[l]][,1]<-nnrest[[l]][,1]<-tresd[[l]][,1]<-nnresd[[l]][,1]<-x
	trest[[l]][,2]<-nnrest[[l]][,2]<-yt
	tresd[[l]][,2]<-nnresd[[l]][,2]<-yd
	## without ensembling
	for(k in 1:(N-bn)){
		n<-bn+k-1
		datt<-as.data.frame(cbind(yt,x))
		datd<-as.data.frame(cbind(yd,x))
		## without ensembling - tree and treed
		tree1t<-rpart(yt ~ x, data=datt, subset=1:n, method="anova") 
		opt <- which.min(tree1t$cptable[,"xerror"])
		cp <- tree1t$cptable[opt, "CP"]
		tree1t_prune<-prune(tree1t, cp = cp)
		tmodt[[l]][[k]][[1]]<-tree1t_prune
		trest[[l]][n+1,3]<-predict(tree1t_prune,newdata=datt,type="vector")[n+1]
		if(k==1) trest[[l]][1:n,3]<-predict(tree1t_prune,newdata=datt,type="vector")[1:n]
		## without ensembling - tree and doppler
		tree1d<-rpart(yd ~ x, data=datd, subset=1:n, method="anova") 
		opt <- which.min(tree1d$cptable[,"xerror"])
		cp <- tree1d$cptable[opt, "CP"]
		tree1d_prune<-prune(tree1d, cp = cp)
		tmodd[[l]][[k]][[1]]<-tree1d_prune
		tresd[[l]][n+1,3]<-predict(tree1d_prune,newdata=datd,type="vector")[n+1]
		if(k==1) tresd[[l]][1:n,3]<-predict(tree1d_prune,newdata=datd,type="vector")[1:n]
		## without ensembling - neural net and treed
		nn1t<-train(yt ~ x, data=datt, "nnet", subset=1:n, tuneGrid=mygrid, linout=TRUE, maxit=10000, trace=F)
		nnmodt[[l]][[k]][[1]]<-nn1t
		nnrest[[l]][n+1,3]<-predict(nn1t$finalModel, newdata=datt)[n+1]
		if(k==1) nnrest[[l]][1:n,3]<-predict(nn1t$finalModel,newdata=datt)[1:n]
		## without ensembling - neural net and doppler
		nn1d<-train(yd ~ x, data=datd, "nnet", subset=1:n, tuneGrid=mygrid, linout=TRUE, maxit=10000, trace=F)
		nnmodd[[l]][[k]][[1]]<-nn1d
		nnresd[[l]][n+1,3]<-predict(nn1d$finalModel, newdata=datd)[n+1]
		if(k==1) nnresd[[l]][1:n,3]<-predict(nn1d$finalModel, newdata=datd)[1:n]
	}
		## with ensembling: BMA
	for(k in 1:(N-bn)){
		n<-bn+k-1
		datt<-as.data.frame(cbind(yt,x))
		datd<-as.data.frame(cbind(yd,x))
		## with ensembling BMA- tree and treed
		tree2t<-bart(x[1:n], yt[1:n], x.test=x[1:(n+1)], ntree=100) 
		tmodt[[l]][[k]][[2]]<-tree2t
		trest[[l]][n+1,4]<-tree2t$yhat.test.mean[n+1]
		## with ensembling BMA - tree and doppler
		tree2d<-bart(x[1:n], yd[1:n], x.test=x[1:(n+1)], ntree=100) 
		tmodd[[l]][[k]][[2]]<-tree2d
		tresd[[l]][n+1,4]<-tree2d$yhat.test.mean[n+1]
		## with ensembling BMA - neural net and treed
		## assumes equal prior probability each model and equal conditional prior probability 
		## each child model given the base models; see Chitsazan 2015
		## weighted sum of squared errors (Qp,q) and the BIC complexity term (N ln 2π + mp ln N)
		### use e^{-0.5BIC}/sum[e^{-0.5BIC}] as weights
		post<-pred<-NULL
		for(j in 1:dim(mygrid)[1]){
			nn2t<-train(yt ~ x, data=datt, "nnet", subset=1:n, tuneGrid=mygrid[j,], linout=TRUE, maxit=10000, trace=F)
			post<-c(post,exp(-0.5*(nn2t$results[3]+n*log(2*pi)+(1+mygrid[j,2])*log(n))))
			pred<-c(pred,predict(nn2t$finalModel, newdata=datt)[n+1])
		}
				
		post<-unlist(post)/sum(unlist(post))
		nnmodt[[l]][[k]][[2]]<-cbind(post,pred)
		nnrest[[l]][n+1,4]<-sum(post*pred)
		## with ensembling BMA - neural net and doppler
		post<-pred<-NULL
		for(j in 1:dim(mygrid)[1]){
			nn2d<-train(yd ~ x, data=datd, "nnet", subset=1:n, tuneGrid=mygrid[j,], linout=TRUE, maxit=10000, trace=F)
			post<-c(post,nn2d$results[3]+n*log(2*pi)+(1+mygrid[j,2])*log(n))
			pred<-c(pred,predict(nn2d$finalModel, newdata=datd)[n+1])
		}
		post<-unlist(post)/sum(unlist(post))
		nnmodd[[l]][[k]][[2]]<-cbind(post,pred)
		nnresd[[l]][n+1,4]<-sum(post*pred)
	}
		## with ensembling: bagging
	for(k in 1:(N-bn)){
		n<-bn+k-1
		datt<-as.data.frame(cbind(yt,x))
		datd<-as.data.frame(cbind(yd,x))
		## with ensembling bagging- tree and treed
		trControl<-trainControl(method="oob", selectionFunction="oneSE")
		tree3t<-train(yt ~ x, data=datt, method="treebag", subset=1:n, metric="RMSE", keepX=TRUE, nbagg=10, trControl=trControl) 
		tmodt[[l]][[k]][[3]]<-tree3t
		trest[[l]][n+1,5]<-predict(tree3t,datt)[n+1]
		## with ensembling bagging- tree and doppler
		tree3d<-train(yd ~ x, data=datd, method="treebag", subset=1:n, metric="RMSE", keepX=TRUE, nbagg=10, trControl=trControl) 
		tmodd[[l]][[k]][[3]]<-tree3d
		tresd[[l]][n+1,5]<-predict(tree3d,datd)[n+1]
		## with ensembling bagging- neural net and treed
		nn3t<-train(yt ~ x, data=datt, method="avNNet", subset=1:n, tuneGrid=mygrid2, linout=TRUE, maxit=10000, trace=F) 
		nnmodt[[l]][[k]][[3]]<-nn3t
		nnrest[[l]][n+1,5]<-predict(nn3t,datt)[n+1]
		## with ensembling bagging- neural net and doppler
		nn3d<-train(yd ~ x, data=datd, method="avNNet", subset=1:n, tuneGrid=mygrid2, linout=TRUE, maxit=10000, trace=F) 
		nnmodd[[l]][[k]][[3]]<-nn3d
		nnresd[[l]][n+1,5]<-predict(nn3d,datd)[n+1]
	}
		## with ensembling: stacking
	for(k in 1:(N-bn)){
		n<-bn+k-1
		datt<-as.data.frame(cbind(yt,x))
		datd<-as.data.frame(cbind(yd,x))
		## with ensembling stacking- tree and treed
		tmp<-NULL
		mod1<-rpart(yt ~ x, data=datt, subset=1:n, method="anova")
		tmp<-cbind(tmp,predict(mod1,datt)[1:(n+1)])
		mod2<-rpart(yt ~ x, data=datt, subset=1:n, method="anova", control=rpart.control(minsplit=10, minbucket=3, cp=0.05))
		if(!identical(predict(mod1,datt)[1:n],predict(mod2,datt)[1:n])) tmp<-cbind(tmp,predict(mod2,datt)[1:(n+1)])
		mod3<-rpart(yt ~ x, data=datt, subset=1:n, method="anova", control=rpart.control(minsplit=5, minbucket=3))
		if(sum(as.numeric(apply(tmp,2,function(x){identical(predict(mod3,datt)[1:n],x)})))<1) tmp<-cbind(tmp,predict(mod3,datt)[1:(n+1)])
		mod4<-ctree(yt ~ x, data=datt, subset=1:n)
		if(sum(as.numeric(apply(tmp,2,function(x){identical(predict(mod4,datt)[1:n],x)})))<1) tmp<-cbind(tmp,predict(mod4,datt)[1:(n+1)])
		mod5<-ctree(yt ~ x, data=datt, subset=1:n, controls=ctree_control(stump=TRUE))
		if(sum(as.numeric(apply(tmp,2,function(x){identical(predict(mod5,datt)[1:n],x)})))<1) tmp<-cbind(tmp,predict(mod5,datt)[1:(n+1)])
		#
		datt2<-as.data.frame(cbind(yt[1:(n+1)],tmp))
		tree4t<-lm(V1 ~ ., data=datt2, subset=1:n) 
		tmodt[[l]][[k]]<-tree4t
		trest[[l]][n+1,6]<-predict(tree4t,datt2)[n+1]
		## with ensembling stacking- tree and doppler
		tmp<-NULL
		mod1<-rpart(yd ~ x, data=datd, subset=1:n, method="anova")
		tmp<-cbind(tmp,predict(mod1,datd)[1:(n+1)])
		mod2<-rpart(yd ~ x, data=datd, subset=1:n, method="anova", control=rpart.control(minsplit=10, minbucket=3, cp=0.05))
		if(!identical(predict(mod1,datd)[1:n],predict(mod2,datd)[1:n])) tmp<-cbind(tmp,predict(mod2,datd)[1:(n+1)])
		mod3<-rpart(yd ~ x, data=datd, subset=1:n, method="anova", control=rpart.control(minsplit=5, minbucket=3))
		if(sum(as.numeric(apply(tmp,2,function(x){identical(predict(mod3,datt)[1:n],x)})))<1) tmp<-cbind(tmp,predict(mod3,datd)[1:(n+1)])
		mod4<-ctree(yd ~ x, data=datd, subset=1:n)
		if(length(nodes(mod4, unique(where(mod4))))>1) if(sum(as.numeric(apply(tmp,2,function(x){identical(predict(mod4,datt)[1:n],x)})))<1) tmp<-cbind(tmp,predict(mod4,datd)[1:(n+1)])
		mod5<-ctree(yd ~ x, data=datd, subset=1:n, controls=ctree_control(stump=TRUE))
		if(length(nodes(mod5, unique(where(mod5))))>1) if(sum(as.numeric(apply(tmp,2,function(x){identical(predict(mod5,datt)[1:n],x)})))<1) tmp<-cbind(tmp,predict(mod5,datd)[1:(n+1)])
		#
		datd2<-as.data.frame(cbind(yd[1:(n+1)],tmp))
		tree4d<-lm(V1 ~ ., data=datd2, subset=1:n) 
		tmodd[[l]][[k]]<-tree4d
		tresd[[l]][n+1,6]<-predict(tree4d,datd2)[n+1]
		## with ensembling stacking- neuralnet and treed
		tmp<-NULL
		for(j in 1:dim(mygrid)[1]){
			mod<-train(yt ~ x, data=datt, "nnet", subset=1:n, tuneGrid=mygrid[j,], linout=TRUE, maxit=10000, trace=F)
			if(j<2) tmp<-cbind(tmp,predict(mod,datt)[1:(n+1)])
			if(j>1) if(sum(as.numeric(apply(tmp,2,function(x){identical(predict(mod,datt)[1:n],x)})))<1) tmp<-cbind(tmp,predict(mod,datt)[1:(n+1)])
		}
		datt2<-as.data.frame(cbind(yt[1:(n+1)],tmp))
		nn4t<-lm(V1 ~ ., data=datt2, subset=1:n) 
		nnmodt[[l]][[k]]<-nn4t
		nnrest[[l]][n+1,6]<-predict(nn4t,datt2)[n+1]
		## with ensembling stacking- neuralnet and doppler
		tmp<-NULL
		for(j in 1:dim(mygrid)[1]){
			mod<-train(yd ~ x, data=datd, "nnet", subset=1:n, tuneGrid=mygrid[j,], linout=TRUE, maxit=10000, trace=F)
			if(j<2) tmp<-cbind(tmp,predict(mod,datd)[1:(n+1)])
			if(j>1) if(sum(as.numeric(apply(tmp,2,function(x){identical(predict(mod,datd)[1:n],x)})))<1) tmp<-cbind(tmp,predict(mod,datd)[1:(n+1)])
		}
		datd2<-as.data.frame(cbind(yd[1:(n+1)],tmp))
		nn4d<-lm(V1 ~ ., data=datd2, subset=1:n) 
		nnmodd[[l]][[k]]<-nn4d
		nnresd[[l]][n+1,6]<-predict(nn4d,datd2)[n+1]
	}

		## with ensembling - boosting	
		## because this takes FOREVER, only run 5 at a time and read in data from median
		## see ch11.6boosting.R for code
		
		## with ensembling - median method
	for(k in 1:(N-bn)){
		n<-bn+k-1
		datt<-as.data.frame(cbind(yt,x))
		datd<-as.data.frame(cbind(yd,x))
		## with ensembling median - tree and treed
		## approx BIC with N/sig*[mse + logN*(d/N)*sig]
		## use e^{-0.5BIC}/sum[e^{-0.5BIC}] as weights
		rp.ctrl<-expand.grid(minsplit=c(5,10,15), cp=c(0.01, 0.05), xval=c(3,5,10), maxdepth=c(2,6))
		post<-pred<-NULL
		for(kk in 1:dim(rp.ctrl)[1]){
			tree6t<-rpart(yd ~ x, data=datt, subset=1:n, method="anova", control=rpart.control(minsplit=rp.ctrl[kk,1], cp=rp.ctrl[kk,2], xval=rp.ctrl[kk,3], maxdepth=rp.ctrl[kk,4]))
			mse<-sum(tree6t$frame[(tree6t$frame$var=="<leaf>"),]$dev/tree6t$frame[(tree6t$frame$var=="<leaf>"),]$wt)
			post<-c(post,exp(-0.5*((n*mse)+log(n)*(dim(tree6t$frame)[1]/n))))
			pred<-c(pred,predict(tree6t, newdata=datt)[n+1])
		}
		post<-unlist(post)/sum(unlist(post))
		tmodt[[l]][[k]]<-cbind(post,pred)
		## select models
		tmp = myweightmedian(post, pred)
		wmBICnewy = tmp[1]
		## compute median prediction
		trest[[l]][n+1,8]<-wmBICnewy
		## with ensembling median - tree and doppler
		post<-pred<-NULL
		for(kk in 1:dim(rp.ctrl)[1]){
			tree6d<-rpart(yd ~ x, data=datd, subset=1:n, method="anova", control=rpart.control(minsplit=rp.ctrl[kk,1], cp=rp.ctrl[kk,2], xval=rp.ctrl[kk,3], maxdepth=rp.ctrl[kk,4]))
			mse<-sum(tree6d$frame[(tree6d$frame$var=="<leaf>"),]$dev/tree6d$frame[(tree6d$frame$var=="<leaf>"),]$wt)
			post<-c(post,exp(-0.5*((n*mse)+log(n)*(dim(tree6d$frame)[1]/n))))
			pred<-c(pred,predict(tree6d, newdata=datd)[n+1])
		}
		post<-unlist(post)/sum(unlist(post))
		tmodd[[l]][[k]]<-cbind(post,pred)
		## select models
		tmp = myweightmedian(post, pred)
		wmBICnewy = tmp[1]
		## compute median prediction
		tresd[[l]][n+1,8]<-wmBICnewy
		## with ensembling median - neural net and treed
		post<-pred<-NULL
		for(j in 1:dim(mygrid)[1]){
			nn6t<-train(yt ~ x, data=datt, "nnet", subset=1:n, tuneGrid=mygrid[j,], linout=TRUE, maxit=10000, trace=F)
			post<-c(post,exp(-0.5*(nn6t$results[3]+n*log(2*pi)+(1+mygrid[j,2])*log(n))))
			pred<-c(pred,predict(nn6t$finalModel, newdata=datt)[n+1])
		}				
		post<-unlist(post)/sum(unlist(post))
		nnmodt[[l]][[k]]<-cbind(post,pred)
		## select models
		tmp = myweightmedian(post, pred)
		wmBICnewy = tmp[1]
		## compute median prediction
		nnrest[[l]][n+1,8]<-wmBICnewy
		## with ensembling median - neural net and doppler
		post<-pred<-NULL
		for(j in 1:dim(mygrid)[1]){
			nn6d<-train(yd ~ x, data=datd, "nnet", subset=1:n, tuneGrid=mygrid[j,], linout=TRUE, maxit=10000, trace=F)
			post<-c(post,exp(-0.5*(nn6d$results[3]+n*log(2*pi)+(1+mygrid[j,2])*log(n))))
			pred<-c(pred,predict(nn6d$finalModel, newdata=datt)[n+1])
		}				
		post<-unlist(post)/sum(unlist(post))
		nnmodd[[l]][[k]]<-cbind(post,pred)
		## select models
		tmp = myweightmedian(post, pred)
		wmBICnewy = tmp[1]
		## compute median prediction
		nnresd[[l]][n+1,8]<-wmBICnewy
	}	
}

