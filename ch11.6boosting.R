#	ch 11.6 
#	adjunct for boosting

# load libraries, load median workspace, then ... 
test1<-trest
test2<-tresd

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
}
# just run trees ...
for (l in 1:R)
{
		## with ensembling - boosting	
		## because this takes FOREVER, only run 5 at a time
		x<-test1[[l]][,1]
		yt<-test1[[l]][,2]
		yd<-test2[[l]][,2]
		trest[[l]][,1]<-nnrest[[l]][,1]<-tresd[[l]][,1]<-nnresd[[l]][,1]<-x
		trest[[l]][,2]<-nnrest[[l]][,2]<-yt
		tresd[[l]][,2]<-nnresd[[l]][,2]<-yd
	for(k in 1:(N-bn)){
		n<-bn+k-1
		datt<-as.data.frame(cbind(yt,x))
		datd<-as.data.frame(cbind(yd,x))
		## with ensembling boosting- tree and treed
		trControl<-trainControl(method="none")
		tree5t<-train(yt ~ x, data=datt, method="blackboost", subset=1:n, metric="RMSE", tuneLength=1, trControl=trControl) 
		tmodt[[l]][[k]]<-tree5t
		trest[[l]][n+1,7]<-predict(tree5t,datt)[n+1]
		## with ensembling boosting- tree and doppler
		trControl<-trainControl(method="none")
		tree5d<-train(yd ~ x, data=datd, method="blackboost", subset=1:n, metric="RMSE", tuneLength=1, trControl=trControl) 
		tmodd[[l]][[k]]<-tree5d
		tresd[[l]][n+1,7]<-predict(tree5d,datd)[n+1]
	}
}
# just run nns ....
for (l in 1:R)
{
		## with ensembling - boosting	
		## because this takes FOREVER, only run 5 at a time and read in data from median
		## case for runs above 5
		## test1<-trest
		## test2<-tresd
		x<-test1[[l]][,1]
		yt<-test1[[l]][,2]
		yd<-test2[[l]][,2]
		trest[[l]][,1]<-nnrest[[l]][,1]<-tresd[[l]][,1]<-nnresd[[l]][,1]<-x
		trest[[l]][,2]<-nnrest[[l]][,2]<-yt
		tresd[[l]][,2]<-nnresd[[l]][,2]<-yd
		for(k in 1:(N-bn)){
		n<-bn+k-1
		datt<-as.data.frame(cbind(yt,x))
		datd<-as.data.frame(cbind(yd,x))

		## with ensembling boosting- neuralnet and treed
		v=.05 
		mods<-vector("list",dim(mygrid)[1])
		YP<-NULL
		datt3<-datt
		df<-datt3$yt
		for(t in 1:dim(mygrid)[1]){
  			mod<-train(yt ~ x, data=datt3, "nnet", subset=1:n, tuneGrid=mygrid[9,], linout=TRUE, maxit=10000, trace=F)
  			mods[[t]]<-mod
  			yp=predict(mod,newdata=datt3)
  			if(t==1){
  				df=yp
  				YP=yp
  			}
  			if(t>1){
  				df=df - v*yp
  				YP=cbind(YP,v*yp)
  			}
  			datt3$yt<-df
		}
		nnmodt[[l]][[k]]<-mods
		nnrest[[l]][n+1,7]<-apply(YP,1,sum)[n+1]
		## with ensembling boosting- neural net and doppler
		v=.1
		mods<-vector("list",dim(mygrid)[1])
		YP<-NULL
		datd3<-datd
		df<-datd3$yd
		for(t in 1:dim(mygrid)[1]){
  			mod<-train(yd ~ x, data=datd3, "nnet", subset=1:n, tuneGrid=mygrid[9,], linout=TRUE, maxit=10000, trace=F)
  			mods[[t]]<-mod
  			yp=predict(mod,newdata=datd3)
  			if(t==1){
  				df=yp
  				YP=yp
  			}
  			if(t>1){
  				df=df - v*yp
  				YP=cbind(YP,v*yp)
  			}
  			datd3$yd<-df
		}
		nnmodd[[l]][[k]]<-mods
		nnresd[[l]][n+1,7]<-apply(YP,1,sum)[n+1]
	}
}
