#	ch11.6 code for making result plots
# 	assume that code in ch11.6computing.R and ch11.6boosting.R have been run and 
#	results of each method saved in different workspaces

## examine results of simulation studies
##	sanity check; plots; summaries

## for boosting ...
load("/Book/ch11/sim_boost_1-5.RData")
l<-5
n<-199
tmp<-sort(tresd[[l]][,1], index.return=T)
plot(tmp$x, tresd[[l]][tmp$ix,2],type="b", xlab="x", ylab="doppler")
# why column 3 and column 7 is empty?
lines(tmp$x, tresd[[l]][tmp$ix,3],col="red")
lines(tmp$x, nnresd[[l]][tmp$ix,3],col="blue")
legend('bottomright', legend=c('data','tree','nn'),col=c('black','red','blue'),lty=1)

## bagging
load("Book/ch11/sim_bagging_res.RData")
l<-1
n<-199
tmp<-sort(tresd[[l]][,1], index.return=T)
plot(tmp$x, tresd[[l]][tmp$ix,2],type="b", xlab="x", ylab="doppler")
lines(tmp$x, tresd[[l]][tmp$ix,5],col="red")
lines(tmp$x, nnresd[[l]][tmp$ix,5],col="blue")
legend('topright', legend=c('data','tree','nn'),col=c('black','red','blue'),lty=1)

## BMA
load("Book/ch11/sim_bma_res.RData")
l<-1
n<-199
tmp<-sort(tresd[[l]][,1], index.return=T)
plot(tmp$x, tresd[[l]][tmp$ix,2],type="b", xlab="x", ylab="doppler")
lines(tmp$x, tresd[[l]][tmp$ix,4],col="red")
lines(tmp$x, nnresd[[l]][tmp$ix,4],col="blue")
legend('topright', legend=c('data','tree','nn'),col=c('black','red','blue'),lty=1)

## generate error plots with bounds
load("Book/ch11/sim_median.RData")
errtrest<-matrix(999,200,30)
for(i in 1:R) errtrest[,i]<-abs(trest[[i]][,2]-trest[[i]][,8])
sdtrest<-apply(errtrest,1,function(x){sqrt(var(x))})
errnnrest<-matrix(999,200,30)
for(i in 1:R) errnnrest[,i]<-abs(nnrest[[i]][,2]-nnrest[[i]][,8])
sdnnrest<-apply(errnnrest,1,function(x){sqrt(var(x))})
errtresthi<-apply(errtrest,1,mean)[51:200]+sdtrest[51:200]
errtrestlo<-apply(errtrest,1,mean)[51:200]-sdtrest[51:200]
errnnresthi<-apply(errnnrest,1,mean)[51:200]+sdnnrest[51:200]
errnnrestlo<-apply(errnnrest,1,mean)[51:200]-sdnnrest[51:200]
plot(51:200,apply(errtrest,1,mean)[51:200],type="l",ylim=c(-2,6), main="median predictive error: treed", xlab="x",ylab="average predictive L1 error +/- SD")
polygon(c(51:200,rev(51:200)),c(errtresthi,rev(errtrestlo)),col='grey',border=T,lty=3)
lines(51:200,apply(errtrest,1,mean)[51:200])
lines(51:200,apply(errnnrest,1,mean)[51:200])
polygon(c(51:200,rev(51:200)),c(errnnresthi,rev(errnnrestlo)),col=rgb(0,0,0,max=255,alpha=125),border=T,lty=2)
lines(51:200,apply(errnnrest,1,mean)[51:200])
legend('bottomright',legend=c('treed','neural net'),lty=c(1,6))
# save as median_treed_error_mean+sd

errtresd<-matrix(999,200,30)
for(i in 1:R) errtresd[,i]<-abs(tresd[[i]][,2]-tresd[[i]][,8])
sdtresd<-apply(errtresd,1,function(x){sqrt(var(x))})
errnnresd<-matrix(999,200,30)
for(i in 1:R) errnnresd[,i]<-abs(nnresd[[i]][,2]-nnresd[[i]][,8])
sdnnresd<-apply(errnnresd,1,function(x){sqrt(var(x))})
errtresdhi<-apply(errtresd,1,mean)[51:200]+sdtresd[51:200]
errtresdlo<-apply(errtresd,1,mean)[51:200]-sdtresd[51:200]
errnnresdhi<-apply(errnnresd,1,mean)[51:200]+sdnnresd[51:200]
errnnresdlo<-apply(errnnresd,1,mean)[51:200]-sdnnresd[51:200]
plot(51:200,apply(errtresd,1,mean)[51:200],type="l",ylim=c(-2,4), main="median predictive error: doppler", xlab="x",ylab="average predictive L1 error +/- SD")
polygon(c(51:200,rev(51:200)),c(errtresdhi,rev(errtresdlo)),col='grey',border=T,lty=3)
lines(51:200,apply(errtresd,1,mean)[51:200])
lines(51:200,apply(errnnresd,1,mean)[51:200],lty=6)
polygon(c(51:200,rev(51:200)),c(errnnresdhi,rev(errnnresdlo)),col=rgb(0,0,0,max=255,alpha=125),border=T,lty=2)
lines(51:200,apply(errnnresd,1,mean)[51:200],lty=6)
legend('bottomright',legend=c('treed','neural net'),lty=c(1,6))
# save as median_doppler_error_mean+sd

par(mfrow=c(2,2))
boxplot(errtrest, main="median predictive error: trees and treed",cex.main=0.6,ylim=c(0,9))
boxplot(errnnrest, main="median predictive error: neural net and treed",cex.main=0.6,ylim=c(0,9))
boxplot(errtresd, main="median predictive error: trees and doppler",cex.main=0.6,ylim=c(0,5))
boxplot(errnnresd, main="median predictive error: neural net and doppler",cex.main=0.6,ylim=c(0,5))
# save as boxplot_median_errors

load("Book/ch11/sim_bma_res.RData")
R<-30
errtrest<-matrix(999,200,30)
for(i in 1:R) errtrest[,i]<-abs(trest[[i]][,2]-trest[[i]][,4])
sdtrest<-apply(errtrest,1,function(x){sqrt(var(x))})
errnnrest<-matrix(999,200,30)
for(i in 1:R) errnnrest[,i]<-abs(nnrest[[i]][,2]-nnrest[[i]][,4])
sdnnrest<-apply(errnnrest,1,function(x){sqrt(var(x))}) 
errtresthi<-apply(errtrest,1,mean)[51:200]+sdtrest[51:200]
errtrestlo<-apply(errtrest,1,mean)[51:200]-sdtrest[51:200]
errnnresthi<-apply(errnnrest,1,mean)[51:200]+sdnnrest[51:200]
errnnrestlo<-apply(errnnrest,1,mean)[51:200]-sdnnrest[51:200]
plot(51:200,apply(errtrest,1,mean)[51:200],type="l",ylim=c(-2,8), main="bma predictive error: treed", xlab="x",ylab="average predictive L1 error +/- SD")
polygon(c(51:200,rev(51:200)),c(errtresthi,rev(errtrestlo)),col='grey',border=T,lty=3)
lines(51:200,apply(errtrest,1,mean)[51:200])
lines(51:200,apply(errnnrest,1,mean)[51:200],lty=6)
polygon(c(51:200,rev(51:200)),c(errnnresthi,rev(errnnrestlo)),col=rgb(0,0,0,max=255,alpha=125),border=T,lty=2)
lines(51:200,apply(errnnrest,1,mean)[51:200],lty=6)
legend('bottomright',legend=c('treed','neural net'),lty=c(1,6))
# save as bma_treed_error_mean+sd

errtresd<-matrix(999,200,30)
for(i in 1:R) errtresd[,i]<-abs(tresd[[i]][,2]-tresd[[i]][,4])
sdtresd<-apply(errtresd,1,function(x){sqrt(var(x))})
errnnresd<-matrix(999,200,30)
for(i in 1:R) errnnresd[,i]<-abs(nnresd[[i]][,2]-nnresd[[i]][,4])
sdnnresd<-apply(errnnresd,1,function(x){sqrt(var(x))})
errtresdhi<-apply(errtresd,1,mean)[51:200]+sdtresd[51:200]
errtresdlo<-apply(errtresd,1,mean)[51:200]-sdtresd[51:200]
errnnresdhi<-apply(errnnresd,1,mean)[51:200]+sdnnresd[51:200]
errnnresdlo<-apply(errnnresd,1,mean)[51:200]-sdnnresd[51:200]
plot(51:200,apply(errtresd,1,mean)[51:200],type="l",ylim=c(-2,4), main="bma predictive error: doppler", xlab="x",ylab="average predictive L1 error +/- SD")
polygon(c(51:200,rev(51:200)),c(errtresdhi,rev(errtresdlo)),col='grey',border=T,lty=3)
lines(51:200,apply(errtresd,1,mean)[51:200])
lines(51:200,apply(errnnresd,1,mean)[51:200],lty=6)
polygon(c(51:200,rev(51:200)),c(errnnresdhi,rev(errnnresdlo)),col=rgb(0,0,0,max=255,alpha=125),border=T,lty=2)
lines(51:200,apply(errnnresd,1,mean)[51:200],lty=6)
legend('bottomright',legend=c('treed','neural net'),lty=c(1,6))
# save as bma_doppler_error_mean+sd

par(mfrow=c(2,2))
boxplot(errtrest, main="bma predictive error: trees and treed",cex.main=0.6,ylim=c(0,16))
boxplot(errnnrest, main="bma predictive error: neural net and treed",cex.main=0.6,ylim=c(0,16))
boxplot(errtresd, main="bma predictive error: trees and doppler",cex.main=0.6,ylim=c(0,4))
boxplot(errnnresd, main="bma predictive error: neural net and doppler",cex.main=0.6,ylim=c(0,4))
# save as boxplot_bma_errors

load("Book/ch11/sim_noens_res.RData")
R<-30
errtrest<-matrix(999,200,30)
for(i in 1:R) errtrest[,i]<-abs(trest[[i]][,2]-trest[[i]][,3])
sdtrest<-apply(errtrest,1,function(x){sqrt(var(x))})
errnnrest<-matrix(999,200,30)
for(i in 1:R) errnnrest[,i]<-abs(nnrest[[i]][,2]-nnrest[[i]][,3])
sdnnrest<-apply(errnnrest,1,function(x){sqrt(var(x))}) 
errtresthi<-apply(errtrest,1,mean)[51:200]+sdtrest[51:200]
errtrestlo<-apply(errtrest,1,mean)[51:200]-sdtrest[51:200]
errnnresthi<-apply(errnnrest,1,mean)[51:200]+sdnnrest[51:200]
errnnrestlo<-apply(errnnrest,1,mean)[51:200]-sdnnrest[51:200]
plot(51:200,apply(errtrest,1,mean)[51:200],type="l",ylim=c(-2,4), main="nonensembled predictive error: treed", xlab="n",ylab="average predictive L1 error +/- SD")
polygon(c(51:200,rev(51:200)),c(errtresthi,rev(errtrestlo)),col='grey',border=T,lty=3)
lines(51:200,apply(errtrest,1,mean)[51:200])
lines(51:200,apply(errnnrest,1,mean)[51:200],lty=6)
polygon(c(51:200,rev(51:200)),c(errnnresthi,rev(errnnrestlo)),col=rgb(0,0,0,max=255,alpha=125),border=T,lty=2)
lines(51:200,apply(errnnrest,1,mean)[51:200],lty=6)
 legend('bottomright',legend=c('trees','neural net'),lty=c(1,6))
# save as noens_treed_error_mean+sd

errtresd<-matrix(999,200,30)
for(i in 1:R) errtresd[,i]<-abs(tresd[[i]][,2]-tresd[[i]][,3])
sdtresd<-apply(errtresd,1,function(x){sqrt(var(x))})
errnnresd<-matrix(999,200,30)
for(i in 1:R) errnnresd[,i]<-abs(nnresd[[i]][,2]-nnresd[[i]][,3])
sdnnresd<-apply(errnnresd,1,function(x){sqrt(var(x))})
errtresdhi<-apply(errtresd,1,mean)[51:200]+sdtresd[51:200]
errtresdlo<-apply(errtresd,1,mean)[51:200]-sdtresd[51:200]
errnnresdhi<-apply(errnnresd,1,mean)[51:200]+sdnnresd[51:200]
errnnresdlo<-apply(errnnresd,1,mean)[51:200]-sdnnresd[51:200]
plot(51:200,apply(errtresd,1,mean)[51:200],type="l",ylim=c(-2,4), main="nonensembled predictive error: doppler", xlab="n",ylab="average predictive L1 error +/- SD")
polygon(c(51:200,rev(51:200)),c(errtresdhi,rev(errtresdlo)),col='grey',border=T,lty=3)
lines(51:200,apply(errtresd,1,mean)[51:200])
lines(51:200,apply(errnnresd,1,mean)[51:200],lty=6)
polygon(c(51:200,rev(51:200)),c(errnnresdhi,rev(errnnresdlo)),col=rgb(0,0,0,max=255,alpha=125),border=T,lty=2)
lines(51:200,apply(errnnresd,1,mean)[51:200],lty=6)
legend('bottomright',legend=c('trees','neural net'),lty=c(1,6))
# save as noens_doppler_error_mean+sd

par(mfrow=c(2,2))
boxplot(errtrest, main="nonensembled predictive error: trees and treed",cex.main=0.6,ylim=c(0,11))
boxplot(errnnrest, main="nonensembled predictive error: neural net and treed",cex.main=0.6,ylim=c(0,11))
boxplot(errtresd, main="nonensembled predictive error: trees and doppler",cex.main=0.6,ylim=c(0,4))
boxplot(errnnresd, main="nonensembled predictive error: neural net and doppler",cex.main=0.6,ylim=c(0,4))
# save as boxplot_noens_errors

load("Book/ch11/sim_bagging_res.RData")
R<-30
errtrest<-matrix(999,200,30)
for(i in 1:R) errtrest[,i]<-abs(trest[[i]][,2]-trest[[i]][,5])
sdtrest<-apply(errtrest,1,function(x){sqrt(var(x))})
errnnrest<-matrix(999,200,30)
for(i in 1:R) errnnrest[,i]<-abs(nnrest[[i]][,2]-nnrest[[i]][,5])
sdnnrest<-apply(errnnrest,1,function(x){sqrt(var(x))}) 
errtresthi<-apply(errtrest,1,mean)[51:200]+sdtrest[51:200]
errtrestlo<-apply(errtrest,1,mean)[51:200]-sdtrest[51:200]
errnnresthi<-apply(errnnrest,1,mean)[51:200]+sdnnrest[51:200]
errnnrestlo<-apply(errnnrest,1,mean)[51:200]-sdnnrest[51:200]
plot(51:200,apply(errtrest,1,mean)[51:200],type="l",ylim=c(-2,4), main="bagging predictive error: treed", xlab="n",ylab="average predictive L1 error +/- SD")
polygon(c(51:200,rev(51:200)),c(errtresthi,rev(errtrestlo)),col='grey',border=T,lty=3)
lines(51:200,apply(errtrest,1,mean)[51:200])
lines(51:200,apply(errnnrest,1,mean)[51:200],lty=6)
polygon(c(51:200,rev(51:200)),c(errnnresthi,rev(errnnrestlo)),col=rgb(0,0,0,max=255,alpha=125),border=T,lty=2)
lines(51:200,apply(errnnrest,1,mean)[51:200],lty=6)
 legend('bottomright',legend=c('trees','neural net'),lty=c(1,6))
# save as bagging_treed_error_mean+sd

errtresd<-matrix(999,200,30)
for(i in 1:R) errtresd[,i]<-abs(tresd[[i]][,2]-tresd[[i]][,5])
sdtresd<-apply(errtresd,1,function(x){sqrt(var(x))})
errnnresd<-matrix(999,200,30)
for(i in 1:R) errnnresd[,i]<-abs(nnresd[[i]][,2]-nnresd[[i]][,5])
sdnnresd<-apply(errnnresd,1,function(x){sqrt(var(x))})
errtresdhi<-apply(errtresd,1,mean)[51:200]+sdtresd[51:200]
errtresdlo<-apply(errtresd,1,mean)[51:200]-sdtresd[51:200]
errnnresdhi<-apply(errnnresd,1,mean)[51:200]+sdnnresd[51:200]
errnnresdlo<-apply(errnnresd,1,mean)[51:200]-sdnnresd[51:200]
plot(51:200,apply(errtresd,1,mean)[51:200],type="l",ylim=c(-2,4), main="bagging predictive error: doppler", xlab="n",ylab="average predictive L1 error +/- SD")
polygon(c(51:200,rev(51:200)),c(errtresdhi,rev(errtresdlo)),col='grey',border=T,lty=3)
lines(51:200,apply(errtresd,1,mean)[51:200])
lines(51:200,apply(errnnresd,1,mean)[51:200],lty=6)
polygon(c(51:200,rev(51:200)),c(errnnresdhi,rev(errnnresdlo)),col=rgb(0,0,0,max=255,alpha=125),border=T,lty=2)
lines(51:200,apply(errnnresd,1,mean)[51:200],lty=6)
legend('bottomright',legend=c('trees','neural net'),lty=c(1,6))
# save as noens_doppler_error_mean+sd

par(mfrow=c(2,2))
boxplot(errtrest, main="bagging predictive error: trees and treed",cex.main=0.6,ylim=c(0,13))
boxplot(errnnrest, main="bagging predictive error: neural net and treed",cex.main=0.6,ylim=c(0,13))
boxplot(errtresd, main="bagging predictive error: trees and doppler",cex.main=0.6,ylim=c(0,4))
boxplot(errnnresd, main="bagging predictive error: neural net and doppler",cex.main=0.6,ylim=c(0,4))
# save as boxplot_bagging_errors

load("Book/ch11/sim_stacking_res.RData")
R<-30
errtrest<-matrix(999,200,30)
for(i in 1:R) errtrest[,i]<-abs(trest[[i]][,2]-trest[[i]][,6])
sdtrest<-apply(errtrest,1,function(x){sqrt(var(x))})
errnnrest<-matrix(999,200,30)
for(i in 1:R) errnnrest[,i]<-abs(nnrest[[i]][,2]-nnrest[[i]][,6])
sdnnrest<-apply(errnnrest,1,function(x){sqrt(var(x))}) 
errtresthi<-apply(errtrest,1,mean)[51:200]+sdtrest[51:200]
errtrestlo<-apply(errtrest,1,mean)[51:200]-sdtrest[51:200]
errnnresthi<-apply(errnnrest,1,mean)[51:200]+sdnnrest[51:200]
errnnrestlo<-apply(errnnrest,1,mean)[51:200]-sdnnrest[51:200]
plot(51:200,apply(errtrest,1,mean)[51:200],type="l",ylim=c(-2,4), main="stacking predictive error: treed", xlab="n",ylab="average predictive L1 error +/- SD")
polygon(c(51:200,rev(51:200)),c(errtresthi,rev(errtrestlo)),col='grey',border=T,lty=3)
lines(51:200,apply(errtrest,1,mean)[51:200])
lines(51:200,apply(errnnrest,1,mean)[51:200],lty=6)
polygon(c(51:200,rev(51:200)),c(errnnresthi,rev(errnnrestlo)),col=rgb(0,0,0,max=255,alpha=125),border=T,lty=2)
lines(51:200,apply(errnnrest,1,mean)[51:200],lty=6)
 legend('bottomright',legend=c('trees','neural net'),lty=c(1,6))
# save as stacking_treed_error_mean+sd

errtresd<-matrix(999,200,30)
for(i in 1:R) errtresd[,i]<-abs(tresd[[i]][,2]-tresd[[i]][,6])
sdtresd<-apply(errtresd,1,function(x){sqrt(var(x))})
errnnresd<-matrix(999,200,30)
for(i in 1:R) errnnresd[,i]<-abs(nnresd[[i]][,2]-nnresd[[i]][,6])
sdnnresd<-apply(errnnresd,1,function(x){sqrt(var(x))})
errtresdhi<-apply(errtresd,1,mean)[51:200]+sdtresd[51:200]
errtresdlo<-apply(errtresd,1,mean)[51:200]-sdtresd[51:200]
errnnresdhi<-apply(errnnresd,1,mean)[51:200]+sdnnresd[51:200]
errnnresdlo<-apply(errnnresd,1,mean)[51:200]-sdnnresd[51:200]
plot(51:200,apply(errtresd,1,mean)[51:200],type="l",ylim=c(-2,4), main="stacking predictive error: doppler", xlab="n",ylab="average predictive L1 error +/- SD")
polygon(c(51:200,rev(51:200)),c(errtresdhi,rev(errtresdlo)),col='grey',border=T,lty=3)
lines(51:200,apply(errtresd,1,mean)[51:200])
lines(51:200,apply(errnnresd,1,mean)[51:200],lty=6)
polygon(c(51:200,rev(51:200)),c(errnnresdhi,rev(errnnresdlo)),col=rgb(0,0,0,max=255,alpha=125),border=T,lty=2)
lines(51:200,apply(errnnresd,1,mean)[51:200],lty=6)
legend('bottomright',legend=c('trees','neural net'),lty=c(1,6))
# save as stacking_doppler_error_mean+sd

par(mfrow=c(2,2))
boxplot(errtrest, main="stacking predictive error: trees and treed",cex.main=0.6,ylim=c(0,14))
boxplot(errnnrest, main="stacking predictive error: neural net and treed",cex.main=0.6,ylim=c(0,14))
boxplot(errtresd, main="stacking predictive error: trees and doppler",cex.main=0.6,ylim=c(0,9))
boxplot(errnnresd, main="stacking predictive error: neural net and doppler",cex.main=0.6,ylim=c(0,9))
# save as boxplot_stacking_errors

load("/Users/jclarke/Documents/predictive/Book/ch11/sim_boost_trees_res.RData")
load("/Users/jclarke/Documents/predictive/Book/ch11/sim_boost_nn_res.RData")
nnrest<-allnnrest
nnresd<-allnnresd
R<-30
errtrest<-matrix(999,200,30)
for(i in 1:R) errtrest[,i]<-abs(trest[[i]][,2]-trest[[i]][,7])
sdtrest<-apply(errtrest,1,function(x){sqrt(var(x))})
errnnrest<-matrix(999,200,30)
for(i in 1:R) errnnrest[,i]<-abs(nnrest[[i]][,2]-nnrest[[i]][,7])
sdnnrest<-apply(errnnrest,1,function(x){sqrt(var(x))}) 
errtresthi<-apply(errtrest,1,mean)[51:200]+sdtrest[51:200]
errtrestlo<-apply(errtrest,1,mean)[51:200]-sdtrest[51:200]
errnnresthi<-apply(errnnrest,1,mean)[51:200]+sdnnrest[51:200]
errnnrestlo<-apply(errnnrest,1,mean)[51:200]-sdnnrest[51:200]
plot(51:200,apply(errtrest,1,mean)[51:200],type="l",ylim=c(-2,4), main="boosting predictive error: treed", xlab="n",ylab="average predictive L1 error +/- SD")
polygon(c(51:200,rev(51:200)),c(errtresthi,rev(errtrestlo)),col='grey',border=T,lty=3)
lines(51:200,apply(errtrest,1,mean)[51:200])
lines(51:200,apply(errnnrest,1,mean)[51:200],lty=6)
polygon(c(51:200,rev(51:200)),c(errnnresthi,rev(errnnrestlo)),col=rgb(0,0,0,max=255,alpha=125),border=T,lty=2)
lines(51:200,apply(errnnrest,1,mean)[51:200],lty=6)
 legend('bottomright',legend=c('trees','neural net'),lty=c(1,6))
# save as boosting_treed_error_mean+sd

errtresd<-matrix(999,200,30)
for(i in 1:R) errtresd[,i]<-abs(tresd[[i]][,2]-tresd[[i]][,7])
sdtresd<-apply(errtresd,1,function(x){sqrt(var(x))})
errnnresd<-matrix(999,200,30)
for(i in 1:R) errnnresd[,i]<-abs(nnresd[[i]][,2]-nnresd[[i]][,7])
sdnnresd<-apply(errnnresd,1,function(x){sqrt(var(x))})
errtresdhi<-apply(errtresd,1,mean)[51:200]+sdtresd[51:200]
errtresdlo<-apply(errtresd,1,mean)[51:200]-sdtresd[51:200]
errnnresdhi<-apply(errnnresd,1,mean)[51:200]+sdnnresd[51:200]
errnnresdlo<-apply(errnnresd,1,mean)[51:200]-sdnnresd[51:200]
plot(51:200,apply(errtresd,1,mean)[51:200],type="l",ylim=c(-2,4), main="boosting predictive error: doppler", xlab="n",ylab="average predictive L1 error +/- SD")
polygon(c(51:200,rev(51:200)),c(errtresdhi,rev(errtresdlo)),col='grey',border=T,lty=3)
lines(51:200,apply(errtresd,1,mean)[51:200])
lines(51:200,apply(errnnresd,1,mean)[51:200],lty=6)
polygon(c(51:200,rev(51:200)),c(errnnresdhi,rev(errnnresdlo)),col=rgb(0,0,0,max=255,alpha=125),border=T,lty=2)
lines(51:200,apply(errnnresd,1,mean)[51:200],lty=6)
legend('bottomright',legend=c('trees','neural net'),lty=c(1,6))
# save as boosting_doppler_error_mean+sd

par(mfrow=c(2,2))
boxplot(errtrest, main="boosting predictive error: trees and treed",cex.main=0.6,ylim=c(0,9))
boxplot(errnnrest, main="boosting predictive error: neural net and treed",cex.main=0.6,ylim=c(0,9))
boxplot(errtresd, main="boosting predictive error: trees and doppler",cex.main=0.6,ylim=c(0,6))
boxplot(errnnresd, main="boosting predictive error: neural net and doppler",cex.main=0.6,ylim=c(0,6))
# save as boxplot_boosting_errors

