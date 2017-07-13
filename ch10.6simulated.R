## ch10.6 computing
# simulation with doppler function
# 

library(wavethresh)
# get doppler function
tt <- sort(runif(100))
dopp <- doppler(tt)
#plot(tt, dopp, type="l")
# get legendre polynomials up to order 10
# http://www.r-bloggers.com/fitting-legendre-orthogonal-polynomials-in-r/
library(orthopolynom)
leg10coef <- legendre.polynomials(n=10, normalized=TRUE)
# To use this polynomial in a model, we need to create a design matrix with sensible 
# column names and without the intercept:
x <- tt
leg10 <- as.matrix(as.data.frame(polynomial.values(polynomials=leg10coef,x=x)))
colnames(leg10) <- paste("coef",0:10,sep="")
leg10 <- leg10[, 2:ncol(leg10)]
# values of noise
sig<-c(0.1)
dopp1<-dopp+rnorm(length(tt),0,sig)

## trees
library(rpart)
library(rpart.plot)
set.seed(1001)
dat<-as.data.frame(cbind(dopp1,leg10))
fit1<-rpart(dopp1 ~ ., data=dat, method="anova") 
rpart.plot(fit1, type=1,main="Regression: Doppler", ycompress=T, extra=1, branch=1, varlen=0, digits=3, round=0.5, shadow.col="gray", box.palette="Grays", fallen.leaves=FALSE, split.cex=1.4, split.suffix=" ?", split.box.col="lightgray",  split.border.col="darkgray", split.round=.5,yesno.yshift=0.6, boxes.include.gap=TRUE,eq=" ",lt=" < ", ge=" >= ")

fit1pred<-predict(fit1)
plot(tt, dopp, type="l", xlab="x", ylab="doppler", ylim=c(-0.7,0.7), lwd=2, cex=1.2, cex.lab=1.2, cex.axis=1.2)
lines(tt,dopp1,lty=2, lwd=2)
lines(tt,fit1pred,lty=3, lwd=2)
legend('topleft',legend=c("function","data","fitted"), lty=1:3, cex=1.2)
# save as rpart_fitted_v2.pdf

## analyses with NN
# use caret package and vary learning parameter as well as size
library(caret)
set.seed(2020)
mygrid<-expand.grid(.decay=c(0.5, 0.1, 5e-4), .size=c(6,8,10,12))
fit1nn.caret<-train(dopp1 ~ ., data=dat, "nnet", tuneGrid=mygrid, linout=TRUE, maxit=10000, trace=F)
fit1nn.caret
#Neural Network 
#100 samples
# 10 predictor
#No pre-processing
#Resampling: Bootstrapped (25 reps) 
#Summary of sample sizes: 100, 100, 100, 100, 100, 100, ... 
#Resampling results across tuning parameters:
#  decay   size  RMSE          Rsquared      RMSE SD        Rsquared SD  
#  0.0005   6    0.3322081248  0.4658452358  0.18529412868  0.21980680598
#  0.0005   8    0.3548313925  0.4474177617  0.19179833102  0.19648652043
#  0.0005  10    0.3392132294  0.4485362017  0.11025100409  0.13411003326
#  0.0005  12    0.3632213360  0.4420866538  0.16278556681  0.20320292577
#  0.1000   6    0.2639084128  0.3678660845  0.03264790461  0.12533362590
#  0.1000   8    0.2634568654  0.3713608962  0.03210656395  0.12164322472
#  0.1000  10    0.2635622307  0.3710524734  0.03329431403  0.12834273739
#  0.1000  12    0.2641481582  0.3676404096  0.03385648015  0.13167276409
#  0.5000   6    0.2812871610  0.2796842847  0.02124841633  0.09587794023
#  0.5000   8    0.2812008889  0.2801300509  0.02124693813  0.09582232787
#  0.5000  10    0.2811459886  0.2804102401  0.02125479926  0.09583631279
#  0.5000  12    0.2811076435  0.2806198379  0.02126225800  0.09581593606
#RMSE was used to select the optimal model using  the smallest value.
#The final values used for the model were size = 8 and decay = 0.1. 
plot(fit1nn.caret, par.settings=list(superpose.line=list(col=c('black','grey25','grey50'),lwd=2),superpose.symbol=list(col='black',cex=1.5),axis.text=list(cex=1.2),par.xlab.text=list(cex=1.2),par.ylab.text=list(cex=1.2)))
# save as nnet_caret_train_v2.pdf

library(devtools)
source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')
plot(fit1nn.caret$finalModel, circle.cex=7, cex.val=1.2, circle.col="grey75")
# save as nnet_caret_best_v2.pdf

fit1nn.pred<-predict(fit1nn.caret$finalModel, newdata=dat)
plot(tt, dopp, type="l", xlab="x", ylab="doppler", ylim=c(-0.7,0.7), lwd=2, cex=1.2, cex.lab=1.2, cex.axis=1.2)
lines(tt,dopp1,lty=2, lwd=2)
lines(tt,fit1nn.pred,lty=3, lwd=2)
legend('topleft',legend=c("function","data","fitted"), lty=1:3, cex=1.2)
# save as nnet_caret_fitted_v2.pdf

## RVM
library(kernlab)
set.seed(1001)
dat<-as.data.frame(cbind(dopp1,leg10))
fit1rvm<-rvm(dopp1 ~ ., data=dat, type="regression") 
Using automatic sigma estimation (sigest) for RBF or laplace kernel 
fit1rvmpred<-predict(fit1rvm)
plot(tt, dopp, type="l", xlab="x", ylab="doppler", ylim=c(-0.7,0.7), lwd=2, cex=1.2, cex.lab=1.2, cex.axis=1.2)
lines(tt,dopp1,lty=2, lwd=2)
lines(tt,fit1rvmpred,lty=3, lwd=2)
legend('topleft',legend=c("function","data","fitted"), lty=1:3, cex=1.2)
# save as rvm_fitted_v2.pdf

# indices of relevance vectors
RVindex(fit1rvm)

## SVM
library(e1071)
set.seed(1001)
dat<-as.data.frame(cbind(dopp1,leg10))
obj <- tune(svm, dopp1 ~ .,  data = dat, ranges = list(epsilon = seq(0,1,0.1), cost = 2^(2:9)))
print(obj)
#Parameter tuning of ‘svm’:
#- sampling method: 10-fold cross validation 
#- best parameters:
# epsilon cost
#     0.8  512
#- best performance: 0.04955503066 
myColours<-colorRampPalette(brewer.pal(9,"YlOrRd"), space = "Lab")
plot(obj, color.palette=myColours)
# convert to greyscale and save as svm_train_v2.pdf

eps<-0.8
cost<-512
fit1svm<-svm(dopp1 ~ ., data=dat, type="eps-regression", epsilon=eps, cost=cost) 
fit1svmpred<-predict(fit1svm)
plot(tt, dopp, type="l", xlab="x", ylab="doppler", ylim=c(-0.7,0.7), lwd=2, cex=1.2, cex.lab=1.2, cex.axis=1.2)
lines(tt,dopp1,lty=2, lwd=2)
lines(tt,fit1svmpred,lty=3, lwd=2)
points(tt[fit1svm$index],rep(-0.6,length(fit1svm$index)), lwd=2)
legend('topleft',legend=c("function","data","fitted"), lty=1:3, cex=1.2)
# save as svm_fitted_2_v2.pdf

## SCAD
library(ncvreg)
set.seed(4001)
fit1scad<-ncvreg(X=leg10, y=dopp1, family="gaussian", penalty="SCAD")
myColours<-colorRampPalette(brewer.pal(9,"YlOrRd"), space = "Lab")
plot(fit1scad, cex.lab=1.6,cex.axis=1.2,lwd=2,col=myColours(8),shade=FALSE)
# use acrobat to save grayscale as scad_lambda_plot_v2.pdf

lambda<-fit1scad$lambda[which(fit1scad$loss==min(fit1scad$loss))[1]]
fit1scad
#$beta
#                    0.0001
#(Intercept)  0.31950198752
#coef1       -0.35059049967
#coef2        0.06997041999
#coef3        0.16449679350
#coef4        0.03098710080
#coef5       -0.17514667495
#coef6       -0.06315106350
#coef7        0.09019374949
#coef8       -0.17101581786
#coef9        0.00000000000
#coef10       0.22314066160
#$lambda
#[1] 0.0001322916358
fit1scadpred<-predict(fit1scad, X=leg10, type="response", lambda=lambda)
plot(tt, dopp, type="l", xlab="x", ylab="doppler", ylim=c(-0.7,0.7), lwd=2, cex=1.2, cex.lab=1.2, cex.axis=1.2)
lines(tt,dopp1,lty=2,lwd=2)
lines(tt,fit1scadpred,lty=3,lwd=2)
legend('topleft',legend=c("function","data","fitted"), lty=1:3, cex=1.2)
# save as scad_fitted_v2.pdf

