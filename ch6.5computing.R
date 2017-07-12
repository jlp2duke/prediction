#	computation for longitudinal chapter

#	linear model everything independent freq and bayes (baseline)
#	frequentist repeated measures ANOVA
#	linear model with compound symmetry freq and bayes
#	linear mixed models with Y, X, Year as fixed and elevation and Station random
#		freq and bayes

# input data for sites with >=25 observations from 1971-1999 inclusive
#	for single missing data was linearly interpolated
#	for >1 missing adjacent, use lognormal

dat<-read.csv('Gage_Locations_Yearly_data_update.csv', header=T, row.names=NULL)
dim(dat)
#[1] 29 33
nstat<-dim(dat)[1]
dat2<-dat
for(i in 1:dim(dat)[1]){
	nmiss<-sum(is.na(dat[i,5:33]))
	if(nmiss>0){
		mloc<-which(is.na(dat[i,]))
		mu<-mean(as.numeric(dat[i,5:33]), na.rm=T)
		sd<-sqrt(var(as.numeric(dat[i,5:33]), na.rm=T))
		#	see R page for lognormal and backsolve
		lmu<-log(mu^2/sqrt(sd+mu^2))
		lvar<-sqrt(log(1+(sd/mu^2)))
		dat2[i,mloc]<-rlnorm(nmiss,meanlog=lmu, sdlog=sqrt(lvar))
	}
}
plot(c(1971:1999),dat2[1,5:33],type="n",ylim=c(150,1500), xlab="year",ylab="Precipitation",main="Precipitation by Station")
for(i in 1:29) lines(c(1971:1999),dat2[i,5:33],type="l",lty=2)
lines(lowess(dat3$Year,dat3$Precip),lwd=3)

# convert to long format
dat3<-reshape(dat2, direction="long", varying=list(names(dat)[5:33]), v.names="Precip", idvar="Station", timevar="Year", times=1971:1999)
# plot by station and year
library(lattice)
xyplot(Precip ~ Year | Station, data=dat3, type = c("g","p","r"), xlab = list(label="Year",cex=1.4), ylab = list(label="Annual Precipitation",cex=1.4), aspect = "xy",par.strip.text = list(cex = 0.75),scales=list(cex=1.1))
# save as gage_xyplot_v2.pdf

##	linear (glm with identity link or fixed effects linear); see Sec 6.3
##		choose fixed effects (this is fixed effects longitudinal) (6.37)
##		measurements over time; use t=1,...,15 as training, predict t=16 (cross-sectional)
##		then predict t=17, 18, ...
library(MASS)
dat4<-dat3[dat3$Year<1986,]
gage.lm.0<-lm(Precip ~ Elevation + UTM_x + UTM_y + Year, data=dat4)
# including Station induces singularities
summary(gage.lm.0)
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-449.21 -146.81  -19.43  122.00  659.32 
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  9.769e+03  4.727e+03   2.067  0.03935 *  
#Elevation    5.591e-02  1.957e-02   2.857  0.00449 ** 
#UTM_x       -8.004e-04  3.339e-04  -2.397  0.01694 *  
#UTM_y        1.786e-05  2.339e-06   7.632  1.5e-13 ***
#Year        -4.569e+00  2.389e+00  -1.913  0.05644 .  
#Residual standard error: 215.3 on 430 degrees of freedom
#Multiple R-squared:  0.1466,	Adjusted R-squared:  0.1386 
#F-statistic: 18.46 on 4 and 430 DF,  p-value: 5.168e-14
## predict for next year (Year=1986) for all Stations
pred.lm.0<-predict(gage.lm.0, newdata=dat3[dat3$Year==1987,], type="response")
# put in loop for predictions
pred.lm.0<-vector("list",14)
cpe.lm.0<-matrix(0,29,14)
for(j in 1:14){
	dat4<-dat3[dat3$Year<(1985+j),]
	gage.lm.0<-lm(Precip ~ Elevation + UTM_x + UTM_y + Year, data=dat4)
	pred.lm.0[[j]]<-predict(gage.lm.0, newdata=dat3[dat3$Year==(1985+j),], type="response", interval="prediction")
	cpe.lm.0[,j]<-abs(dat3[dat3$Year==(1985+j),6]-pred.lm.0[[j]][,1])
}
pred.lm.0.all<-rbind(pred.lm.0[[1]],pred.lm.0[[2]],pred.lm.0[[3]],pred.lm.0[[4]],pred.lm.0[[5]],pred.lm.0[[6]],pred.lm.0[[7]],pred.lm.0[[8]],pred.lm.0[[9]],pred.lm.0[[10]],pred.lm.0[[11]],pred.lm.0[[12]],pred.lm.0[[13]],pred.lm.0[[14]])
pred.lm.0.all2<-rbind(matrix(NA,435,3),pred.lm.0.all)
dat5<-dat3
dat5$Station<-c("ARCHERS","ARDENCAPLE","KARAMA","EMBORI","ENASOIT","GATHIURU","ISIOLO","JACOBSON","KAMWAKI","LAMURIA","LOLDOTO","LOLMARIK","LORUKU","MPALA","MUGIE","MUTARA","MWEA","NANYUKI","NARO MORU","NICOLSON","OL DONYO","OL JOGI","OL MYSOR","OL PEJETA","SEGERA","SOLIO","MARMANET","TIMAU","TRENCH")
xyplot(Precip + pred.lm.0.all2[,1] + pred.lm.0.all2[,2] + pred.lm.0.all2[,3] ~ Year | Station, data=dat5, type=c("l","l","l","l"), col.line="black", xlab = list(label="Year",cex=1.4), ylab = list(label="Annual Precipitation",cex=1.4), aspect = "xy",par.strip.text = list(cex = 0.6),scales=list(cex=1.1,rot=c(45,45)))
# save as gage_lm_0_xyplot_v2.pdf

### Bayesian linear model, no random effects (?)
library(MCMCglmm)
dat4<-dat3[dat3$Year<1986,]
prior.gage.ble.0 <- list(B = list(mu = rep(0, 5), V = diag(5)*10e+08), R = list(V = 1, nu = 0.002))
gage.ble.0<-MCMCglmm(Precip ~ Elevation + UTM_x + UTM_y + Year, family="gaussian", data=dat4, prior=prior.gage.ble.0)
summary(gage.ble.0)
# Iterations = 3001:12991
# Thinning interval  = 10
# Sample size  = 1000 
# DIC: 5914.912 
# R-structure:  ~units
#      post.mean l-95% CI u-95% CI eff.samp
#units     46507    40810    52857     1000
# Location effects: Precip ~ Elevation + UTM_x + UTM_y + Year 
#             post.mean   l-95% CI   u-95% CI eff.samp  pMCMC
#(Intercept)  8.931e+03  7.735e+02  1.834e+04     1000  0.052
#Elevation    5.769e-02  2.354e-02  9.600e-02     1000  0.002
#UTM_x       -8.213e-04 -1.499e-03 -2.314e-04     1000  0.014
#UTM_y        1.780e-05  1.343e-05  2.223e-05     1000 <0.001
#Year        -4.145e+00 -8.999e+00 -1.279e-01     1000  0.066
#(Intercept) .  
#Elevation   ** 
#UTM_x       *  
#UTM_y       ***
#Year        .  
plot(gage.ble.0$Sol)
# put in loop for predictions
pred.ble.0<-vector("list",14)
cpe.ble.0<-matrix(0,29,14)
for(j in 1:14){
	x2<-as.matrix(dat3[dat3$Year<=(1985+j),2:5])
	y<-dat3[dat3$Year<=(1985+j),6]
	yNa<-y
	whichNa<-which(dat3$Year==(1985+j))
	yNa[whichNa]<-NA
	gage.ble.0<-MCMCglmm(yNa ~ Elevation + UTM_x + UTM_y + Year, family="gaussian", data=data.frame(yNa,x2), prior=prior.gage.ble.0)
	pred.ble.0[[j]]<-predict(gage.ble.0, type="response", interval="prediction")[whichNa,]
	cpe.ble.0[,j]<-abs(dat3[dat3$Year==(1985+j),6]-pred.ble.0[[j]][,1])
}
pred.ble.0.all<-rbind(pred.ble.0[[1]],pred.ble.0[[2]],pred.ble.0[[3]],pred.ble.0[[4]],pred.ble.0[[5]],pred.ble.0[[6]],pred.ble.0[[7]],pred.ble.0[[8]],pred.ble.0[[9]],pred.ble.0[[10]],pred.ble.0[[11]],pred.ble.0[[12]],pred.ble.0[[13]],pred.ble.0[[14]])
pred.ble.0.all2<-rbind(matrix(NA,435,3),pred.ble.0.all)
dat5<-dat3
dat5$Station<-c("ARCHERS","ARDENCAPLE","KARAMA","EMBORI","ENASOIT","GATHIURU","ISIOLO","JACOBSON","KAMWAKI","LAMURIA","LOLDOTO","LOLMARIK","LORUKU","MPALA","MUGIE","MUTARA","MWEA","NANYUKI","NARO MORU","NICOLSON","OL DONYO","OL JOGI","OL MYSOR","OL PEJETA","SEGERA","SOLIO","MARMANET","TIMAU","TRENCH")
xyplot(Precip + pred.ble.0.all2[,1] + pred.ble.0.all2[,2] + pred.ble.0.all2[,3] ~ Year | Station, data=dat5, type=c("l","l","l","l"), col.line="black", xlab = list(label="Year",cex=1.4), ylab = list(label="Annual Precipitation",cex=1.4), aspect = "xy",par.strip.text = list(cex = 0.6),scales=list(cex=1.1,rot=c(45,45)))
# save as gage_ble_0_xyplot_v2.pdf

# bayesian linear models with conjugate priors (R BLR package)
# for each run do CPE ...
# BLR parameters
set.seed(1080)
seed<-1080
nIter<-5500
burnIn<-500
priorC<-list(varE=list(df=3,S=1),varBR=list(df=3,S=1))
library(BLR)
pred.ble.1<-vector("list", 14)
cpe.ble.1<-matrix(0,29,14)
for(j in 1:14){
	x2<-as.matrix(dat3[dat3$Year<=(1985+j),2:5])
	y<-dat3[dat3$Year<=(1985+j),6]
	yNa<-y
	whichNa<-which(dat3$Year==(1985+j))
	yNa[whichNa]<-NA
	# suppress output to screen
	capture.output(gage.ble.1<-BLR(y=yNa,XR=x2,prior=priorC,nIter=nIter,burnIn=burnIn,saveAt=paste("example_",seed,sep="")),file=NULL)
	pred.ble.1[[j]]<-cbind(gage.ble.1$yHat[gage.ble.1$whichNa],gage.ble.1$yHat[gage.ble.1$whichNa]-1.96*gage.ble.1$SD.yHat[gage.ble.1$whichNa],gage.ble.1$yHat[gage.ble.1$whichNa]+1.96*gage.ble.1$SD.yHat[gage.ble.1$whichNa])
	cpe.ble.1[,j]<-abs(dat3[dat3$Year==(1985+j),6]-pred.ble.1[[j]][,1])
}
pred.ble.1.all<-rbind(pred.ble.1[[1]],pred.ble.1[[2]],pred.ble.1[[3]],pred.ble.1[[4]],pred.ble.1[[5]],pred.ble.1[[6]],pred.ble.1[[7]],pred.ble.1[[8]],pred.ble.1[[9]],pred.ble.1[[11]],pred.ble.1[[11]],pred.ble.1[[12]],pred.ble.1[[13]],pred.ble.1[[14]])
pred.ble.1.all2<-rbind(matrix(NA,435,3),pred.ble.1.all)
xyplot(Precip + pred.ble.1.all2[,1] + pred.ble.1.all2[,2] + pred.ble.1.all2[,3] ~ Year | Station, data=dat3, type=c("l","l","l","l"),xlab = "Year", ylab = "Annual Precipitation", aspect = "xy",par.strip.text = list(cex = 0.35))

# frequentist repeated measures aov
gage.lm.1<-aov(Precip ~ Elevation+UTM_x+UTM_y+Year+Error(Station/Year), data=dat3)
summary(gage.lm.1)
#Error: Station
#          Df   Sum Sq Mean Sq F value  Pr(>F)   
#Elevation  1   481850  481850   0.894 0.35357   
#UTM_x      1    49414   49414   0.092 0.76462   
#UTM_y      1  6266496 6266496  11.620 0.00222 **
#Residuals 25 13481753  539270                   
#Error: Station:Year
#          Df Sum Sq Mean Sq F value  Pr(>F)   
#Year       1 331138  331138   9.904 0.00389 **
#Residuals 28 936188   33435                   
#Error: Within
#           Df   Sum Sq Mean Sq F value Pr(>F)
#Residuals 783 27833595   35547               
## NO PREDICT OR FITTED FOR aov

## try AR(1) covariance structure with heterogeneous variances
library(nlme)
gage.lm.1<-gls(Precip ~ Elevation+UTM_x+UTM_y+Year, data=dat3, corr = corAR1(, form = ~ 1 | Station), weight = varIdent(form = ~ 1 | Year))
summary(gage.lm.1)
#Generalized least squares fit by REML
#  Model: Precip ~ Elevation + UTM_x + UTM_y + Year 
#  Data: dat3 
#       AIC      BIC    logLik
#  11336.97 11502.47 -5633.486
#Correlation Structure: AR(1)
# Formula: ~1 | Station 
# Parameter estimate(s):
#     Phi 
#0.474586 
#Variance function:
# Structure: Different standard deviations per stratum
# Formula: ~1 | Year 
# Parameter estimates:
#     1971      1972      1973      1974      1975      1976 
#1.0000000 0.9247595 0.9488802 0.8448122 0.6406882 1.0605252 
#     1977      1978      1979      1980      1981      1982 
#1.4915247 0.9632790 0.8680292 1.8015799 1.2938355 0.9451257 
#     1983      1984      1985      1986      1987      1988 
#0.8166830 1.6951460 0.8511889 1.1531761 1.1932078 1.7351402 
#     1989      1990      1991      1992      1993      1994 
#0.9647264 1.3181491 1.0001912 0.8362428 0.8573492 1.0702822 
#     1995      1996      1997      1998      1999 
#0.9419957 0.9866689 2.6886546 1.4356010 1.2878391 
#Coefficients:
#               Value Std.Error   t-value p-value
#(Intercept) 8449.370 2434.9961  3.469973  0.0005
#Elevation      0.072    0.0198  3.615301  0.0003
#UTM_x         -0.001    0.0003 -1.608564  0.1081
#UTM_y          0.000    0.0000  8.559134  0.0000
#Year          -3.943    1.2262 -3.215815  0.0014
# Correlation: 
#          (Intr) Elevtn UTM_x  UTM_y 
#Elevation -0.026                     
#UTM_x     -0.043  0.215              
#UTM_y     -0.005  0.128  0.001       
#Year      -0.999  0.000  0.000  0.000
#Standardized residuals:
#       Min         Q1        Med         Q3        Max 
#-3.0140588 -0.7073856 -0.0469502  0.6685699  3.3925019 
#Residual standard error: 197.5285 
#Degrees of freedom: 841 total; 836 residual
## gls does NOT give predictive intervals, only point predictors, so use different package
library(rms)
# put in loop for predictions
pred.lm.1<-vector("list", 14)
cpe.lm.1<-matrix(0,29,14)
for(j in 1:14){
	dat4<-dat3[dat3$Year<(1985+j),]
	ddist <- datadist(dat4$Elevation, dat4$UTM_x, dat4$UTM_y, dat4$Year, dat4$Station)
	options(datadist='ddist')
	gage.lm.1<-Gls(Precip ~ Elevation+UTM_x+UTM_y+Year, data=dat4, corr = corAR1(, form = ~ 1 | Station), weight = varIdent(form = ~ 1 | Year))
	pred.lm.1[[j]]<-predict(gage.lm.1, dat3[dat3$Year==(1985+j),], se.fit=T)
	pred.lm.1[[j]]<-cbind(pred.lm.1[[j]]$linear.predictor, pred.lm.1[[j]]$se.fit, pred.lm.1[[j]]$linear.predictor-1.96*pred.lm.1[[j]]$se.fit, pred.lm.1[[j]]$linear.predictor+1.96*pred.lm.1[[j]]$se.fit)
	cpe.lm.1[,j]<-abs(dat3[dat3$Year==(1985+j),6]-pred.lm.1[[j]][,1])
}
pred.lm.1.all<-rbind(pred.lm.1[[1]],pred.lm.1[[2]],pred.lm.1[[3]],pred.lm.1[[4]],pred.lm.1[[5]],pred.lm.1[[6]],pred.lm.1[[7]],pred.lm.1[[8]],pred.lm.1[[9]],pred.lm.1[[11]],pred.lm.1[[11]],pred.lm.1[[12]],pred.lm.1[[13]],pred.lm.1[[14]])
pred.lm.1.all2<-rbind(matrix(NA,435,4),pred.lm.1.all)
xyplot(Precip + pred.lm.1.all2[,1] + pred.lm.1.all2[,3] + pred.lm.1.all2[,4] ~ Year | Station, data=dat3, type=c("l","l","l","l"),xlab = "Year", ylab = "Annual Precipitation", aspect = "xy",par.strip.text = list(cex = 0.35))

## frequentist
## what if we use gls and specify compound symmetry?
library(nlme)
gage.lm.2<-gls(Precip ~ Elevation+UTM_x+UTM_y+Year, data=dat3, corr = corCompSymm(, form= ~ 1 | Station))
summary(gage.lm.2)
#Generalized least squares fit by REML
#  Model: Precip ~ Elevation + UTM_x + UTM_y + Year 
#  Data: dat3 
#       AIC      BIC    logLik
#  11315.13 11348.23 -5650.564
#Correlation Structure: Compound symmetry
# Formula: ~1 | Station 
# Parameter estimate(s):
#     Rho 
#0.332931 
#Coefficients:
#                Value Std.Error   t-value p-value
#(Intercept) -4184.934 1565.6709 -2.672933  0.0077
#Elevation       0.061    0.0484  1.248872  0.2121
#UTM_x           0.000    0.0008 -0.270404  0.7869
#UTM_y           0.000    0.0000  3.353975  0.0008
#Year            2.400    0.7763  3.091835  0.0021
# Correlation: 
#          (Intr) Elevtn UTM_x  UTM_y 
#Elevation -0.099                     
#UTM_x     -0.164  0.215              
#UTM_y     -0.019  0.128  0.001       
#Year      -0.984  0.000  0.000  0.000
#Standardized residuals:
#       Min         Q1        Med         Q3        Max 
#-2.3347025 -0.6877943 -0.1127778  0.5922008  3.6623281 
#Residual standard error: 230.6035 
#Degrees of freedom: 841 total; 836 residual
## gls does NOT give predictive intervals, only point predictors, so use different package
library(rms)
library(nlme)
# put in loop for predictions
pred.lm.2<-vector("list", 14)
cpe.lm.2<-matrix(0,29,14)
for(j in 1:14){
	dat4<-dat3[dat3$Year<(1985+j),]
	ddist <- datadist(dat4$Elevation, dat4$UTM_x, dat4$UTM_y, dat4$Year, dat4$Station)
	options(datadist='ddist')
	gage.lm.2<-Gls(Precip ~ Elevation+UTM_x+UTM_y+Year, data=dat4, corr = corCompSymm(, form= ~ 1 | Station))
	pred.lm.2[[j]]<-predict(gage.lm.2, dat3[dat3$Year==(1985+j),], se.fit=T)
	pred.lm.2[[j]]<-cbind(pred.lm.2[[j]]$linear.predictor, pred.lm.2[[j]]$se.fit, pred.lm.2[[j]]$linear.predictor-1.96*pred.lm.2[[j]]$se.fit, pred.lm.2[[j]]$linear.predictor+1.96*pred.lm.2[[j]]$se.fit)
	cpe.lm.2[,j]<-abs(dat3[dat3$Year==(1985+j),6]-pred.lm.2[[j]][,1])
}
pred.lm.2.all<-rbind(pred.lm.2[[1]],pred.lm.2[[2]],pred.lm.2[[3]],pred.lm.2[[4]],pred.lm.2[[5]],pred.lm.2[[6]],pred.lm.2[[7]],pred.lm.2[[8]],pred.lm.2[[9]],pred.lm.2[[11]],pred.lm.2[[11]],pred.lm.2[[12]],pred.lm.2[[13]],pred.lm.2[[14]])
pred.lm.2.all2<-rbind(matrix(NA,435,4),pred.lm.2.all)
dat5<-dat3
dat5$Station<-c("ARCHERS","ARDENCAPLE","KARAMA","EMBORI","ENASOIT","GATHIURU","ISIOLO","JACOBSON","KAMWAKI","LAMURIA","LOLDOTO","LOLMARIK","LORUKU","MPALA","MUGIE","MUTARA","MWEA","NANYUKI","NARO MORU","NICOLSON","OL DONYO","OL JOGI","OL MYSOR","OL PEJETA","SEGERA","SOLIO","MARMANET","TIMAU","TRENCH")
xyplot(Precip + pred.lm.2.all2[,1] + pred.lm.2.all2[,3] + pred.lm.2.all2[,4] ~ Year | Station, data=dat5, type=c("l","l","l","l"),col.line="black", xlab = list(label="Year",cex=1.4), ylab = list(label="Annual Precipitation",cex=1.4), aspect = "xy",par.strip.text = list(cex = 0.6),scales=list(cex=1.1,rot=c(45,45)))
# save as gage_lm_2_xyplot_v2.pdf

# Bayes
# compound symmetry (NOT CORRECT - this is compound symmetry on fixed effects not random effects)
# small DIC is good!
set.seed(10101)
library(MCMCglmm)
PBV.yfixed <- diag(5) * 10000 + matrix(1,nrow=5,ncol=5)*10000/2 - diag(5)
prior.gage.ble.2 <- list(B = list(mu = rep(0, 5),V = PBV.yfixed),R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu=0.002)))
gage.ble.2<-MCMCglmm(Precip ~ Elevation+UTM_x+UTM_y+Year, random=~Station, family="gaussian", data=dat3, prior=prior.gage.ble.2, pr = TRUE, saveX = TRUE, saveZ = TRUE)
# DIC: 11232.78 
# G-structure:  ~Station
#        post.mean l-95% CI u-95% CI eff.samp
#Station     18646     9654    31153     1106
# R-structure:  ~units
#      post.mean l-95% CI u-95% CI eff.samp
#units     35891    32729    39601     1165
# Location effects: Precip ~ Elevation + UTM_x + UTM_y + Year 
#             post.mean   l-95% CI   u-95% CI eff.samp pMCMC
#(Intercept) -1.981e+01 -2.214e+02  1.833e+02     1199 0.834
#Elevation    4.787e-02 -5.118e-02  1.479e-01     1000 0.338
#UTM_x       -6.282e-04 -2.324e-03  1.060e-03     1000 0.446
#UTM_y        1.915e-05  8.290e-06  3.049e-05     1000 0.002
#Year         3.728e-01  7.221e-02  6.751e-01     1000 0.012             
#(Intercept)   
#Elevation     
#UTM_x         
#UTM_y       **
#Year        * 
# fitted values after marginalizing over random effects
gage.ble.2.fitted<-predict(gage.ble.2, marginal = ~Station, type = "response", interval="prediction")
pop.int<-posterior.mode(gage.ble.2$Sol[,1])
pop.year<-posterior.mode(gage.ble.2$Sol[,2])
pop.UTM_x<-posterior.mode(gage.ble.2$Sol[,3])
pop.UTM_y<-posterior.mode(gage.ble.2$Sol[,4])
pop.Elevation<-posterior.mode(gage.ble.2$Sol[,5])
station.int<-posterior.mode(gage.ble.2$Sol[,c(6:34)])
# use posterior.mode to generate predictions for next time point
W.1<-cBind(gage.ble.2$X, gage.ble.2$Z)
pred.ble.2<-W.1%*%posterior.mode(gage.ble.2$Sol)
# put in loop for predictions
pred.ble.2<-vector("list",14)
cpe.ble.2<-matrix(99,29,14)
for(j in 1:14){
	dat4<-dat3[dat3$Year<=(1985+j),1:5]
	y<-dat3[dat3$Year<=(1985+j),6]
	yNa<-y
	whichNa<-which(dat3$Year==(1985+j))
	yNa[whichNa]<-NA
	# suppress output to screen
	dat4<-data.frame(yNa,dat4)
	gage.ble.2<-MCMCglmm(yNa ~ Elevation+UTM_x+UTM_y+Year, random=~Station, family="gaussian", data=dat4, prior=prior.gage.ble.2, pr = TRUE, saveX = TRUE, saveZ = TRUE)
#	W.1<-cBind(gage.ble.2$X, gage.ble.2$Z)
	pred.ble.2[[j]]<-predict(gage.ble.2, marginal = ~Station, type = "response", interval="prediction")
	pred.ble.2[[j]]<-pred.ble.2[[j]][((14+j)*29+1):((15+j)*29),]
#	pred.ble.2[[j]]<-as.vector(W.1[((14+j)*29+1):((15+j)*29),]%*%posterior.mode(gage.ble.2$Sol))
	cpe.ble.2[,j]<-abs(dat3[dat3$Year==(1985+j),6]-pred.ble.2[[j]][,1])
}
pred.ble.2.all<-rbind(pred.ble.2[[1]],pred.ble.2[[2]],pred.ble.2[[3]],pred.ble.2[[4]],pred.ble.2[[5]],pred.ble.2[[6]],pred.ble.2[[7]],pred.ble.2[[8]],pred.ble.2[[9]],pred.ble.2[[11]],pred.ble.2[[11]],pred.ble.2[[12]],pred.ble.2[[13]],pred.ble.2[[14]])
pred.ble.2.all2<-rbind(matrix(NA,435,3),pred.ble.2.all)
xyplot(Precip + pred.ble.2.all2[,1] + pred.ble.2.all2[,2] + pred.ble.2.all2[,3] ~ Year | Station, data=dat3, type=c("l","l","l","l"),xlab = "Year", ylab = "Annual Precipitation", aspect = "xy",par.strip.text = list(cex = 0.35))

library(nlme)
## FYI nlme vs. lme4
#As chl pointed out, the main difference is what kind of variance-covariance structure you 
#can specify for the random effects. In lme4 you can specify either:
#
#    diagonal covariance structures (i.e., enforce mutually uncorrelated random effects via 
#    syntax like ~ (1 | group)+ (0 + x1 | group) + (0 + x2 | group))
#    or unstructured covariance matrices 
#    (i.e. all correlations are estimated, ~ (1 + x1 + x2 | group))
#    or partially diagonal, partially unstructured covariance 
#    (y ~ (1 + x1 | group) + (0 + x2 | group)), where you would estimate a correlation 
#    between the random intercept and random slope for x1, but no correlations between 
#    the random slope for x2 and the random intercept and between the random slope for x2 
#    and the random slope for x1).
#	nlme offers a much broader class of covariance structures for the random effects. 
#	My experience is that the flexibility of lme4 is sufficient for most applications, however.

#	Can use gls to do repeated measures and specify structured covariance matrix
#	I think can do this using lme as well (NOT lme4)

#	linear mixed models (section 4)
#		Y, fixed effects, random effects (U_i)
#		pool data to estimate \beta and V_i by REML
#		Proc Mix? see p. 201 (6.60)
#		R nlme

library(lme4)

#Note that none of the following approaches takes the uncertainty
#of the random effects parameters into account
# try first model, random intercept+slope, effects may be correlated
gage.lme.1 <- lmer(Precip ~ Elevation+UTM_x+UTM_y+Year+(1|Station), dat3)
summary(gage.lme.1)
#Linear mixed model fit by REML ['lmerMod']
#Formula: Precip ~ Elevation + UTM_x + UTM_y + Year + (1 | Station) 
#   Data: dat3 
#REML criterion at convergence: 11301.13 
#Random effects:
# Groups   Name        Variance Std.Dev.
# Station  (Intercept) 17705    133.1   
# Residual             35473    188.3   
#Number of obs: 841, groups: Station, 29
#Fixed effects:
#              Estimate Std. Error t value
#(Intercept) -4.185e+03  1.566e+03  -2.673
#Elevation    6.050e-02  4.845e-02   1.249
#UTM_x       -2.235e-04  8.264e-04  -0.270
#UTM_y        1.942e-05  5.791e-06   3.354
#Year         2.400e+00  7.763e-01   3.092
#Correlation of Fixed Effects:
#	(Intr) Elevtn UTM_x  UTM_y 
#Elevation -0.099                     
#UTM_x     -0.164  0.215              
#UTM_y     -0.019  0.128  0.001       
#Year      -0.984  0.000  0.000  0.000
# put in loop for predictions
pred.lme.1<-vector("list",14)
for(j in 1:14) pred.lme.1[[j]]<-matrix(99,29,3)
cpe.lme.1<-matrix(99,29,14)
for(j in 1:14){
	dat4<-dat3[dat3$Year<(1985+j),]
	gage.lme.1<-lmer(Precip ~ Elevation+UTM_x+UTM_y+Year+(1|Station), data=dat4)
	mm<-as.matrix(cbind(1,dat3[dat3$Year==(1985+j),2:5]))
	pred.lme.1[[j]][,1]<-(mm%*%fixef(gage.lme.1) + ranef(gage.lme.1)$Station)$"(Intercept)"
	pvar1 <- diag(mm %*% tcrossprod(vcov(gage.lme.1),mm))	# fixed effects only
	tvar1 <- pvar1+VarCorr(gage.lme.1)$Station[1]	# FE plus RE variance
	pred.lme.1[[j]][,2]<-pred.lme.1[[j]][,1]-1.96*sqrt(tvar1)
	pred.lme.1[[j]][,3]<-pred.lme.1[[j]][,1]+1.96*sqrt(tvar1)
	cpe.lme.1[,j]<-abs(dat3[dat3$Year==(1985+j),6]-pred.lme.1[[j]][,1])
}
pred.lme.1.all<-rbind(pred.lme.1[[1]],pred.lme.1[[2]],pred.lme.1[[3]],pred.lme.1[[4]],pred.lme.1[[5]],pred.lme.1[[6]],pred.lme.1[[7]],pred.lme.1[[8]],pred.lme.1[[9]],pred.lme.1[[11]],pred.lme.1[[11]],pred.lme.1[[12]],pred.lme.1[[13]],pred.lme.1[[14]])
pred.lme.1.all2<-rbind(matrix(NA,435,3),pred.lme.1.all)
dat5<-dat3
dat5$Station<-c("ARCHERS","ARDENCAPLE","KARAMA","EMBORI","ENASOIT","GATHIURU","ISIOLO","JACOBSON","KAMWAKI","LAMURIA","LOLDOTO","LOLMARIK","LORUKU","MPALA","MUGIE","MUTARA","MWEA","NANYUKI","NARO MORU","NICOLSON","OL DONYO","OL JOGI","OL MYSOR","OL PEJETA","SEGERA","SOLIO","MARMANET","TIMAU","TRENCH")
xyplot(Precip + pred.lme.1.all2[,1] + pred.lme.1.all2[,2] + pred.lme.1.all2[,3] ~ Year | Station, data=dat5, type=c("l","l","l","l"),col.line="black", xlab = list(label="Year",cex=1.4), ylab = list(label="Annual Precipitation",cex=1.4), aspect = "xy",par.strip.text = list(cex = 0.6),scales=list(cex=1.1,rot=c(45,45)))

##	Bayesian LMM(section 2): 
##	
library(INLA)
## include random effect/repeated measures
form.inla <- Precip ~ Elevation + UTM_x + UTM_y + Year + f(Station, model="iid")
gage.inla.1 <- inla(form.inla, data=dat3, family="gaussian") 
summary(gage.inla.1)
#Call:
#"inla(formula = form.inla, family = \"gaussian\", data = dat3)"
#Time used:
#Pre-processing    Running inla Post-processing 
#         0.6033          1.4859          0.0556 
#          Total 
#         2.1448 
#Fixed effects:
#                  mean        sd 0.025quant   0.5quant
#(Intercept) -3801.2788 2032.3959 -7791.5574 -3801.3360
#Elevation       0.0618    0.0168     0.0288     0.0618
#UTM_x          -0.0002    0.0003    -0.0008    -0.0002
#UTM_y           0.0000    0.0000     0.0000     0.0000
#Year            2.2067    1.0218     0.2006     2.2067
#            0.975quant       mode kld
#(Intercept)   185.6696 -3801.2788   0
#Elevation       0.0948     0.0618   0
#UTM_x           0.0003    -0.0002   0
#UTM_y           0.0000     0.0000   0
#Year            4.2111     2.2067   0
#Random effects:
#Name	  Model
# Station   IID model 
#Model hyperparameters:
#                                        mean   sd    
#Precision for the Gaussian observations  0.000  0.000
#Precision for Station                   129.05  20.76
#                                        0.025quant 0.5quant
#Precision for the Gaussian observations  0.000      0.000  
#Precision for Station                    95.42     126.48  
#                                        0.975quant mode  
#Precision for the Gaussian observations  0.000      0.000
#Precision for Station                   176.41     120.93
#Expected number of effective parameters(std dev): 4.746(2e-04)
#Number of equivalent replicates : 177.19 
#Marginal Likelihood:  -5868.30 
# put in loop for predictions
pred.inla.1<-vector("list", 14)
cpe.inla.1<-matrix(99,29,14)
for(j in 1:14){
	x2<-as.matrix(dat3[dat3$Year<=(1985+j),2:5])
	r<-as.matrix(dat3[dat3$Year<=(1985+j),1])
	y<-dat3[dat3$Year<=(1985+j),6]
	yNa<-y
	whichNa<-which(dat3$Year==(1985+j))
	yNa[whichNa]<-NA
	form.inla <- yNa ~ x2 + f(r, model="iid")
	gage.inla.1 <- inla(form.inla, data=data.frame(yNa,x2,r), family="gaussian", control.predictor=list(compute=TRUE)) 
	pred.inla.1[[j]]<-gage.inla.1$summary.fitted.values[whichNa,]
	cpe.inla.1[,j]<-abs(dat3[dat3$Year==(1985+j),6]-pred.inla.1[[j]][,1])
}
pred.inla.1.all<-rbind(pred.inla.1[[1]],pred.inla.1[[2]],pred.inla.1[[3]],pred.inla.1[[4]],pred.inla.1[[5]],pred.inla.1[[6]],pred.inla.1[[7]],pred.inla.1[[8]],pred.inla.1[[9]],pred.inla.1[[11]],pred.inla.1[[11]],pred.inla.1[[12]],pred.inla.1[[13]],pred.inla.1[[14]])
pred.inla.1.all2<-rbind(matrix(NA,435,3),cbind(pred.inla.1.all$mean,pred.inla.1.all$"0.025quant",pred.inla.1.all$"0.975quant"))
xyplot(Precip + pred.inla.1.all2[,1] + pred.inla.1.all2[,2] + pred.inla.1.all2[,3] ~ Year | Station, data=dat3, type=c("l","l","l","l"),xlab = "Year", ylab = "Annual Precipitation", aspect = "xy",par.strip.text = list(cex = 0.35))

## same model but with MCMCglmm; inverse Wishart prior for the (co)variances and a normal prior
for fixed effects
library(MCMCglmm)
# specify prior for variance, prior for mean, prior for random effects
prior<-list(R = list(V = 1, nu = 0.002), G = list(G1=list(V = 1000, nu=0.002)))
gage.mcmc.1<-MCMCglmm(Precip ~ Elevation + UTM_x + UTM_y + Year, random = ~Station, family="gaussian", data=dat3, prior=prior, verbose=FALSE, pr=TRUE, nitt=100000, burnin=20000)
plot(gage.mcmc.1$VCV)
summary(gage.mcmc.1)
# Iterations = 20001:99991
# Thinning interval  = 10
# Sample size  = 8000 
# DIC: 11228.12 
# G-structure:  ~Station
#        post.mean l-95% CI u-95% CI eff.samp
#Station     18803     8805    30916     8000
# R-structure:  ~units
#      post.mean l-95% CI u-95% CI eff.samp
#units     35593    32207    39108     8000
# Location effects: Precip ~ Elevation + UTM_x + UTM_y + Year 
#             post.mean   l-95% CI   u-95% CI eff.samp
#(Intercept) -4.111e+03 -7.048e+03 -9.206e+02     8708
#Elevation    6.153e-02 -3.317e-02  1.654e-01     8000
#UTM_x       -2.318e-04 -1.956e-03  1.447e-03     8000
#UTM_y        1.956e-05  7.373e-06  3.124e-05     8000
#Year         2.362e+00  8.946e-01  3.891e+00     8752
#              pMCMC   
#(Intercept) 0.00675 **
#Elevation   0.20800   
#UTM_x       0.79500   
#UTM_y       0.00125 **
#Year        0.00275 **
pred.mcmc.1<-predict(gage.mcmc.1, marginal = ~Station, type="response", interval="prediction")
# put in loop for predictions
pred.mcmc.1<-vector("list", 14)
cpe.mcmc.1<-matrix(99,29,14)
for(j in 1:14){
	dat4<-dat3[dat3$Year<=(1985+j),1:5]
	y<-dat3[dat3$Year<=(1985+j),6]
	yNa<-y
	whichNa<-which(dat3$Year==(1985+j))
	yNa[whichNa]<-NA
	# suppress output to screen
	dat4<-data.frame(yNa,dat4)
	prior<-list(R = list(V = 1, nu = 0.002), B = list(mu=rep(0,5), V=diag(5)*100), G = list(G1=list(V = 100, nu=0.002)))
	gage.mcmc.1<-MCMCglmm(yNa ~ Elevation + UTM_x + UTM_y + Year, random = ~Station, family="gaussian", data=dat4, prior=prior, verbose=FALSE, pr=TRUE, nitt=100000, burnin=20000)
	pred.mcmc.1[[j]]<-predict(gage.mcmc.1, marginal = ~Station, type="response", interval="prediction")[whichNa,]
	cpe.mcmc.1[,j]<-abs(dat3[dat3$Year==(1985+j),6]-pred.mcmc.1[[j]][,1])
}
pred.mcmc.1.all<-rbind(pred.mcmc.1[[1]],pred.mcmc.1[[2]],pred.mcmc.1[[3]],pred.mcmc.1[[4]],pred.mcmc.1[[5]],pred.mcmc.1[[6]],pred.mcmc.1[[7]],pred.mcmc.1[[8]],pred.mcmc.1[[9]],pred.mcmc.1[[11]],pred.mcmc.1[[11]],pred.mcmc.1[[12]],pred.mcmc.1[[13]],pred.mcmc.1[[14]])
pred.mcmc.1.all2<-rbind(matrix(NA,435,3),pred.mcmc.1.all)
dat5<-dat3
dat5$Station<-c("ARCHERS","ARDENCAPLE","KARAMA","EMBORI","ENASOIT","GATHIURU","ISIOLO","JACOBSON","KAMWAKI","LAMURIA","LOLDOTO","LOLMARIK","LORUKU","MPALA","MUGIE","MUTARA","MWEA","NANYUKI","NARO MORU","NICOLSON","OL DONYO","OL JOGI","OL MYSOR","OL PEJETA","SEGERA","SOLIO","MARMANET","TIMAU","TRENCH")
xyplot(Precip + pred.mcmc.1.all2[,1] + pred.mcmc.1.all2[,2] + pred.mcmc.1.all2[,3] ~ Year | Station, data=dat5, type=c("l","l","l","l"),col.line="black", xlab = list(label="Year",cex=1.4), ylab = list(label="Annual Precipitation",cex=1.4), aspect = "xy",par.strip.text = list(cex = 0.6),scales=list(cex=1.1,rot=c(45,45)))
# save as gage_mcmc_1_xyplot_v2.pdf

