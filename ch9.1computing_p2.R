#		ch9.1 computing part 2

set.seed(4050)
library(MASS)
library(leaps)
library(plyr)

nreps<-100
nobs<-250
burnin<-50
npred<-nobs-burnin
predstep<-5
nvar<-4
sd<-5
alpha.to.leave<-0.05
alpha.to.enter<-0.05

# items for storing results
predP<-predS<-predAIC<-predBIC<-predN<-matrix(0,nreps,npred)
errP<-errS<-errAIC<-errBIC<-errN<-matrix(0,nreps,npred)
modP<-vector(mode="list", length=nreps)
for(ii in 1:nreps) modP[[ii]]<-vector(mode="list", length=(npred/predstep))
modN<-modS<-modAIC<-modBIC<-modP

i<-1
for(i in 1:nreps){
	ii<-jj<-1
	data<-matrix(0,nobs,nvar)
	for(ii in 1:nvar) data[,ii]<-rnorm(nobs, mean = 0, sd = sd)
	y<-0
	# true model for MDS
	# y = logX1 + sqrt(X2) + X3 + X4^2 + log(X1)*sqrt(X2) + sqrt(X2)*X3 +X3*X4^2 + X4*log(X1) + eps
	for(jj in c(1)) y<-log(abs(data[,jj]))
	for(jj in c(2)) y<-y+sign(data[,jj])*sqrt(abs(data[,jj]))
	for(jj in c(3)) y<-y+data[,jj]
	for(jj in c(4)) y<-y+(data[,jj])^2
	for(jj in c(5)) y<-y+log(abs(data[,(jj-4)]))*sign(data[,(jj-3)])*sqrt(abs(data[,(jj-3)]))
	for(jj in c(6)) y<-y+data[,(jj-3)]*sign(data[,(jj-4)])*sqrt(abs(data[,(jj-4)]))
	for(jj in c(7)) y<-y+data[,(jj-4)]*(data[,(jj-3)])^2
	for(jj in c(8)) y<-y+sign(data[,(jj-7)])*log(abs(data[,(jj-7)]))*(data[,(jj-4)])^2
	y<-y+rnorm(nobs, mean = 0, sd = sd/2)
	datasub<-as.data.frame(cbind(y,data))
	colnames(datasub)<-c("Y",paste("X",1:nvar,sep=""))
	j<-burnin
	k<-1
	while(j<nobs){
		tmp<-sample(c(1:4))
		tmp2<-paste("X",tmp,sep="")
		fmla<-as.formula(paste("Y ~ ",paste(tmp2, collapse="+"),sep=""))
		full<-fmla
#		full<-as.formula("Y ~ X11+X12+X13")
#		full<-as.formula("Y ~ X1+X2+X3+X4")
		null<-as.formula("Y ~ 1")
		# stepwise F test
		stepF<-stepwise(dataframe=datasub[1:j,], full.model=full, initial.model=full, alpha.to.enter, alpha.to.leave)
		# predict
		predP[i,j-burnin+1:5]<-predict(stepF, newdata=datasub[j+1:5,])
		# error
		errP[i,j-burnin+1:5]<-predP[i,j-burnin+1:5]-datasub[j+1:5,1]
		# save
		tmp<-attr(terms(stepF),"term.labels")	
		if(length(tmp)>0){
			modP[[i]][[k]]<-as.numeric(sapply(tmp,function(X){unlist(strsplit(X,"X"))})[2,])
		} else { modP[[i]][[k]]<-0}	
		# best subsets adjR2
		full<-as.formula("Y ~ .")
		bestF<-regsubsets(x=full, data=datasub[1:j,], weights=NULL, nbest=1, nvmax=30, force.in=NULL, force.out=NULL, intercept=TRUE, method=c("exhaustive"), really.big=TRUE)
		bestFind<-summary(bestF)$which[which.max(summary(bestF)$adjr2),]
		bestFmod<-as.formula(paste("Y",paste(colnames(datasub)[-1][bestFind[-1]],collapse="+"),sep="~"))
		bestFlm<-lm(bestFmod,data=datasub,subset=1:j)
		# predict
		predS[i,j-burnin+1:5]<-predict(bestFlm, newdata=datasub[j+1:5,])
		# error
		errS[i,j-burnin+1:5]<-predS[i,j-burnin+1:5]-datasub[j+1:5,1]
		# save
		tmp<-attr(terms(bestFlm),"term.labels")	
		if(length(tmp)>0){
			modS[[i]][[k]]<-as.numeric(sapply(tmp,function(X){unlist(strsplit(X,"X"))})[2,])	
		} else { modS[[i]][[k]]<-0}	
		# best BIC (thanks to Geyer, Spring 2015, Stat5102 notes)
		# see http://www.stat.umn.edu/geyer/5102/examp/select.html
		sout<-summary(bestF)
		p <- apply(sout$which, 1, sum)
		bic<-sout$bic - log(j) * p + log(j) * (p/2)
		ibest <- seq(along = sout$bic)[sout$bic == min(sout$bic)]
		foo <- sout$which[ibest, ]
		form <- names(foo)[foo][-1]
		form <- paste(form, collapse = " + ")
		form <- paste("y ~", form)
		bestBIC <- lm(as.formula(form), data=datasub,subset=1:j)
		# predict
		predBIC[i,j-burnin+1:5]<-predict(bestBIC, newdata=datasub[j+1:5,])
		# error
		errBIC[i,j-burnin+1:5]<-predBIC[i,j-burnin+1:5]-datasub[j+1:5,1]
		# save
		tmp<-attr(terms(bestBIC),"term.labels")	
		if(length(tmp)>0){
			modBIC[[i]][[k]]<-as.numeric(sapply(tmp,function(X){unlist(strsplit(X,"X"))})[2,])	
		} else { modBIC[[i]][[k]]<-0}	
		# best AIC
		p <- apply(sout$which, 1, sum)
		aic <- sout$bic - log(j) * p + 2 * p
		ibest <- seq(along = aic)[aic == min(aic)]
		foo <- sout$which[ibest, ]
		form <- names(foo)[foo][-1]
		form <- paste(form, collapse = " + ")
		form <- paste("y ~", form)
		bestAIC <- lm(as.formula(form), data=datasub,subset=1:j)
		# predict
		predAIC[i,j-burnin+1:5]<-predict(bestAIC, newdata=datasub[j+1:5,])
		# error
		errAIC[i,j-burnin+1:5]<-predAIC[i,j-burnin+1:5]-datasub[j+1:5,1]
		# save		
		tmp<-attr(terms(bestAIC),"term.labels")	
		if(length(tmp)>0){
			modAIC[[i]][[k]]<-as.numeric(sapply(tmp,function(X){unlist(strsplit(X,"X"))})[2,])	
		} else { modAIC[[i]][[k]]<-0}	
		# best using compromise btwn AIC and BIC
		p <- apply(sout$which, 1, sum)
		new <- sout$bic - log(j) * p + j^(1/2) * p
		ibest <- seq(along = new)[new == min(new)]
		foo <- sout$which[ibest, ]
		form <- names(foo)[foo][-1]
		form <- paste(form, collapse = " + ")
		form <- paste("y ~", form)
		bestN <- lm(as.formula(form), data=datasub,subset=1:j)
		# predict
		predN[i,j-burnin+1:5]<-predict(bestN, newdata=datasub[j+1:5,])
		# error
		errN[i,j-burnin+1:5]<-predN[i,j-burnin+1:5]-datasub[j+1:5,1]
		# save
		tmp<-attr(terms(bestN),"term.labels")	
		if(length(tmp)>0){
			modN[[i]][[k]]<-as.numeric(sapply(tmp,function(X){unlist(strsplit(X,"X"))})[2,])		
		} else { modN[[i]][[k]]<-0}	
		# update
		j<-j+predstep
		k<-k+1
	}
}	
# CPE plots for 5 methods...
nval<-c(burnin+1:npred)
cpeP<-t(apply(errP,1,function(x){cumsum(abs(x))}))
resP<-matrix(0,2,npred)
resP[1,]<-apply(cpeP,2,median)
resP[2,]<-apply(cpeP,2,var)/nval
resP[2,]<-sqrt(resP[2,])
plot(1:npred,resP[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="stepwise")
lines(1:npred,resP[1,]+2*resP[2,],col="violet")
lines(1:npred,resP[1,]-2*resP[2,],col="violet")
cpeS<-t(apply(errS,1,function(x){cumsum(abs(x))}))
resS<-matrix(0,2,npred)
resS[1,]<-apply(cpeS,2,median)
resS[2,]<-apply(cpeS,2,var)/nval
resS[2,]<-sqrt(resS[2,])
lines(1:npred,resS[1,],lty=2)
#plot(1:npred,resS[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="basic model")
lines(1:npred,resS[1,]+2*resS[2,],col="red")
lines(1:npred,resS[1,]-2*resS[2,],col="red")
cpeAIC<-t(apply(errAIC,1,function(x){cumsum(abs(x))}))
resAIC<-matrix(0,2,npred)
resAIC[1,]<-apply(cpeAIC,2,median)
resAIC[2,]<-apply(cpeAIC,2,var)/nval
resAIC[2,]<-sqrt(resAIC[2,])
lines(1:npred,resAIC[1,],lty=2)
#plot(1:npred,resAIC[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="stepwise")
lines(1:npred,resAIC[1,]+2*resAIC[2,],col="green")
lines(1:npred,resAIC[1,]-2*resAIC[2,],col="green")
cpeBIC<-t(apply(errBIC,1,function(x){cumsum(abs(x))}))
resBIC<-matrix(0,2,npred)
resBIC[1,]<-apply(cpeBIC,2,median)
resBIC[2,]<-apply(cpeBIC,2,var)/nval
resBIC[2,]<-sqrt(resBIC[2,])
lines(1:npred,resBIC[1,],lty=2)
#plot(1:npred,resBIC[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="stepwise")
lines(1:npred,resBIC[1,]+2*resBIC[2,],col="blue")
lines(1:npred,resBIC[1,]-2*resBIC[2,],col="blue")
cpeN<-t(apply(errN,1,function(x){cumsum(abs(x))}))
resN<-matrix(0,2,npred)
resN[1,]<-apply(cpeN,2,median)
resN[2,]<-apply(cpeN,2,var)/nval
resN[2,]<-sqrt(resN[2,])
lines(1:npred,resN[1,],lty=2)
#plot(1:npred,resN[1,],type="l",xlab="predictive step",ylab="cumulative predictive error", main="stepwise")
lines(1:npred,resN[1,]+2*resN[2,],col="orange")
lines(1:npred,resN[1,]-2*resN[2,],col="orange")
#legend('topleft',legend=c('subset','AIC','BIC','n^1/2'),col=c("red","green","blue","orange"),lty=1)
legend('topleft',legend=c('stepwise','subset','AIC','BIC','n^1/2'),col=c("violet","red","green","blue","orange"),lty=1)

# CPE plots for 5 methods...only last 100
vplot<-c(100:200)
nval<-c(burnin+1:npred)
cpeP<-t(apply(errP,1,function(x){cumsum(abs(x))}))
resP<-matrix(0,2,npred)
resP[1,]<-apply(cpeP,2,median)
resP[2,]<-apply(cpeP,2,var)/nval
resP[2,]<-sqrt(resP[2,])
plot(vplot,resP[1,vplot],type="l",lty=1,xlab="predictive step",ylab="cumulative predictive error", main="stepwise")
#plot(vplot,resP[1,vplot],type="l",col="violet",xlab="predictive step",ylab="cumulative predictive error", main="stepwise")
#lines(vplot,resP[1,vplot]+2*resP[2,vplot],col="violet")
#lines(vplot,resP[1,vplot]-2*resP[2,vplot],col="violet")
cpeS<-t(apply(errS,1,function(x){cumsum(abs(x))}))
resS<-matrix(0,2,npred)
resS[1,]<-apply(cpeS,2,median)
resS[2,]<-apply(cpeS,2,var)/nval
resS[2,]<-sqrt(resS[2,])
lines(vplot,resS[1,vplot], lty=2)
#lines(vplot,resS[1,vplot], col="red")
#plot(vplot,resS[1,vplot],type="l",xlab="predictive step",ylab="cumulative predictive error", main="basic model")
#lines(vplot,resS[1,vplot]+2*resS[2,vplot],col="red")
#lines(vplot,resS[1,vplot]-2*resS[2,vplot],col="red")
cpeAIC<-t(apply(errAIC,1,function(x){cumsum(abs(x))}))
resAIC<-matrix(0,2,npred)
resAIC[1,]<-apply(cpeAIC,2,median)
resAIC[2,]<-apply(cpeAIC,2,var)/nval
resAIC[2,]<-sqrt(resAIC[2,])
lines(vplot,resAIC[1,vplot],lty=3)
#lines(vplot,resAIC[1,vplot],col="green")
#plot(vplot,resAIC[1,vplot],type="l",xlab="predictive step",ylab="cumulative predictive error", main="stepwise")
#lines(vplot,resAIC[1,vplot]+2*resAIC[2,vplot],col="green")
#lines(vplot,resAIC[1,vplot]-2*resAIC[2,vplot],col="green")
cpeBIC<-t(apply(errBIC,1,function(x){cumsum(abs(x))}))
resBIC<-matrix(0,2,npred)
resBIC[1,]<-apply(cpeBIC,2,median)
resBIC[2,]<-apply(cpeBIC,2,var)/nval
resBIC[2,]<-sqrt(resBIC[2,])
lines(vplot,resBIC[1,vplot],lty=4)
#lines(vplot,resBIC[1,vplot],col="blue")
#plot(vplot,resBIC[1,vplot],type="l",xlab="predictive step",ylab="cumulative predictive error", main="stepwise")
#lines(vplot,resBIC[1,vplot]+2*resBIC[2,vplot],col="blue")
#lines(vplot,resBIC[1,vplot]-2*resBIC[2,vplot],col="blue")
cpeN<-t(apply(errN,1,function(x){cumsum(abs(x))}))
resN<-matrix(0,2,npred)
resN[1,]<-apply(cpeN,2,median)
resN[2,]<-apply(cpeN,2,var)/nval
resN[2,]<-sqrt(resN[2,])
lines(vplot,resN[1,vplot],lty=5)
#lines(vplot,resN[1,vplot],col="black")
#plot(vplot,resN[1,vplot],type="l",xlab="predictive step",ylab="cumulative predictive error", main="stepwise")
#lines(vplot,resN[1,vplot]+2*resN[2,vplot],col="orange")
#lines(vplot,resN[1,vplot]-2*resN[2,vplot],col="orange")
#legend('topleft',legend=c('subset','AIC','BIC','n^1/2'),col=c("red","green","blue","orange"),lty=1)
legend('topleft',legend=c('stepwise','subset','AIC','BIC','n^1/2'),lty=1:5)
#legend('topleft',legend=c('stepwise','subset','AIC','BIC','n^1/2'),col=c("violet","red","green","blue","orange"),lty=1:5)

# histogram of terms selected in final models
# account for intercept only model
final<-length(modP[[1]])
termStep<-NULL
for(k in 1:nreps){
	termStep<-c(termStep,modP[[k]][[final]]+1)
}
termSub<-NULL
for(k in 1:nreps){
	termSub<-c(termSub,modS[[k]][[final]]+1)
}
termAIC<-NULL
for(k in 1:nreps){
	termAIC<-c(termAIC,modAIC[[k]][[final]]+1)
}
termBIC<-NULL
for(k in 1:nreps){
	termBIC<-c(termBIC,modBIC[[k]][[final]]+1)
}
termN<-NULL
for(k in 1:nreps){
	termN<-c(termN,modN[[k]][[final]]+1)
}
par(mfrow=c(3,2))
hist(termStep, breaks=seq(0.5,nvar+1.5,by=1), xlim=c(1,nvar+1),ylim=c(0,100))
hist(termSub, breaks=seq(0.5,nvar+1.5,by=1),xlim=c(1,nvar+1),ylim=c(0,100))
hist(termAIC, breaks=seq(0.5,nvar+1.5,by=1),xlim=c(1,nvar+1),ylim=c(0,100))
hist(termBIC, breaks=seq(0.5,nvar+1.5,by=1),xlim=c(1,nvar+1),ylim=c(0,100))
hist(termN, breaks=seq(0.5,nvar+1.5,by=1),xlim=c(1,nvar+1),ylim=c(0,100))

# MDS of models selected
final<-length(modP[[1]])
modelsP<-matrix(0,nreps,nvar+1)
for(k in 1:nreps){
	modelsP[k,c(modP[[k]][[final]]+1)]<-1
}
modelsS<-matrix(0,nreps,nvar+1)
for(k in 1:nreps){
	modelsS[k,c(modS[[k]][[final]]+1)]<-1
}
modelsAIC<-matrix(0,nreps,nvar+1)
for(k in 1:nreps){
	modelsAIC[k,c(modAIC[[k]][[final]]+1)]<-1
}
modelsBIC<-matrix(0,nreps,nvar+1)
for(k in 1:nreps){
	modelsBIC[k,c(modBIC[[k]][[final]]+1)]<-1
}
modelsN<-matrix(0,nreps,nvar+1)
for(k in 1:nreps){
	modelsN[k,c(modN[[k]][[final]]+1)]<-1
}
# make size proportional to frequency
#tmp<-as.data.frame(modelsN)
#tmp<-as.data.frame(modelsBIC)
#tmp<-as.data.frame(modelsAIC)
#tmp<-as.data.frame(modelsS)
tmp<-as.data.frame(modelsP)
tmp2<-count(tmp, vars=c("V1","V2","V3","V4","V5"))
cnt<-tmp2$freq
size<-cnt/sum(cnt)
tmp3<-as.matrix(tmp2[,1:(nvar+1)])
d<-dist(tmp3)
fit<-cmdscale(d,eig=TRUE,k=2)
# 
x<-fit$points[,1]
y<-fit$points[,2]
plot(x,y,xlab="Coordinate1", ylab="Coordinate2",main="MDS representation of sampling distribution of stepwise")
colnames(tmp3)<-c("intercept","X1","X2","X3","X4")
for(i in 1:dim(fit$points)[1]){
	points(x=fit$points[i,1],y=fit$points[i,2],cex=cnt[i])
	if((fit$points[i,2]>=0)&(fit$points[i,1]>0)){
	 text(x=fit$points[i,1]-0.05,y=fit$points[i,2]-0.05,paste(colnames(tmp3)[which(tmp3[i,]>0)],collapse="+"),cex=0.7)
	}
	if((fit$points[i,2]>=0)&(fit$points[i,1]<=0)){
	 text(x=fit$points[i,1]+0.05,y=fit$points[i,2]-0.05,paste(colnames(tmp3)[which(tmp3[i,]>0)],collapse="+"),cex=0.7)
	}	
	if((fit$points[i,2]<0)&(fit$points[i,1]<=0)){
	 text(x=fit$points[i,1]+0.05,y=fit$points[i,2]+0.05,paste(colnames(tmp3)[which(tmp3[i,]>0)],collapse="+"),cex=0.7)
	}
	if((fit$points[i,2]<0)&(fit$points[i,1]>0)){
	 text(x=fit$points[i,1]-0.05,y=fit$points[i,2]+0.05,paste(colnames(tmp3)[which(tmp3[i,]>0)],collapse="+"),cex=0.7)
	}
}

# courtesy of Paul Rubin 
# see http://orinanobworld.blogspot.com/2011/02/stepwise-regression-in-r.html
# This is an R function to perform stepwise regression based on a "nested model" F test for inclusion/exclusion
    # of a predictor.  To keep it simple, I made no provision for forcing certain variables to be included in
    # all models, did not allow for specification of a data frame, and skipped some consistency checks (such as whether
    # the initial model is a subset of the full model).
    #
    # One other note: since the code uses R's drop1 and add1 functions, it respects hierarchy in models. That is,
    # regardless of p values, it will not attempt to drop a term while retaining a higher order interaction
    # involving that term, nor will it add an interaction term if the lower order components are not all present.
    # (You can of course defeat this by putting interactions into new variables and feeding it what looks like
    # a first-order model.)
    #
    # Consider this to be "beta" code (and feel free to improve it).  I've done very limited testing on it.
    #
    # Author: Paul A. Rubin (rubin@msu.edu)
    #
    # adapted original stepwise code for only backward, forward, data frame, and return chosen
    
stepwise <- function(dataframe, full.model, initial.model, alpha.to.enter, alpha.to.leave) {
     # full.model is the model containing all possible terms
      # initial.model is the first model to consider
      # alpha.to.enter is the significance level above which a variable may enter the model
      # alpha.to.leave is the significance level below which a variable may be deleted from the model
      # (Useful things for someone to add: specification of a data frame; a list of variables that must be included)
      full <- lm(full.model, data=dataframe);  # fit the full model
      msef <- (summary(full)$sigma)^2;  # MSE of full model
      n <- length(full$residuals);  # sample size
      allvars <- attr(full$terms, "predvars");  # this gets a list of all predictor variables
      current <- lm(initial.model, data=dataframe);  # this is the current model
      while (TRUE) {  # process each model until we break out of the loop
        temp <- summary(current);  # summary output for the current model
        rnames <- rownames(temp$coefficients);  # list of terms in the current model
#        print(temp$coefficients);  # write the model description
        p <- dim(temp$coefficients)[1];  # current model's size
        mse <- (temp$sigma)^2;  # MSE for current model
        cp <- (n-p)*mse/msef - (n-2*p);  # Mallow's cp
#        fit <- sprintf("\nS = %f, R-sq = %f, R-sq(adj) = %f, C-p = %f",
#                       temp$sigma, temp$r.squared, temp$adj.r.squared, cp);
#        write(fit, file="");  # show the fit
#        write("=====", file="");  # print a separator
        if (p > 1) {  # don't try to drop a term if only one is left
          d <- drop1(current, test="F");  # looks for significance of terms based on F tests
          pmax <- max(d[-1,6]);  # maximum p-value of any term (have to skip the intercept to avoid an NA value)
          if (pmax > alpha.to.leave) {
            # we have a candidate for deletion
            var <- rownames(d)[d[,6] == pmax];  # name of variable to delete
            if (length(var) > 1) {
              # if an intercept is present, it will be the first name in the list
              # there also could be ties for worst p-value
              # taking the second entry if there is more than one is a safe solution to both issues
              var <- var[2];
            }
#            write(paste("--- Dropping", var, "\n"), file="");  # print out the variable to be dropped
            f <- formula(current);  # current formula
            f <- as.formula(paste(f[2], "~", paste(f[3], var, sep=" - ")));  # modify the formula to drop the chosen variable (by subtracting it)
            current <- lm(f, data=dataframe);  # fit the modified model
            next;  # return to the top of the loop
          }
        }
        # if we get here, we failed to drop a term; try adding one
        # note: add1 throws an error if nothing can be added (current == full), which we trap with tryCatch
        a <- tryCatch(add1(current, scope=full, test="F"), error=function(e) NULL);  # looks for significance of possible additions based on F tests
        if (is.null(a)) {
          break;  # there are no unused variables (or something went splat), so we bail out
        }
        pmin <- min(a[-1,6]);  # minimum p-value of any term (skipping the intercept again)
        if (pmin < alpha.to.enter) {
          # we have a candidate for addition to the model
          var <- rownames(a)[a[,6] == pmin];  # name of variable to add
          if (length(var) > 1) {
            # same issue with ties, intercept as above
            var <- var[2];
          }
#          write(paste("+++ Adding", var, "\n"), file="");  # print the variable being added
          f <- formula(current);  # current formula
          f <- as.formula(paste(f[2], "~", paste(f[3], var, sep=" + ")));  # modify the formula to add the chosen variable
          current <- lm(f, data=dataframe);  # fit the modified model
          next;  # return to the top of the loop
        }
        # if we get here, we failed to make any changes to the model; time to punt
                break;
      } 
        return(current);
}		


