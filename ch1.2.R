# music21 data analysis

## classify mozart and haydn
moz<-read.csv('mozart_out.csv', header=F, row.names=1)
dim(moz)
#[1] 100 633
haydn<-read.csv('haydn_xml_out.csv', header=F, row.names=1)
dim(haydn)
#[1] 244 633
alldat<-rbind(moz,haydn)
allvar<-c(rep(-99,dim(alldat)[2]))
for(i in 1:dim(alldat)[2]) allvar[i]<-var(alldat[,i])
sum(allvar<0.2)
#[1] 602
allvar0<-which(allvar<0.2)
allnew<-alldat[,setdiff(1:dim(alldat)[2],allvar0)]
# zero for renaissance, one for baroque
class<-c(rep(0,dim(moz)[1]),rep(1,dim(haydn)[1]))
## now remove any variables that are highly correlated (>0.9)
allnew3<-allnew
ind<-1
while(ind>0){
	step<-0
	tmp<-cor(allnew3)
	check<-NULL
	for(i in 1:(dim(tmp)[1]-1)){
		for(j in (i+1):dim(tmp)[1]){
			if(abs(tmp[i,j])>0.9){
			 	check<-rbind(check,c(i,j))
			 	step<-step+1
			}
		}
	}
	if(step>1) allnew3<-allnew3[,-check[1,2]]
	if(step==1) allnew3<-allnew3[,-check[2]]
	if(step==0) ind=0
}
dim(allnew3)
#[1] 344  25

# select 50% analysis set by random proportional selection 
allnew3<-cbind(class,allnew3)
set.seed(1045)
set0sel<-sample.int(dim(moz)[1],0.5*0.29*dim(allnew3)[1])
length(set0sel)
#[1] 49
set1sel<-sample.int(dim(haydn)[1],0.5*0.71*dim(allnew3)[1])
length(set1sel)
#[1] 122
set0sel<-sort(set0sel)
set1sel<-sort(set1sel)
trainnew<-allnew3[c(set0sel,(set1sel+dim(moz)[1])),]
dim(trainnew)
#[1] 171 26
testnew<-allnew3[setdiff(1:dim(allnew3)[1],c(set0sel,(set1sel+dim(moz)[1]))),]
dim(testnew)
#[1] 173	26

# remove variables with sqrt(var)<0.2 in training from both train and test
trainvar<-rep(0,dim(trainnew)[2])
for(i in 1:dim(trainnew)[2]) trainvar[i]<-var(trainnew[,i])
trainvar0<-which(trainvar<0.2)
trainnew2<-trainnew[,setdiff(1:dim(trainnew)[2],trainvar0)]
testnew2<-testnew[,setdiff(1:dim(testnew)[2],trainvar0)]

# remove highly correlated variables in training from both train and test
## now remove any variables that are highly correlated (>0.9)
trainnew3<-trainnew2
testnew3<-testnew2
ind<-1
while(ind>0){
	step<-0
	tmp<-cor(trainnew3)
	check<-NULL
	for(i in 1:(dim(tmp)[1]-1)){
		for(j in (i+1):dim(tmp)[1]){
			if(abs(tmp[i,j])>0.9){
			 	check<-rbind(check,c(i,j))
			 	step<-step+1
			}
		}
	}
	if(step>1){ 
		trainnew3<-trainnew3[,-check[1,2]]
		testnew3<-testnew3[,-check[1,2]]
		}
	if(step==1){
		 trainnew3<-trainnew3[,-check[2]]
		 testnew3<-testnew3[,-check[2]]
		 }
	if(step==0) ind=0
}
dim(trainnew3)
#[1] 171  25

glm.out <- glm(class ~ ., family=binomial(logit), data=trainnew3)
summary(glm.out)
Call:
glm(formula = class ~ ., family = binomial(logit), data = trainnew3)
Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-2.9625  -0.1299   0.2514   0.5415   1.7236  
Coefficients:
              Estimate Std. Error z value Pr(>|z|)    
(Intercept) -50.861826  15.453814  -3.291 0.000998 ***
V130          0.113581   0.606543   0.187 0.851457    
V131         -0.758335   0.341439  -2.221 0.026351 *  
V135         -0.207635   0.304853  -0.681 0.495810    
V145         -0.039893   0.012775  -3.123 0.001792 ** 
V146         -0.049845   0.106908  -0.466 0.641039    
V415         -0.123732   0.073386  -1.686 0.091789 .  
V417          0.074978   0.107542   0.697 0.485682    
V424          0.017394   0.013125   1.325 0.185096    
V425         -1.142426   0.397304  -2.875 0.004035 ** 
V426          0.894759   0.307795   2.907 0.003649 ** 
V428         -0.285261   0.667468  -0.427 0.669104    
V431          0.327891   0.427792   0.766 0.443395    
V438         -0.035291   0.062473  -0.565 0.572143    
V439         -0.061913   0.221943  -0.279 0.780275    
V440          0.314889   0.324349   0.971 0.331631    
V441          0.201974   0.113566   1.778 0.075328 .  
V442         -1.019466   0.633040  -1.610 0.107304    
V443          0.031295   0.113748   0.275 0.783221    
V445          0.891988   0.229067   3.894 9.86e-05 ***
V449         -0.105012   0.080591  -1.303 0.192566    
V605         -0.129924   0.080497  -1.614 0.106523    
V609         -0.004087   0.011952  -0.342 0.732387    
V610          0.127104   0.068623   1.852 0.063995 .  
V632         -3.216045 219.371595  -0.015 0.988303    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Dispersion parameter for binomial family taken to be 1)
    Null deviance: 204.87  on 170  degrees of freedom
Residual deviance: 121.42  on 146  degrees of freedom
AIC: 171.42
Number of Fisher Scoring iterations: 16

table(fitted=round(fitted(glm.out)),true=trainnew3$class)
      true
fitted   0   1
     0  34  10
     1  15 112
glm.pred<-predict(glm.out, newdata=testnew3[,2:25],type="response")
table(pred=round(glm.pred),true=testnew3$class)
    true
pred  0  1
   0 24 25
   1 27 97

library(fpc)
clust<-pamk(trainnew3[,2:25],critout=TRUE)
2  clusters  0.4141576 
3  clusters  0.3200186 
4  clusters  0.4003161 
5  clusters  0.3905701 
6  clusters  0.3752969 
7  clusters  0.3413487 
8  clusters  0.3426997 
9  clusters  0.3481117 
10  clusters  0.2786395 

# svm
library(kernlab)
svp<-ksvm(class ~ ., data=trainnew3,type="C-svc",kernel="rbfdot", C=1)
svp
Support Vector Machine object of class "ksvm" 
SV type: C-svc  (classification) 
 parameter : cost C = 1 
Gaussian Radial Basis kernel function. 
 Hyperparameter : sigma =  0.0274483459932135 
Number of Support Vectors : 121 
Objective Function Value : -75.8909 
Training error : 0.175439 

svp.pred<-predict(svp, testnew3[,2:25])
table(pred = svp.pred, true = testnew3$class)
    true
pred   0   1
   0   8   0
   1  43 122

# rpart
 library(rpart)
rpart.model<-rpart(class ~ ., data=trainnew3, method="class")
summary(rpart.model)
Call:
rpart(formula = class ~ ., data = trainnew3, method = "class")
  n= 171 
          CP nsplit rel error   xerror      xstd
1 0.07142857      0 1.0000000 1.000000 0.1206657
2 0.04081633      4 0.7142857 1.244898 0.1278401
3 0.03061224      5 0.6734694 1.326531 0.1295435
4 0.01000000      7 0.6122449 1.326531 0.1295435
Variable importance
V445 V443 V145 V146 V610 V441 V417 V449 V609 V605 V131 V130 V439 V425 
  17   11   11    8    8    8    8    7    6    6    4    3    1    1 
Node number 1: 171 observations,    complexity param=0.07142857
  predicted class=1  expected loss=0.2865497  P(node) =1
    class counts:    49   122
   probabilities: 0.287 0.713 
  left son=2 (106 obs) right son=3 (65 obs)
  Primary splits:
      V445 < 64.47253 to the left,  improve=5.604341, (0 missing)
      V145 < 145.9068 to the right, improve=4.752623, (0 missing)
      V146 < 33.375   to the right, improve=4.752623, (0 missing)
      V443 < 51.5     to the left,  improve=3.772785, (0 missing)
      V431 < 1.5      to the right, improve=3.315292, (0 missing)
  Surrogate splits:
      V443 < 57.5     to the left,  agree=0.661, adj=0.108, (0 split)
      V130 < 2.350319 to the right, agree=0.655, adj=0.092, (0 split)
      V131 < 0.5      to the right, agree=0.649, adj=0.077, (0 split)
      V145 < 19.17725 to the right, agree=0.649, adj=0.077, (0 split)
      V605 < 8.5      to the right, agree=0.649, adj=0.077, (0 split)
Node number 2: 106 observations,    complexity param=0.07142857
  predicted class=1  expected loss=0.3867925  P(node) =0.619883
    class counts:    41    65
   probabilities: 0.387 0.613 
  left son=4 (50 obs) right son=5 (56 obs)
  Primary splits:
      V449 < 6        to the right, improve=4.443019, (0 missing)
      V417 < 4.375    to the left,  improve=3.820056, (0 missing)
      V145 < 120.0387 to the right, improve=3.316208, (0 missing)
      V443 < 52.5     to the left,  improve=2.689750, (0 missing)
      V441 < 48.5     to the left,  improve=2.590267, (0 missing)
  Surrogate splits:
      V130 < 2.79324  to the left,  agree=0.660, adj=0.28, (0 split)
      V439 < 6        to the right, agree=0.623, adj=0.20, (0 split)
      V445 < 63.50926 to the right, agree=0.623, adj=0.20, (0 split)
      V145 < 104.7413 to the right, agree=0.604, adj=0.16, (0 split)
      V425 < 2.5      to the left,  agree=0.594, adj=0.14, (0 split)
Node number 3: 65 observations
  predicted class=1  expected loss=0.1230769  P(node) =0.380117
    class counts:     8    57
   probabilities: 0.123 0.877 
Node number 4: 50 observations,    complexity param=0.07142857
  predicted class=0  expected loss=0.46  P(node) =0.2923977
    class counts:    27    23
   probabilities: 0.540 0.460 
  left son=8 (10 obs) right son=9 (40 obs)
  Primary splits:
      V145 < 106.9582 to the right, improve=3.240000, (0 missing)
      V131 < 1.5      to the right, improve=2.755584, (0 missing)
      V609 < 180.5    to the right, improve=2.671978, (0 missing)
      V146 < 28.625   to the right, improve=1.867778, (0 missing)
      V449 < 9.5      to the left,  improve=1.831790, (0 missing)
  Surrogate splits:
      V609 < 175.5    to the right, agree=0.92, adj=0.6, (0 split)
      V441 < 52.5     to the right, agree=0.88, adj=0.4, (0 split)
      V605 < 17.5     to the right, agree=0.88, adj=0.4, (0 split)
      V146 < 28.875   to the right, agree=0.86, adj=0.3, (0 split)
      V610 < 39.5     to the right, agree=0.86, adj=0.3, (0 split)
Node number 5: 56 observations,    complexity param=0.03061224
  predicted class=1  expected loss=0.25  P(node) =0.3274854
    class counts:    14    42
   probabilities: 0.250 0.750 
  left son=10 (26 obs) right son=11 (30 obs)
  Primary splits:
      V445 < 62.95865 to the left,  improve=4.343590, (0 missing)
      V443 < 53.5     to the left,  improve=4.200000, (0 missing)
      V146 < 27.08333 to the left,  improve=2.907692, (0 missing)
      V417 < 4.522727 to the left,  improve=2.907692, (0 missing)
      V441 < 49.5     to the left,  improve=2.583653, (0 missing)
  Surrogate splits:
      V443 < 53.5     to the left,  agree=0.768, adj=0.500, (0 split)
      V145 < 34.46919 to the left,  agree=0.714, adj=0.385, (0 split)
      V605 < 12.5     to the left,  agree=0.714, adj=0.385, (0 split)
      V610 < 26.5     to the left,  agree=0.714, adj=0.385, (0 split)
      V146 < 26.03125 to the left,  agree=0.679, adj=0.308, (0 split)
Node number 8: 10 observations
  predicted class=0  expected loss=0.1  P(node) =0.05847953
    class counts:     9     1
   probabilities: 0.900 0.100 
Node number 9: 40 observations,    complexity param=0.07142857
  predicted class=1  expected loss=0.45  P(node) =0.2339181
    class counts:    18    22
   probabilities: 0.450 0.550 
  left son=18 (24 obs) right son=19 (16 obs)
  Primary splits:
      V441 < 48.5     to the left,  improve=3.675000, (0 missing)
      V610 < 31.5     to the left,  improve=3.000000, (0 missing)
      V417 < 6.214286 to the left,  improve=2.667384, (0 missing)
      V609 < 113.5    to the left,  improve=2.526817, (0 missing)
      V415 < 11.5     to the left,  improve=2.495341, (0 missing)
  Surrogate splits:
      V443 < 52.5     to the left,  agree=0.875, adj=0.688, (0 split)
      V609 < 116.5    to the left,  agree=0.825, adj=0.562, (0 split)
      V417 < 2.982143 to the left,  agree=0.800, adj=0.500, (0 split)
      V610 < 28       to the left,  agree=0.775, adj=0.438, (0 split)
      V146 < 26.625   to the left,  agree=0.750, adj=0.375, (0 split)
Node number 10: 26 observations,    complexity param=0.03061224
  predicted class=1  expected loss=0.4615385  P(node) =0.1520468
    class counts:    12    14
   probabilities: 0.462 0.538 
  left son=20 (19 obs) right son=21 (7 obs)
  Primary splits:
      V417 < 5.25     to the left,  improve=1.9456330, (0 missing)
      V145 < 34.68095 to the right, improve=1.8754580, (0 missing)
      V131 < 1.5      to the right, improve=1.1655010, (0 missing)
      V146 < 27.08333 to the left,  improve=0.8480769, (0 missing)
      V441 < 49.5     to the left,  improve=0.8480769, (0 missing)
  Surrogate splits:
      V443 < 51.5     to the left,  agree=0.846, adj=0.429, (0 split)
      V146 < 27.08333 to the left,  agree=0.808, adj=0.286, (0 split)
      V130 < 3.106136 to the left,  agree=0.769, adj=0.143, (0 split)
      V131 < 0.5      to the right, agree=0.769, adj=0.143, (0 split)
      V605 < 15.5     to the left,  agree=0.769, adj=0.143, (0 split)
Node number 11: 30 observations
  predicted class=1  expected loss=0.06666667  P(node) =0.1754386
    class counts:     2    28
   probabilities: 0.067 0.933 
Node number 18: 24 observations,    complexity param=0.04081633
  predicted class=0  expected loss=0.375  P(node) =0.1403509
    class counts:    15     9
   probabilities: 0.625 0.375 
  left son=36 (12 obs) right son=37 (12 obs)
  Primary splits:
      V131 < 1.5      to the right, improve=2.083333, (0 missing)
      V415 < 16.5     to the right, improve=1.500000, (0 missing)
      V146 < 24.875   to the right, improve=1.180070, (0 missing)
      V130 < 2.44245  to the right, improve=0.762605, (0 missing)
      V440 < 1.5      to the left,  improve=0.762605, (0 missing)
  Surrogate splits:
      V146 < 25.25    to the right, agree=0.750, adj=0.500, (0 split)
      V417 < 1.875    to the right, agree=0.750, adj=0.500, (0 split)
      V443 < 48.5     to the right, agree=0.750, adj=0.500, (0 split)
      V610 < 17.5     to the right, agree=0.708, adj=0.417, (0 split)
      V145 < 38.71413 to the right, agree=0.667, adj=0.333, (0 split)
Node number 19: 16 observations
  predicted class=1  expected loss=0.1875  P(node) =0.09356725
    class counts:     3    13
   probabilities: 0.188 0.812 
Node number 20: 19 observations
  predicted class=0  expected loss=0.4210526  P(node) =0.1111111
    class counts:    11     8
   probabilities: 0.579 0.421 
Node number 21: 7 observations
  predicted class=1  expected loss=0.1428571  P(node) =0.04093567
    class counts:     1     6
   probabilities: 0.143 0.857 
Node number 36: 12 observations
  predicted class=0  expected loss=0.1666667  P(node) =0.07017544
    class counts:    10     2
   probabilities: 0.833 0.167 
Node number 37: 12 observations
  predicted class=1  expected loss=0.4166667  P(node) =0.07017544
    class counts:     5     7
   probabilities: 0.417 0.583 
rpart.pred<-predict(rpart.model,testnew3[,2:25],type="class")
table(pred = rpart.pred, true = testnew3[,1])
    true
pred   0   1
   0  19  13
   1  32 109
library(rpart.plot)
prp(rpart.model)
prp(rpart.model, under.cex=1.2, split.cex=1.2)

# random forest
library(randomForest)
library(Hmisc)
set.seed(105)
rf<-randomForest(class ~ ., data=trainnew3, ntree=100, mtry=5, importance=T, proximity=T)
print(rf)
Call:
 randomForest(formula = class ~ ., data = trainnew3, ntree = 100, mtry = 5, importance = T, proximity = T) 
               Type of random forest: classification
                     Number of trees: 100
No. of variables tried at each split: 5
        OOB estimate of  error rate: 25.73%
Confusion matrix:
   0   1 class.error
0 14  35  0.71428571
1  9 113  0.07377049
rf.pred<-predict(rf, newdata=testnew3)
table(pred=rf.pred, true=testnew3$class)
    true
pred   0   1
   0  17   6
   1  34 116

## k-NN (k=1 default)
library(RWeka)
install.packages("RWeka", dependencies=T)
knn<-IBk(class ~ ., data=trainnew3)
 summary(knn)
=== Summary ===
Correctly Classified Instances         171              100      %
Incorrectly Classified Instances         0                0      %
Kappa statistic                          1     
Mean absolute error                      0.0058
Root mean squared error                  0.0058
Relative absolute error                  1.4101 %
Root relative squared error              1.2784 %
Coverage of cases (0.95 level)         100      %
Mean rel. region size (0.95 level)      50      %
Total Number of Instances              171     
=== Confusion Matrix ===
   a   b   <-- classified as
  49   0 |   a = 0
   0 122 |   b = 1
knn.pred<-predict(knn, newdata=testnew3)
table(pred=knn.pred, true=testnew3$class)
    true
pred   0   1
   0  24  21
   1  27 101
# for info about options for IBk use WOW("IBk")
 knn2<-IBk(class ~ ., data=trainnew3, control=Weka_control(K = 2, I=TRUE))
summary(knn2)
=== Summary ===
Correlation coefficient                  1     
Mean absolute error                      0.0024
Root mean squared error                  0.0045
Relative absolute error                  0.5946 %
Root relative squared error              0.986  %
Total Number of Instances              171     
knn2.fit<-predict(knn2, newdata=trainnew3[,2:25], type="class")
 table(fit=round(knn2.fit), true=trainnew3$class)
   true
fit   0   1
  0  49   0
  1   0 122
knn2.pred<-predict(knn2, newdata=testnew3[,2:25], type="class")
> table(fit=round(knn2.pred), true=testnew3$class)
   true
fit   0   1
  0  24  21
  1  27 101
knn3<-IBk(class ~ ., data=trainnew3, control=Weka_control(K = 3, I=TRUE))
knn3.fit<-predict(knn3, newdata=trainnew3[,2:25], type="class")
table(fit=round(knn3.fit), true=trainnew3$class)
   true
fit   0   1
  0  49   0
  1   0 122
knn3.pred<-predict(knn3, newdata=testnew3[,2:25], type="class")
table(fit=round(knn3.pred), true=testnew3$class)
   true
fit   0   1
  0  15  11
  1  36 111

# gradient boosted tree
library(gbm)
gbm.out<-gbm(class ~ ., data=trainnew3, n.trees=1000, shrinkage=0.05, interaction.depth=2, train.fraction=1, n.minobsinnode=5, cv.folds=5, verbose=T)
best.iter <- gbm.perf(gbm.out,method="cv")
f.fit<-predict(gbm.out, trainnew3, best.iter, type="response")
table(fitted=round(f.fit), true=trainnew3$class)
      true
fitted   0   1
     0  34   3
     1  15 119
f.predict <- predict(gbm.out,testnew3,best.iter,type="response")
 table(pred=round(f.predict), true=testnew3$class)
#    true
#pred   0   1
#   0  21  11
#   1  30 111
