# ch10.3 computing

x<-seq(0,12,0.01)
y<-function(x){.5 + 1/(1+ exp(-5*(x -1) )) + 2/(1+ exp(-5*(x-3))) + 3/(1+ exp(-5*(x-5))) + 4/(1+ exp(-5*(x-7))) + 5/(1+ exp(-5*(x-9)))}
z<-y(x)
plot(x,z,xlab="input",ylab="output",cex.axis=1.2, cex.lab=1.2)
# save as firstNNplot_v2.pdf

