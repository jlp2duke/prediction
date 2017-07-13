# ch10.2 computing

library(RColorBrewer)
library(colorRamps)
#layout(rbind(c(1,2),c(3,2)),widths=c(2,3.5),heights=c(2,1))
plot.new()
# basic color blocks
rect(1/2,1/2,0,0,col='red')
rect(1/2,1,0,1/2,col='green')
rect(1/2,1,1/3,1/2,col='blue')
rect(1,1,1/2,1/3,col='yellow')
rect(1,0,1/2,1/3,col='purple')
# add text
mtext(expression(C[1]), side = 1, line = 0, col='blue')
mtext(expression(X[1]), side = 1, line = 2, cex=1.5)
text(-0.05,1/2,expression(C[2]), col='blue', srt=0, xpd=TRUE)
text(-0.15,1/2,expression(X[2]), srt=0, xpd=TRUE, cex=1.5)
mtext(expression(C[3]), col='blue', side = 3, line = 0, adj=1/3)
mtext(expression(X[3]), side = 3, line = 0, cex=1.5)
text(1.05,1/3,expression(C[4]), col='blue', srt=0, xpd=TRUE)
text(1.05,1/2,expression(X[4]), srt=0, xpd=TRUE, cex=1.5)
 
plot.new()
 # root node branches
lines(c(1/2,1/6),c(1-1/10,2/3))
lines(c(1/2,5/6),c(1-1/10,2/3))
# left of root branches
lines(c(0,1/6),c(1/3,2/3-1/10))
lines(c(1/6,1/3),c(2/3-1/10,1/3))
# right of root branches
lines(c(2/3,5/6),c(1/3,2/3-1/10))
lines(c(5/6,1),c(2/3-1/10,1/3))
# bottom branches
lines(c(1/4,1/3),c(0,1/3-1/10))
lines(c(1/3,5/12),c(1/3-1/10,0))
text(1/2,1-1/20,paste(expression(X[1]),"<>",expression(C[1])))
text(1/6,2/3-1/20,paste(expression(X[2]),"<>",expression(C[2])))
text(5/6,2/3-1/20,paste(expression(X[3]),"<>",expression(C[3])))
text(1/3,1/3-1/20,paste(expression(X[4]),"<>",expression(C[4])))	#plot2.pdf

par(mfrow=c(2,1))
y<-sort(rep(seq(10,100,10),2))
x<-sort(rep(seq(5,105,10),2))
x<-x[-c(1,22)]
plot(x,y,type="n",main="little data")
for(i in 1:10){
	points(x[(2*(i-1)+1):(2*i)],y[(2*(i-1)+1):(2*i)])
	lines(x[(2*(i-1)+1):(2*i)],y[(2*(i-1)+1):(2*i)])
}
lines(c(5,105),c(5,105))
y<-sort(rep(seq(10,100,2),2))
x<-sort(rep(seq(9,101,2),2))
x<-x[-c(1,94)]
plot(x,y,type="n",main="lots of data")
for(i in 1:46){
	points(x[(2*(i-1)+1):(2*i)],y[(2*(i-1)+1):(2*i)])
	lines(x[(2*(i-1)+1):(2*i)],y[(2*(i-1)+1):(2*i)])
}
lines(c(5,105),c(5,105))	# plot4.pdf
# save as staircase.pdf

library(lattice)
g <- expand.grid(x = c(seq(-100,0,5),0.1,seq(5,100,5)), y = c(seq(-100,0,5),0.1,seq(5,100,5)))
g$z<-matrix(y,1764,1)
myColours<-colorRampPalette(brewer.pal(9,"YlOrRd"), space = "Lab")
wireframe(z~x*y,data=g, drape = TRUE, colorkey = TRUE, screen = list(z = -60, x = -60),par.settings=list(regions=list(alpha=1,col=myColours(100))))


