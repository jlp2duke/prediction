
tour <- read.csv("tourdefranceclean_WWI.csv", header= T) 
lnspeed<-log(tour$SPEED)
DIST2 <- tour$DIST*tour$DIST
Y2 = tour$YEAR*tour$YEAR
Y3 = (tour$YEAR)^3
DISTY = tour$DIST * tour$YEAR
plot(tour$YEAR, lnspeed)
plot(tour$DIST, lnspeed)

lmmod <- lm(formula = lnspeed ~ tour$YEAR + Y2 + tour$DIST + DIST2 , data = tour)
#scatter plot of the data
library(lattice)
YEAR = tour$YEAR
LENGTH = tour$DIST
LSPEED = lnspeed
cloud(LSPEED ~ YEAR * LENGTH , data = tour, screen = list(z = -200, x = -100, y = -200))
#here is the estimated regression function from lmmod
f <- function(x,y){ 
value <- -41.4 + .04121*x  -.000009573*x^2 + .0004098*y - .00000008263*y^2
return(value)
}
xval = c(seq(1903,2012,1))
yval = c(seq(1400,3500, 1))
pairs = expand.grid(xval, yval)
zval =f(pairs[,1], pairs[,2])
zmat = matrix(zval, 110, 2101)
persp(xval, yval, zmat, theta = 30, phi = 45) # better angle of view


