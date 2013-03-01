library(spatialprobit)
library(McSpatial)

set.seed(9947)
library(maptools)
cmap <- readShapePoly(system.file("maps/CookCensusTracts.shp",
  package="McSpatial"))
cmap <- cmap[cmap$CHICAGO==1&cmap$CAREA!="O'Hare",]
plot(cmap)
wmat <- makew(cmap)$wmat
n = nrow(wmat)
rho <- .4
x <- runif(n,0,10)
ystar <- as.numeric(solve(diag(n) - rho*wmat)%*%(x + rnorm(n,0,2)))
y <- ystar>quantile(ystar,.4)

################################################################################
#
# GMM fit with McSpatial
#
###############################################################################
system.time(fit1 <- spprobit(y~x,  wmat=wmat))
# 0.03

#STANDARD PROBIT ESTIMATES 
#LINEARIZED GMM PROBIT ESTIMATES 
#            Estimate Std. Error   z-value Pr(>|z|)
#(Intercept) -2.27695    0.12206 -18.65500    0e+00
#x            0.50788    0.02434  20.86551    0e+00
#WXB          0.44670    0.09812   4.55247    1e-05
write(wmat, ncolumns=NCOL(wmat), file="C:/Projects/MATLAB/McSpatialTest/wmat.dat")

# fit with spatialprobit
# wmat is a 861 x 861 dense matrix. D.h. spprobit will W als dense matrix haben
# Wenn ich wmat als dense reingebe, dauert es ewig
W <- Matrix(wmat, sparse=TRUE)
fit2 <- sarprobit(y ~ x, W=W, m=20, showProgress=TRUE, 
  start = list(rho = 0.5, beta=c(0, 0)), ndraw=2000, burn.in=500)
summary(fit2)

par(mfrow=c(3,3))
plot(fit2)

# 42 Sekunden
#            Estimate Std. Dev  p-level t-value Pr(>|z|)    
#(Intercept) -2.21934  0.15423  0.00000 -14.389  < 2e-16 ***
#x            0.51319  0.03035  0.00000  16.909  < 2e-16 ***
#rho          0.29025  0.07056  0.00000   4.114 4.27e-05 ***

################################################################################
#
# Kleine Monte-Carlo-Studie: spatialprobit vs. McSpatial
#
################################################################################

# 100 mal wiederholen
setwd("C:/Projects/spatialprobit/")
set.seed(1.2345)
fits.spa <- list()       # spatialprobit MCMC
fits.mcs.gmm <- list()   # McSpatial GMM
fits.mcs.mle   <- list() # McSpatial MLE
for (i in 1:100) {
  cat(format(Sys.time()),i,"/",100,"\n")
  flush.console()
  x <- runif(n,0,10)
  ystar <- as.numeric(solve(diag(n) - rho*wmat)%*%(x + rnorm(n,0,2)))
  y <- ystar>quantile(ystar,.4)
  
  # McSpatial GMM
  fits.mcs.gmm[[i]] <- spprobit(y~x,  wmat=wmat)
  
  # McSpatial MLE
  fits.mcs.mle[[i]] <- spprobitml(y~x,  wmat=wmat)

  # sarprobit
  fits.spa[[i]] <- sarprobit(y ~ x, W=W, m=10, showProgress=FALSE, computeMarginalEffects=FALSE,
    start = list(rho = 0.5, beta=c(0, 0)), ndraw=1000, burn.in=200)
}

params <- c(-2.2, 0.5, 0.4)

# spatialprobit MCMC
coefMat.spa <- t(sapply(fits.spa, coef))
colMeans(coefMat.spa)

# Bias
bias.spa <- colMeans(coefMat.spa) - params
bias.spa

##

# McSpatial GMM
coefMat.mcs.gmm <- t(sapply(fits.mcs.gmm, function(x) x$coef))
colMeans(coefMat.mcs.gmm)

# Bias
bias.mcs.gmm <- colMeans(coefMat.mcs.gmm) - params
bias.mcs.gmm

##

# McSpatial MLE
coefMat.mcs.mle <- t(sapply(fits.mcs.mle, function(x) x$coef))
colMeans(coefMat.mcs.mle)

# Bias
bias.mcs.mle <- colMeans(coefMat.mcs.mle) - params
bias.mcs.mle


par(mfrow=c(2,2))
d1 <- density(coefMat.spa[,1])
d2 <- density(coefMat.mcs.gmm[,1])
d3 <- density(coefMat.mcs.mle[,1])
plot(d1, main=expression(beta[1]), ylim=range(d1$y, d2$y, d3$y))
lines(d2, col="red")
lines(d3, col="green")
abline(v=params[1], lty=2, col="red")
legend("topright", legend=c("spatialprobit MCMC","McSpatial GMM","McSpatial MLE"), 
  col=c("black","red","green"), lty=1, bty="n")

d1 <- density(coefMat.spa[,2])
d2 <- density(coefMat.mcs.gmm[,2])
d3 <- density(coefMat.mcs.mle[,2])
plot(d1, main=expression(beta[2]), ylim=range(d1$y, d2$y, d3$y))
lines(d2, col="red")
lines(d3, col="green")
abline(v=params[2], lty=2, col="red")
legend("topright", legend=c("spatialprobit MCMC","McSpatial GMM","McSpatial MLE"), 
  col=c("black","red","green"), lty=1, bty="n")

d1 <- density(coefMat.spa[,3])
d2 <- density(coefMat.mcs.gmm[,3])
d3 <- density(coefMat.mcs.mle[,3])
plot(d1, main=expression(rho), xlim=range(d1$x, d2$x), ylim=range(d1$y, d2$y, d3$y))
lines(d2, col="red")
lines(d3, col="green")
abline(v=params[3], lty=2, col="red")
legend("topright", legend=c("spatialprobit MCMC","McSpatial GMM","McSpatial MLE"), 
  col=c("black","red","green"), lty=1, bty="n")

save(fits.spa, coefMat.spa,
     fits.mcs.gmm, coefMat.mcs.gmm,  
     fits.mcs.mle, coefMat.mcs.mle,  
     file="McSpatial-test-2013-02-24.RData")

#--------MCMC spatial autoregressive probit--------
#Execution time  = 22.172 secs
#
#N draws         =   1000, N omit (burn-in)=    100
#N observations  =    861, K covariates    =      2
## of 0 Y values =    345, # of 1 Y values =    516
#Min rho         = -1.000, Max rho         =  1.000
#--------------------------------------------------
#
#            Estimate Std. Dev Bayes p-level t-value Pr(>|z|)    
#(Intercept) -2.23936  0.15467       0.00000  -14.48  < 2e-16 ***
#x            0.51774  0.03197       0.00000   16.19  < 2e-16 ***
#rho          0.28647  0.07234       0.00000    3.96  8.1e-05 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

# --> McSpatial GMM scheint in rho verzerrt zu sein.
# --> see also Franzese (2013)

################################################################################
#
# McSpatial ML
#
################################################################################

# ML fit with McSpatial
system.time(
fit2 <- spprobitml(y~x,  wmat=wmat)
)
# 12.20 seconds

W <- Matrix(wmat, sparse=TRUE)
fit2 <- sarprobit(y ~ x, W=W, m=10, showProgress=FALSE, computeMarginalEffects=FALSE,
    start = list(rho = 0.5, beta=c(0, 0)), ndraw=1000, burn.in=200)
fit2$time    
# Execution time  = 18.081 secs

# --> Fazit: McSpatial MLE ist nur unwesentlich schneller als MCMC, möglicherweise aber
# biased in rho


# --> Wahrscheinlich gewinnt McSpatial bei kleinen Samples. Wie sieht es bei
# z.B. n=2000 aus? McSpatial nutzt nur dense matrix, daher wird es Probleme bekommen.

