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
rho = .4
x <- runif(n,0,10)
ystar <- as.numeric(solve(diag(n) - rho*wmat)%*%(x + rnorm(n,0,2)))
y <- ystar>quantile(ystar,.4)

# fit with McSpatial
fit <- spprobit(y~x,  wmat=wmat)

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

################################################################################

# 100 mal wiederholen
setwd("C:/Projects/spatialprobit/")
fits <- list()
for (i in 1:100) {
  cat(format(Sys.time()),i,"/",100,"\n")
  flush.console()
  x <- runif(n,0,10)
  ystar <- as.numeric(solve(diag(n) - rho*wmat)%*%(x + rnorm(n,0,2)))
  y <- ystar>quantile(ystar,.4)

  fits[[i]] <- sarprobit(y ~ x, W=W, m=10, showProgress=FALSE, 
    start = list(rho = 0.5, beta=c(0, 0)), ndraw=1000, burn.in=200)
}
coefMat <- t(sapply(fits, coef))
coefMat

par(mfrow=c(2,2))
plot(density(coefMat[,1]))
abline(v=params[1], lty=2, col="red")
plot(density(coefMat[,2]))
abline(v=params[2], lty=2, col="red")
plot(density(coefMat[,3]))
abline(v=params[3], lty=2, col="red")

save(fits, file="McSpatial-test-2013-01-25.RData")

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

# --> es sieht alles ein wenig verzerrt aus!
# Die Schätzung für rho ist nicht komplett daneben 
# (zumindestens noch innerhalb von 2SD um den wahren Wert 0.4), 
# aber immer irgendwie etwas ab!
# Die Standardfehler dagegen stimmen eigentlich gut überein!
# Ist da beim Ziehen von Rho noch irgendwo ein Fehler???
# TODO: Testen Ziehen von Rho mit M-H vs. anderes Verfahren
# Es hat m.E. nix mit m=1 oder m=10 oder m=20 zu tun --> das scheint ok zu sein.
# rho sieht verzerrt aus!
# TODO: Ziehe aus rho | beta, z, y für bekanntes rho, beta und z

