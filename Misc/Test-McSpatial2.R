# Monte Carlo Study of McSpatial GMM estimation
library(spatialprobit)
library(McSpatial)

setwd("C:/Projects/spatialprobit/")
set.seed(1.2345)

# LeSagePaceExperiment
n <- 400
beta <- c(0, 1, -1)
rho <- 0.75

X <- cbind(intercept = 1, x = rnorm(n), y = rnorm(n))
I_n <- sparseMatrix(i = 1:n, j = 1:n, x = 1)
W <- kNearestNeighbors(x = rnorm(n), y = rnorm(n), k = 6)
wmat <- as.matrix(W)

# 100 mal wiederholen
fits.mcs.gmm <- list()   # McSpatial GMM
fits.mcs.mle   <- list()  # McSpatial MLE
for (i in 1:100) {
  cat(format(Sys.time()),i,"/",100,"\n")
  flush.console()
  eps <- rnorm(n = n, mean = 0, sd = 1)
  z <- solve(qr(I_n - rho * W), X %*% beta + eps)
  y <- as.vector(z >= 0)
  
  # McSpatial GMM
  fits.mcs.gmm[[i]] <- spprobit(y~X-1,  wmat=wmat)
  
  # McSpatial MLE
  fits.mcs.mle[[i]] <- spprobitml(y~X-1,  wmat=wmat)
}

params <- c(beta, rho)

# McSpatial GMM
coefMat.mcs.gmm <- t(sapply(fits.mcs.gmm, function(x) x$coef))
colMeans(coefMat.mcs.gmm)

# Bias
bias.mcs.gmm <- colMeans(coefMat.mcs.gmm) - params
bias.mcs.gmm

# McSpatial MLE
coefMat.mcs.mle <- t(sapply(fits.mcs.mle, function(x) x$coef))
colMeans(coefMat.mcs.mle)

# Bias
bias.mcs.mle <- colMeans(coefMat.mcs.mle) - params
bias.mcs.mle


d2 <- density(coefMat.mcs.gmm[,4])
d3 <- density(coefMat.mcs.mle[,4])
plot(d2, main=expression(rho), col="red", xlim=range(d2$x, d3$x), ylim=range(d2$y, d3$y))
lines(d3, col="green")
abline(v=0.75, lty=2, col="red")
legend("topright", legend=c("spatialprobit MCMC","McSpatial GMM","McSpatial MLE"),
  col=c("black","red","green"), lty=1, bty="n")
