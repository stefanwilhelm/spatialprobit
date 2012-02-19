if (FALSE) {
################################################################################
#
# sarprobit and formula
#
################################################################################

CMK <- data.frame(y=c(0, 0, 1), x=c(1, 2, 3), income=c(1000, 500, 4000))
W <- sparseMatrix(i=1:3, j=c(3, 1, 2), x=c(1, 1, 1))
fit <- sarprobit(y ~ income, W, data=CMK, ndraw=100)

################################################################################
#
# example from LeSage(2009), section 10.1.5
#
################################################################################

# generate random samples from true model
n <- 200             # number of items
beta <- c(0, -1, 1)  # true model parameters k=3 beta=(beta1,beta2,beta3)
rho <- 0.75
# design matrix with two standard normal variates as "coordinates"
X <- cbind(intercept=1, x=rnorm(n), y=rnorm(n))
  
# identity matrix I_n
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
print(object.size(I_n), units="Mb")
  
# build spatial weight matrix W from coordinates in X
W <- kNearestNeighbors(x=rnorm(n), y=rnorm(n), k=6)
  
# create samples from epsilon using independence of distributions (rnorm()) to avoid dense matrix I_n
eps <- rnorm(n=n, mean=0, sd=1)
z <- solve(qr(I_n - rho * W), X %*% beta + eps)
y <- as.vector(z >= 0)  # binary observables, 0 or 1, FALSE or TRUE

# MCMC estimation of spatial probit
Rprof("sar_probit_mcmc.out")
results <- sar_probit_mcmc(y, X, W, ndraw=1000, burn.in=100, thinning=3)
Rprof(NULL)
summaryRprof("sar_probit_mcmc.out")

# visualization of sparsity pattern in W and H = (I_n - rho * W)'(I_n - rho * W)
print(image(W, main = "image(W)")) # print(.) needed for Sweave
S <- (I_n - rho * W)
print(image(S, main = "image(I_n - rho * W)")) # print(.) needed for Sweave
H <- t(S) %*% S
print(image(H, main = "image(H)")) # print(.) needed for Sweave
Sigma <- solve(H)
print(image(Sigma, main = "image(Sigma)")) # print(.) needed for Sweave

B   <- results$B
# parameter estimates for beta and rho
colMeans(B)
 
# standard errors for beta and rho           
apply(B, 2, sd)

cbind(estimates=colMeans(B), SE=apply(B, 2, sd))
  
getParamName <- function(i) {
  if (i==4) {
   name <- substitute(rho)
  } else {
   name <- substitute(beta[i],list(i=i))
  }
  return(name)
}

par(mfrow=c(2,2))
plot.results(results, trueparam = c(beta, rho))

# compile functions
#sar_probit_mcmc_comp <- cmpfun(sar_probit_mcmc)

# compiled MCMC function
#Rprof("sar_probit_mcmc_comp.out")
#results <- sar_probit_mcmc_comp(y, X, W, ndraw=1000, burn.in=100, thinning=3)
#Rprof(NULL)
#summaryRprof("sar_probit_mcmc_comp.out")

################################################################################
#
# example data set from Miguel Godinho de Matos
#
################################################################################

data <- read.table("I:/R/tmvtnorm/doc/Spatial Probit/Miguel/stefan/data/probit_demo.data")     
Y    <- as.matrix(data[,1] )
long <- as.matrix(data[,11])
latt <- as.matrix(data[,12])
X    <- as.matrix(data[,2:10])
W    <- make_neighborsw( latt, long , 15 )   # 15 nearest neighbors
# TO BE FIXED WHEN 2 POINTS HAVE THE SAME DISTANCE
# W    <- buildSpatialWeightMatrix(cbind(latt, long), m=15)   # 15 nearest neighbors
Rprof(filename="sar_probit_mcmc.out")
# ndraw=300, nomit=100
results <- sar_probit_mcmc(Y, X, W, ndraw=1000, burn.in=100, thinning=3, prior=NULL)
Rprof(NULL)
summaryRprof(filename="sar_probit_mcmc.out")

B   <- results$B
# parameter estimates for beta and rho
colMeans(B)
 
# standard errors for beta and rho           
apply(B, 2, sd)

pdf("probit_demo.pdf")
par(mfrow=c(2,2))
plot(results)
dev.off()
}

if (FALSE) {
################################################################################
#
# Debug method for tr(W^i)
#
################################################################################

# generate random samples from true model
n <- 200             # number of items
beta <- c(0, -1, 1)  # true model parameters k=3 beta=(beta1,beta2,beta3)
rho <- 0.75
# design matrix with two standard normal variates as "coordinates"
X <- cbind(intercept=1, x=rnorm(n), y=rnorm(n))
  
# identity matrix I_n
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
print(object.size(I_n), units="Mb")
  
# build spatial weight matrix W from coordinates in X
W <- kNearestNeighbors(x=rnorm(n), y=rnorm(n), k=6)

library(pracma)
Rprof("trW.out")
traceW.i <- tracesWi(W, o=100)
Rprof(NULL)
summaryRprof("trW.out")

################################################################################

undebug(sar_probit_mcmc)
debug(sar_probit_mcmc)
results <- sar_probit_mcmc(y, X, W, ndraw=1000, burn.in=0, thinning=1)
results <- sar_probit_mcmc(y, X, W, ndraw=1000, burn.in=0, thinning=1, showProgress=TRUE)
}