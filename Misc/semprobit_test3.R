################################################################################
#
# Estimate SEM Probit from simulated data
#
################################################################################
library(spatialprobit)
if (FALSE) {
 library(tmvtnorm)
 library(Matrix)    # sparseMatrix
 library(spdep)
 setwd("C:/Projects/R/spatialprobit/R")
 source("sar_base.r")
 source("matrix_operations.r")
 source("stats_distributions.r")
 source("utility_functions.r")
 source("rtnorm.r")
 source("semprobit.r")

}
set.seed(1)

n <- 300
latt <- rnorm(n)
long <- rnorm(n)
W <- kNearestNeighbors(latt,long,6)

# true parameters
sige <- 0.5               # varianz --> std = 0.25
k <- 3
rho <- 0.7
beta <- c(-0.5, 2, -2)

# simulate data
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
AI <- solve(I_n-rho*W)
X <- matrix(rnorm(n*k), n, k)
z <- X%*%beta + AI %*% rnorm(n) * sqrt(sige);
y <- as.numeric(z >= 0)

# estimate SEM probit from simulated data
semp_fit <- sem_probit_mcmc(y, X, W, ndraw=200, burn.in=100, thinning=1, showProgress=TRUE, univariateConditionals=TRUE)
summary(semp_fit)
plot(semp_fit)

semp_fit2 <- sem_probit_mcmc(y, X, W, ndraw=2000, burn.in=100, thinning=1, showProgress=TRUE, univariateConditionals=FALSE)
summary(semp_fit2)
plot(semp_fit2)