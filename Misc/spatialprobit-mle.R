library(spatialprobit)

# Estimate sarprobit with Maximum Likelihood;
# see Marsh et. al (2000), eqn (2.8) for Log-Likelihood function
# u = (I_n - rho * W)^{-1} eps, eps ~ N(0, sige*I_n)
# E[u u'] = sige * [(I_n - rho * W) * (I_n - rho * W)']^{-1}
# here we restrict sige = 1
loglik <- function(params, y, X, W) {
  k <- ncol(X)
  beta <- params[1:k]
  rho <- params[k+1]
  
  ind0 <- which(y == 0)
  ind1 <- which(y == 1)
  
  S <- I_n - rho * W
  D <- diag(1/sqrt(diag( S %*% t(S))))  # D = diag(E[u u'])^{1/2}  (n x n)
  SI <- solve(S)                        # (I_n - rho * W)^{-1}
  
  Xs <- D %*% SI %*% X                  # X^{*} = D * (I_n - rho * W)^{-1} * X
  F <- pnorm(as.double(Xs %*% beta))    # F(X^{*} beta)  # (n x 1)

  lnL <- sum(log(F[ind1])) + sum(log((1 - F[ind0]))) # see Marsh (2000), equation (2.8)
  return(lnL)
}

negloglik <- function(params, y, X, W) {
 -loglik(params, y, X, W)
}

sarprobitlm <- function(y, X, W) {
}

n <- 200

# true parameters
beta <- c(0, 1, -1)
rho <- 0.75

# design matrix with two standard normal variates as "covariates"
X <- cbind(intercept=1, x=rnorm(n), y=rnorm(n))

# sparse identity matrix
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)

# number of nearest neighbors in spatial weight matrix W
m <- 6

# spatial weight matrix with m=6 nearest neighbors
W <- sparseMatrix(i=rep(1:n, each=m),
  j=replicate(n, sample(x=1:n, size=m, replace=FALSE)), x=1/m, dims=c(n, n))

# innovations
eps <- rnorm(n=n, mean=0, sd=1)

# generate data from model
S <- I_n - rho * W
z <- solve(qr(S), X %*% beta + eps)
y <- as.double(z >= 0)  # 0 or 1, FALSE or TRUE

################################################################################
#
# Maximum-Likelihood estimation
#
################################################################################

# true params = c(0, 1, -1, 0.75)
p <- optim(c(0, 0, 0, 0.5), negloglik, hessian=TRUE,
  y=y, X=X, W=W)
# [1]  -0.04222603  1.11179922 -0.89755253  0.19169014

OI<-solve(p$hessian)
se<-sqrt(diag(OI))

res <- cbind(est=p$par, se)
rownames(res) <- c(colnames(X), "rho")

#                  est        se
#intercept -0.04222603 0.0964966
#x          1.11179922 0.1594369
#y         -0.89755253 0.1366937
#rho        0.19169014 0.1960016

################################################################################
#
# Bayesian Fit
#
################################################################################

fit1 <- sarprobit(y ~ X-1, W, showProgress=TRUE)
summary(fit1)

#--------MCMC spatial autoregressive probit--------
#Execution time  = 31.788 secs
#
#N draws         =   1000, N omit (burn-in)=    100
#N observations  =    200, K covariates    =      3
## of 0 Y values =     87, # of 1 Y values =    113
#Min rho         = -1.000, Max rho         =  1.000
#--------------------------------------------------
#
#           Estimate Std. Dev  p-level t-value Pr(>|z|)
#Xintercept -0.03436  0.10444  0.38700  -0.329    0.742
#Xx          1.12200  0.14912  0.00000   7.524 1.77e-12 ***
#Xy         -0.91078  0.13358  0.00000  -6.818 1.06e-10 ***
#rho         0.21239  0.19058  0.13900   1.114    0.266
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# use MLE
params <- as.list(c(beta, rho=rho))
fit.mle <- mle(loglik, start=c(0, 1, -1.1, 0.5))
summary(fit.mle)


