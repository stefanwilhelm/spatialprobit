# Estimating a Probit Model with Spatial Errors (SEMProbit)
# autocorrelation in the error rather than on a lag variable
#
# Referenzen/Anwendungen von Spatial Probit
#
# (1) LeSage (2011)
# (2) Coughlin (2003) Spatial Probit and Geographic Patters (Federal Reserve Bank): http://ideas.repec.org/p/fip/fedlwp/2003-042.html
# (3) Spatial Probit estimation of Freetrade agreements: http://ideas.repec.org/p/gii/giihei/heidwp07-2010.html
# (4) Probit with spatial correlation: http://www.jstor.org/stable/1400629
# http://www.sciencedirect.com/science?_ob=ArticleListURL&_method=list&_ArticleListID=1879027857&_st=13&view=c&_acct=C000026389&_version=1&_urlVersion=0&_userid=525223&md5=0b019d20f88711fe7bd9e24eb4211d66&searchtype=a

#@UNPUBLISHED{Coughlin2003,
#  author = {Cletus C. Coughlin and Thomas A. Garrett and Rubén Hernández-Murillo},
#  title = {Spatial probit and the geographic patterns of state lotteries},
#  note = {Federal Reserve Bank of St. Louis Working Paper 2003-042},
#  year = {2003},
#  owner = {stefan},
#  timestamp = {2013.01.06}
#}
#
#@ARTICLE{Jaimovich2012,
#  author = {Dany Jaimovich},
#  title = {A Bayesian spatial probit estimation of Free Trade Agreement contagion},
#  journal = {Applied Economics Letters},
#  year = {2012},
#  volume = {19},
#  pages = {579-583},
#  number = {6},
#  owner = {stefan},
#  timestamp = {2013.01.06}
#}

# Spatial Errors

# http://www.jstor.org/stable/1400629
#@ARTICLE{Marsh2000,
#  author = {Thomas L. Marsh and Ron C. Mittelhammer and Ray G. Huffaker},
#  title = {Probit with Spatial Correlation by Field Plot: Potato Leafroll Virus
#	Net Necrosis in Potatoes},
#  journal = {Journal of Agricultural, Biological, and Environmental Statistics},
#  year = {2000},
#  volume = {5},
#  pages = {22-36},
#  owner = {stefan},
#  timestamp = {2013.01.06}
#}

#Miguel Godinho de Matos
# A simple extension to this model would be to consider autocorrelation 
# in the error rather than on a lag variable.
# z = xB + u
# u = pWu + e
#
# This formulation allows for a utility based interpretation. 
# Also it requires the truncated multivariate normal for estimation …. as well:
#  z = xB + (I - pW)^(-1)e 
# We should probably add this to the model since it also very often used.
#
# Im SAR Probit ist der DGP:
# z = rho * W * z + X beta + e
# Sz = X beta + e  where S = (I_n - rho W)
# z = (I - pW)^(-1)xB + (I - pW)^(-1)e
# und
# z ~ N((I - pW)^(-1) xB, [(I - pW)'(I - pW)]^(-1))
#
# D.h. für das SEM Probit Model
# z = xB + (I - pW)^(-1)e sollte doch das Ganze
# so aussehen:
# z ~ N(xB, [(I - pW)'(I - pW)]^(-1))
#
# Ausserdem ändert sich die Verteilung für beta p(beta | rho,z)
# und für rho p(rho | beta, z)

# SEM Probit: 
# 1. Erwartungswert von Truncated Multivariate Normal p(z | beta, rho) ändert sich
#    z ~ TN(X beta, H)
# 2. Erwartungswert von von Normalverteilung p(beta | z, rho) ändert sich
# 3. Verteilung von Rho:   p(rho | z, beta) ???
#     SAR Probit:  p(rho | z, beta) ~  |I_n - rho W| exp(-1/2*(S z - X beta)'(S z - X beta))
#                  Sz - Xbeta = e (e ~ N(i.i.d))
#     SEM Probit:  p(rho | z, beta) ~  |I_n - rho W| exp(-1/2*(S z - S X beta)'(S z - S X beta)) ???
#                  z = Xbeta + S^(-1)e
#                  Sz - S X beta = e
# TO BE VERIFIED
# 4. Marginal Effects ändern sich gegenüber SAR Probit: SEM Probit hat kein Spatial Spillover, Marginal Effects wie beim Probit

library(tmvtnorm)
library(Matrix)    # sparseMatrix

if (FALSE) {
 setwd("F:/R/spatialprobit/R")
 source("sar_base.r")
 source("matrix_operations.r")
 source("stats_distributions.r")
 source("utility_functions.r") 
}


# Estimate the probit model with spatial errors (SEM probit)
# z = xB + u
# u = pWu + e
#
# Also it requires the truncated multivariate normal for estimation …. as well:
#  z = xB + (I - pW)^(-1)e 
# where y = 1 if z >= 0 and y = 0 if z < 0 observable
#
# @param y
# @param X
# @param W spatial weight matrix
# @param ndraw number of MCMC iterations
# @param burn.in  number of MCMC burn-in to be discarded
# @param thinning MCMC thinning factor, defaults to 1
# @param m number of burn.in sampling in inner Gibbs sampling
# @param prior
# @param start
# @param m
# @param computeMarginalEffects
# @param showProgress
sem_probit_mcmc <- function(y, X, W, ndraw=1000, burn.in=100, thinning=1, 
  prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12), start=list(rho=0.75, beta=rep(0, ncol(X))),
  m=10, computeMarginalEffects=FALSE, showProgress=FALSE){  

  #start timer
  timet <- Sys.time()
  
  n  <- nrow( X )            # number of observations
  n1 <- nrow( X )          
  n2 <- nrow( W )
  k <- ncol( X )             # number of of parameters/exogenous variables
  I_n <- sparseMatrix(i=1:n, j=1:n, x=1)  # sparse identity matrix
  
  
  #validate inputs
  if( length(c(which(y == 0 ),which(y == 1))) != length( y ) ){
    stop('semprobit: not all y-values are 0 or 1')
  }
  if( n1 != n2 && n1 != n ){
    stop('semprobit: wrong size of spatial weight matrix W')
  }
  # check if we have a constant term in X
  ind <- match( n, apply(X,2,sum))
  if( is.na(ind) ){
    cflag <- 0
    p     <- k
  }else if( ind == 1 ){
    cflag <- 1
    p     <- k - 1
  }else{
    stop('semprobit: intercept term must be in first column of the X-matrix')
  }
  
  # MCMC sampling of beta
  rho  <- start$rho          # start value of row
  beta <- start$beta         # start value of parameters, prior value, we could also sample from beta ~ N(c, T)
  
  # conjugate prior beta ~ N(c, T)
  # parametrize, default to diffuse prior, for beta, e.g. T <- diag(k) * 1e12
  c <- rep(0, k)             # prior distribution of beta ~ N(c, T) : c = 0
  if (is.matrix(prior$T) && ncol(prior$T) == k && isSymmetric(prior$T) && det(prior$T) > 0) {
    T <- prior$T               # prior distribution of beta ~ N(c, T) : T = I_n --> diffuse prior
  } else {
    T <- diag(k)*1e12
  }
  
  Tinv <- solve(T)           # T^{-1}
  S <- I_n - rho * W
  H <- t(S) %*% S            # precision matrix H for beta | rho, z, y
  
  # truncation points for z, depend only on y, can be precalculated
  lower <- ifelse(y > 0, 0,  -Inf)
  upper <- ifelse(y > 0, Inf,   0)
  
  rmin       <- -1   # use -1,1 rho interval as default
  rmax       <-  1
  ldetflag   <-  0   # default to 1999 Pace and Barry MC determinant approx
  tmp <- sar_lndet(ldetflag, W, rmin, rmax)
  detval <- tmp$detval
  
  # Some precalculated quantities for drawing rho
  # rho ~ Beta(a1, a2) prior
  a1         <-  1.0
  a2         <-  1.0
  lnbprior <- log(beta_prior(detval[,1],a1,a2))
  u        <- runif(thinning * ndraw + burn.in)   # u ~ U(0, 1)
  nrho     <- 2001
  nmk      <- (n-k)/2
  detval1  <- detval[,1]  # SW: avoid multiple accesses to detval[,1]
  detval2  <- detval[,2]
  detval1sq <- detval1 * detval1
  yy        <- (detval1[2:nrho] + detval1[1:(nrho-1)])
    
  # matrix to store the beta + rho parameters for each iteration/draw
  B <- matrix(NA, ndraw, k+1)
  
  # progress bar
  if (showProgress) {
    pb <- txtProgressBar(min=0, max=(thinning * ndraw + burn.in), initial=0, style=3)
  }
  
  # immutable matrices
  tX <- t(X)                       # X'               # k x n
  xpx  <- t(X) %*% X               # (X'X)            # k x k
  xpxI <- solve(xpx)               # (X'X)^{-1}       # k x k
  #xxpxIxp <- X %*% xpxI %*% tX    # X(X'X)^(-1)X'    # n x n (argh!)
  xxpxI    <- X %*% xpxI           # X(X'X)^(-1)     # n x k (better, compromise)
  AA       <- solve(xpx + Tinv)    # (X'X + T^{-1})^{-1}
  
  # draw from multivariate normal beta ~ N(c, T). we can precalculate 
  # betadraws ~ N(0, T) before running the chain and later just create beta as
  # beta = c + betadraws ~ N(c, T)
  betadraws <- rmvnorm(n=(burn.in + ndraw * thinning), mean=rep(0, k), sigma=AA)
  
  # names of non-constant parameters
  if(cflag == 0) {
    namesNonConstantParams <- colnames(X)
  } else {
    namesNonConstantParams <- colnames(X)[-1]
  }
    
  # just to set a start value for z
  z <- rep(0, n)
  
  for (i in (1 - burn.in):(ndraw * thinning)) {
  
  # 1. sample from z | rho, beta, y using precision matrix H
  # z ~ N(X beta, H)
  mu <- rep(X %*% beta, n)
  
  # see LeSage (2009) for choice of burn-in size, often m=5 or m=10 is used!
  # we can also use m=1 together with start.value=z, see LeSage (2009), section 10.1.5
  if (m==1) {
    z <- as.double(rtmvnorm.sparseMatrix(n=1, mean=mu, H=H, 
      lower=lower, upper=upper, burn.in=m, start.value=z))
  } else {
    z <- as.double(rtmvnorm.sparseMatrix(n=1, mean=mu, H=H, 
      lower=lower, upper=upper, burn.in=m))
  }
    
  # 2. sample from beta | rho, z, y
  # SAR Probit: c <- AA  %*% (tX %*% S %*% z + Tinv %*% c)
  # SEM Probit: c <- AA  %*% (tX %*% z + Tinv %*% c)
  c <- AA  %*% (tX %*% z + Tinv %*% c)
  T <- AA   # no update basically on T, TODO: check this
  beta <- as.double(c + betadraws[i + burn.in, ])
  
  # 3. sample from rho | beta, z
  #---- DRAW RHO ----
  #see LeSage 2009 chapter 5 - page 132 for the explanation of the
  #code below which is used for numerical integration of the rho prior.
  #I changed from the original code to match the notation of the book
  #using c0 and cd below instead of b0 and bd ....
  #xpz  <- tX %*% z           # X'z
  #Wz   <- as.double(W %*% z) # Wz        # SW: coerce Wz to vector 
  ## (from n x 1 sparse matrix! we do not need a sparse matrix here)
  #xpWz <- tX %*% Wz          # X'Wz      # k x 1
  #e0   <-  z - xxpxI %*% xpz  # z  - X(X'X)^-1X' z
  #ed   <- Wz - xxpxI %*% xpWz # Wz - X(X'X)^(-1)X'Wz
  #epe0 <- as.double(crossprod(e0))  # slightly faster than t(e0) %*% e0
  #eped <- as.double(crossprod(ed))
  #epe0d<- as.double(crossprod(ed, e0))
  #rho  <- draw_rho(detval1,detval2, detval1sq, yy, epe0, eped, epe0d, rho, nmk=nmk, nrho=2001, lnbprior, u=u[i + burn.in])
  rho  <- draw_rho_metropolis(type="SEM", n=1, z, W, X, burn.in=100, start.value=rho, c=1)$rho_t[1]
  
  ############################################################################## 
  
  # update S and H
  S <- I_n - rho * W
  H <- t(S) %*% S      # H = S'S  / SW: crossprod(S) does not seem to work!

  if (i > 0) {
    if (thinning == 1) {
      ind <- i
    }
    else if (i%%thinning == 0) {
      ind <- i%/%thinning
    } else {
      next
    }
    
    B[ind,] <- c(beta, rho)
  
  ##############################################################################
    
  }
    
  if (showProgress) setTxtProgressBar(pb, i + burn.in) # update progress bar
  }
    
  if (showProgress)  close(pb) #close progress bar
  
  # fitted values for estimates (based on z rather than binary y like in fitted(glm.fit))
  # (on reponse scale y vs. linear predictor scale z...)
  beta  <- colMeans(B)[1:k]
  rho   <- colMeans(B)[k+1]
  S     <- (I_n - rho * W)
  fitted.values   <- solve(qr(S), X %*% beta)   # z = (I_n - rho * W)^{-1}(X * beta)
  fitted.response <- as.numeric(fitted.values >= 0) 
  # TODO: linear.predictors  vs. fitted.values
  
  # result
  results       <- NULL
  results$time  <- Sys.time() - timet
  results$nobs  <- n          # number of observations
  results$nvar  <- k          # number of explanatory variables
  results$y     <- y 
  results$zip   <- n - sum(y) # number of zero values in the y-vector
  results$beta  <- colMeans(B)[1:k]
  results$rho   <- colMeans(B)[k+1]
  results$coefficients <- colMeans(B)
  results$fitted.values <- fitted.values
  #results$fitted.reponse <- fitted.reponse  # fitted values on reponse scale (binary y variable)
  results$ndraw <- ndraw
  results$nomit <- burn.in
  results$a1        <- a1
  results$a2        <- a2
  results$rmax      <- rmax 
  results$rmin      <- rmin
  results$tflag     <- 'plevel'
  results$lflag     <- ldetflag
  results$cflag     <- cflag
  results$lndet     <- detval
  results$names     <- c(colnames(X), 'rho')
  results$B         <- B        # (beta, rho) draws
  results$bdraw     <- B[,1:k]  # beta draws
  results$pdraw     <- B[,k+1]  # rho draws
  results$total     <- total
  results$direct    <- direct
  results$indirect  <- indirect
  results$W <- W
  results$X <- X

  #results$predicted <- # prediction required. The default is on the scale of the linear predictors
  class(results)    <- "semprobit"
  return(results)
}

# example:
n <- d <- 50
m <- 3
W <- sparseMatrix(i=rep(1:d, each=m), 
  j=replicate(d, sample(x=1:d, size=m, replace=FALSE)), x=1/m, dims=c(d, d))
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
#I_n <- Diagonal(n)
W <- as.matrix(W)
I_n <- as.matrix(I_n)
eps <- rnorm(n=n, mean=0, sd=0.1)
rho <- 0.75
X   <- cbind(1, runif(n=n, -2, 2))
beta <- c(1, 1)

z <- as.vector(X %*% beta + solve(I_n - rho * W) %*% eps)     # SEM model
y <- as.numeric(z >= 0)

fit <- sem_probit_mcmc(y, X, W, ndraw=100, burn.in=100, thinning=1, showProgress=TRUE)


