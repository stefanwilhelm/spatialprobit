################################################################################
#
# Bayesian estimation of spatial autoregressive probit model (SAR probit)
# using MCMC sampling
#
# Stefan Wilhelm <Stefan.Wilhelm@financial.com>
# using code from Miguel Godinho de Matos <miguelgodinhomatos@cmu.edu>
#
################################################################################

library(tmvtnorm)
library(mvtnorm)
library(Matrix)    # sparseMatrix

# spdep::knearneigh  - k nearest neighbours for spatial weights
# build spatial weight matrix W based on the k nearest neighbors (kNN) (default: k=6)
#
# @param X point coordinates (x, y)
# @param m number of neighbors
# @return sparse matrix (n x n) with nearest neigbors
kNearestNeighbors <- function(x, y, k=6) {
  # number of observations
  n <- length(x)
  
  # spatial weight matrix W based on the k= nearest neighbors
  D <- matrix(NA, n, k)  # (n x k) index matrix to the 6 nearest neigbhors from point i
  for (i in 1:n) {
    px <- x[i]
    py <- y[i]
    # euclidean dist from all points to p
    d <- sqrt((x - px)^2 + (y - py)^2)
    
    # determine the m nearest neighbors (rank 1 is the point itself, ranks 2..(m+1) are the m nearest neighbors
    # TODO: if 2 points have the same distance, e.g. ranks=1, 2, 2.5, 2.5 are odd
    D[i, ] <- which(rank(d) %in% 2:(k+1))
  }
  # sparse matrix representation for spatial weight matrix W
  W <- sparseMatrix(i = rep(1:n, k), j=as.vector(D), x=1/k)
  return(W)
}

# estimated tr(W^i) for i=1,...,100
# see LeSage (2009), chapter 4, p.97/98
#
# "The expectation of the quadratic form u'Au equals tr(A).
#  Since u_i^2 follows a chi^2 distribution with one degree of freedom."
# Pre-calculate traces tr(W^i) i=1,...,100 for the x-impacts calculations
#
# @param W spatial weight matrix (n x n)
# @param o highest order of W^i = W^o
# @param iiter number of MCMC iterations (we run this 50 times)
# @return (n x o) matrix with tr(W^i) in each column, for i=1..,o
tracesWi <- function(W, o=100, iiter=50) {
  n <- nrow(W)
  trW_i <- matrix( data=0, nrow=n, ncol=o )   # n x 100
  for( iii in 1:iiter ){
    u        <- rnorm(n)         #u    -> n x 1      # randn() aus Paket pracma
    wumat    <- matrix(0, n, o)  #wumat-> n x 100
    wumat[,1]<-u             
    wu       <-u                 #wu    -> n x 1
    for(i in 2:o ){
      wu         <- W %*% wu     # (n x n) * (n x 1) = (n x 1) matrix
      wumat[ ,i] <- as.double(wu)
    }
    trW_i        <- trW_i + u * wumat  # n x 100    u' W^i u  (i = 1..100)
  }
  trW_i <- trW_i / iiter
}

# PURPOSE: draw rho from conditional distribution p(rho | beta, z, y)
# ---------------------------------------------------
#  USAGE: 
#
#  results.rmin 
#  results.rmax maximum eigen value
#  results.time execution time
# ---------------------------------------------------
#where epe0 <- t(e0) %*% e0
#      eped <- t(ed) %*% ed
#      epe0d<- t(ed) %*% e0
draw_rho <- function (detval1, detval2, detval1sq, yy, epe0, eped, epe0d, 
    rho, nmk, nrho, lnbprior, u) 
{
    z <- epe0 - 2 * detval1 * epe0d + detval1sq * eped
    z <- -nmk * log(z)
    den <- detval2 + z + lnbprior
    n <- nrho
    adj <- max(den)
    den <- den - adj
    x <- exp(den)
    isum <- sum(yy * (x[2:n] - x[1:(n - 1)])/2)
    z <- abs(x/isum)
    den <- cumsum(z)
    rnd <- u * sum(z)
    ind <- which(den <= rnd)
    idraw <- max(ind)
    if (idraw > 0 && idraw < nrho) {
        results <- detval1[idraw]
    }
    else {
        results <- rho
    }
    return(results)
}

# Bayesian estimation of SAR probit model
#
# @param formula 
sarprobit <- function(formula, W, data, subset, ...) {
  cl <- match.call()                     # cl ist object of class "call"
  mf <- match.call(expand.dots = FALSE)  # mf ist object of class "call"
  m  <- match(c("formula", "data", "subset"), names(mf), 0L)        # m is index vector
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())         # from here mf is a data.frame
  mt <- attr(mf, "terms")                # mt is object of class "terms" and "formula"
  y <- model.response(mf, "numeric")
  if (!is.null(W) && !is.numeric(W) && !inherits(W, "sparseMatrix") && nrow(W) != NROW(y)) 
    stop(gettextf("'W' must be a numeric square matrix, dimension %d should equal %d (number of observations)",
      NROW(W), NROW(y)), domain = NA)
  
  X <- model.matrix(mt, mf, contrasts)
  sar_probit_mcmc(y, X, W, ...)    
}

# Estimate the spatial autoregressive probit model (SAR probit)
# z = rho * W + X \beta + epsilon
# where y = 1 if z >= 0 and y = 0 if z < 0 observable
#
# @param y
# @param X
# @param W spatial weight matrix
# @param ndraw number of MCMC iterations
# @param burn.in  number of MCMC burn-in to be discarded
# @param thinning MCMC thinning factor, defaults to 1
# @param m number of burn.in sampling in inner Gibbs sampling
sar_probit_mcmc <- function(y, X, W, ndraw=1000, burn.in=100, thinning=1, 
  prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12), start=list(rho=0.75, beta=rep(0, ncol(X))),
  m=10, showProgress=FALSE){  

  #start timer
  timet <- Sys.time()
  
  n  <- nrow( X )            # number of observations
  n1 <- nrow( X )          
  n2 <- nrow( W )
  k <- ncol( X )             # number of of parameters/exogenous variables
  I_n <- sparseMatrix(i=1:n, j=1:n, x=1) # sparse identity matrix
  
  
  #validate inputs
  if( length(c(which(y == 0 ),which(y == 1))) != length( y ) ){
    stop('sarp_g: not all y-values are 0 or 1')
  }
  if( n1 != n2 && n1 != n ){
    stop('sarp_g: wrong size of spatial weight matrix W')
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
    stop('sarp_g: intercept term must be in first column of the X-matrix')
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
  
  # Monte Carlo estimation of tr(W^i) for i = 1..o
  trW.i <- tracesWi(W, o=100, iiter=50)
  
  # draw from multivariate normal beta ~ N(c, T). we can precalculate 
  # betadraws ~ N(0, T) befor running the chain and later just create beta as
  # beta = c + betadraws ~ N(c, T)
  betadraws <- rmvnorm(n=(burn.in + ndraw * thinning), mean=rep(0, k), sigma=AA)
  
  # matrices for direct and indirect impacts
  direct       <- matrix(NA, ndraw,p)    # n x p
  indirect     <- matrix(NA, ndraw,p)    # n x p
  total        <- matrix(NA, ndraw,p)    # n x p
  avg_total    <- rep(0,p)               # p x 1
  avg_direct   <- rep(0,p)               # p x 1
  avg_indirect <- rep(0,p)               # p x 1
  
  # just to set a start value for z
  z <- rep(0, n)
  
  for (i in (1 - burn.in):(ndraw * thinning)) {
  
  # 1. sample from z | rho, beta, y using precision matrix H
  # solving equation 
  # (I_n - rho * W) mu = X beta 
  # instead of inverting S = I_n - rho * W as in mu = ( In -  rho W)^{-1} X beta.
  # QR-decomposition for sparse matrices
  QR <- qr(S)  # class "sparseQR"
  mu <- solve(QR, X %*% beta)
  
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
  c <- AA  %*% (tX %*% S %*% z + Tinv %*% c)
  T <- AA   # no update basically on T, TODO: check this
  beta <- as.double(c + betadraws[i + burn.in, ])
  
  # 3. sample from rho | beta, z
  #---- DRAW RHO ----
  #see LeSage 2009 chapter 5 - page 132 for the explanation of the
  #code below which is used for numerical integration of the rho prior.
  #I changed from the original code to match the notation of the book
  #using c0 and cd below instead of b0 and bd ....
  xpz  <- tX %*% z           # X'z
  Wz   <- as.double(W %*% z) # Wz        # SW: coerce Wz to vector 
  # (from n x 1 sparse matrix! we do not need a sparse matrix here)
  xpWz <- tX %*% Wz          # X'Wz      # k x 1
  e0   <-  z - xxpxI %*% xpz  # z  - X(X'X)^-1X' z
  ed   <- Wz - xxpxI %*% xpWz # Wz - X(X'X)^(-1)X'Wz
  epe0 <- as.double(crossprod(e0))  # slightly faster than t(e0) %*% e0
  eped <- as.double(crossprod(ed))
  epe0d<- as.double(crossprod(ed, e0))
  rho  <- draw_rho(detval1,detval2, detval1sq, yy, epe0, eped, epe0d, rho, nmk=nmk, nrho=2001, lnbprior, u=u[i + burn.in])
  
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
  
    # compute effects estimates (direct and indirect impacts) in each MCMC iteration
    o <- 100
    rhovec <- rho^(0:(o-1)) # SW: (100 x 1)   mit [1, rho^1, rho^2 ..., rho^99], see LeSage(2009), eqn (4.145), p.115
    if( cflag == 1 ){ #has intercept
      beff <- beta[-1]      # beff is parameter vector without constant
    }else if(cflag == 0){
      beff <- beta          # no constant in model
    }
  
    # SW: beff is parameter vector without constant!
    # SW: Do we have to compute impacts in every MCMC round or just once at the end?
    # SW: See LeSage (2009), section 5.6.2., p.149/150 for spatial effects estimation in MCMC
    #   direct: M_r(D) = n^{-1} tr(S_r(W))           # SW: efficient approaches available, see chapter 4, pp.114/115
    #    total: M_r(T) = n^{-1} 1'_n S_r(W) 1_n      # SW: Problem: 1'_n S_r(W) 1_n ist dense!
    # indirect: M_r(I) = M_r(T) - M_r(D)
    # SW: See LeSage (2009), section 10.1.6, p.293 for Marginal effects in SAR probit
    pdfz     <- matrix(dnorm(as.double(mu)), ncol=1)  # standard normal pdf phi(mu)
    dir      <- as.real(t(pdfz) %*% trW.i %*% rhovec /n)
    # direct impact : dy_i / d X_ir
    avg_direct     <- dir * beff
      
    #for(r in 1:p ){
    #  tmp               <- apply( dd %*% s * beff[r], 2, sum )  # phi(mu) wird auf die Diagonalelemente draufmultipliziert...
    #  avg_total[r]      <- mean( tmp ) 
    #  total_obs[ ,r]    <- total_obs[,r] + tmp 
    #}
    avg_indirect       <- avg_total - avg_direct    # (r x 1)
    #total[ind, ]      <- avg_total    # an (ndraw-nomit x p) matrix
    direct[ind, ]     <- avg_direct   # an (ndraw-nomit x p) matrix
    #indirect[ind, ]   <- avg_indirect # an (ndraw-nomit x p) matrix
  
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

  #results$predicted <- # prediction required. The default is on the scale of the linear predictors
  class(results)    <- "sarprobit"
  return(results)
}

# compute effects estimates (average direct and indirect impacts) for 
# one set of SAR probit parameters 
# (here we use the point estimates from "sarprobit" object, but we can extend this method to compute impacts for each MCMC iteration)
# Attention: Currently the computation of total and indirect impacts requires a dense (n x n) matrix which will
# work only for small n. If n is large this method will be VERY memory-consuming!
#
#   direct: M_r(D) = n^{-1} tr(S_r(W))           # efficient approaches available, see LeSage (2009), chapter 4, pp.114/115
#    total: M_r(T) = n^{-1} 1'_n S_r(W) 1_n      # Problem: 1'_n S_r(W) 1_n is dense n x n matrix!
# indirect: M_r(I) = M_r(T) - M_r(D)             # Problem: Computation of dense n x n matrix M_r(T) for total effects required!
#  
# @param beta parameter vector
# @param beff parameter vector with constant removed (no intercept parameter)
# @param rho  spatial dependence parameter
# @param W spatial weight matrix
# @param o the maximum o for which to compute tr(W^i) for i = 1..o
# @return
average.impacts <- function(beta, beff=beta, rho, W, X, o=100) {

  p <- length(beff) # number of (non-constant) parameters
  n <- nrow(W) # number of observations
  I_n <- sparseMatrix(i=1:n, j=1:n, x=1) # sparse identity matrix
  
  # solving equation 
  # (I_n - rho * W) mu = X beta 
  # instead of inverting S = I_n - rho * W as in mu = (In -  rho W)^{-1} X beta.
  # QR-decomposition for sparse matrices
  S <- I_n - rho * W
  QR <- qr(S)  # class "sparseQR"
  mu <- solve(QR, X %*% beta)      # n x 1
  
  # S^{-1} = (I_n - rho * W)^{-1}
  SI <- solve(QR)
    
  # Monte Carlo estimation of tr(W^i) for i = 1..o, tr(W^i) is (n x o) matrix
  trW.i <- tracesWi(W, o=o, iiter=50)
  
  # compute effects estimates (direct and indirect impacts) in each MCMC iteration
  rhovec <- rho^(0:(o-1)) # vector (o x 1) with [1, rho^1, rho^2 ..., rho^(o-1)], see LeSage(2009), eqn (4.145), p.115
    
  # See LeSage (2009), section 5.6.2., p.149/150 for spatial effects estimation in MCMC
  #   direct: M_r(D) = n^{-1} tr(S_r(W))           # efficient approaches available, see LeSage (2009), chapter 4, pp.114/115
  #    total: M_r(T) = n^{-1} 1'_n S_r(W) 1_n      # Problem: 1'_n S_r(W) 1_n is dense n x n matrix!
  # indirect: M_r(I) = M_r(T) - M_r(D)             # Problem: Computation of dense n x n matrix M_r(T) for total effects required!
  # See LeSage (2009), section 10.1.6, p.293 for Marginal effects in SAR probit
  pdfz <- matrix(dnorm(as.double(mu)), ncol=1)  # standard normal pdf phi(mu) = phi( (In -  rho W)^{-1} X beta )  # (n x 1)
  dd   <- sparseMatrix(i=1:n, j=1:n, x=pdfz)    # dd is diagonal matrix with pdfz as diagonal (n x n)
  
  dir      <- as.real(t(pdfz) %*% trW.i %*% rhovec /n)  # (1 x n) * (n x o) * (o x 1)
  # direct impact : dy_i / d X_ir = phi((In -  rho W)^{-1} X beta_r) * beta_r ???
  avg_direct     <- dir * beff      # (p x 1)
      
  for(r in 1:p ){
    tmp               <- apply( dd %*% SI * beff[r], 2, sum )  # phi(mu) wird auf die Diagonalelemente draufmultipliziert...
    avg_total[r]      <- mean( tmp ) 
  }
  avg_indirect       <- avg_total - avg_direct    # (p x 1)
   
  return(list(avg_total=avg_total, avg_direct=avg_direct, avg_indirect=avg_indirect))
}

# summary method for class "sarprobit"
summary.sarprobit <- function(object, var_names=NULL, file=NULL, digits = max(3, getOption("digits")-3), ...){
  # check for class "sarprobit"
  if (!inherits(object, "sarprobit")) 
        stop("use only with \"sarprobit\" objects")
        
  nobs      <- object$nobs
  nvar      <- object$nvar
  ndraw     <- object$ndraw
  nomit     <- object$nomit
  draws     <- object$B
  
  #bayesian estimation
  bout_mean <- object$coefficients                         #parameter mean column
  bout_sd   <- apply(draws, 2, sd)                         #parameter sd colum
  # build bayesian significance levels
  # for parameter > 0 count all draws > 0  and determine P(X <= 0)
  # for parameter <= 0 count all draws <0  and determine P(X >= 0)
  bout_sig <- 1 - apply(draws, 2, function(x) { ifelse (mean(x) > 0, sum(x > 0), sum(x < 0)) }) / ndraw
  #standard asymptotic measures
  bout_t    <- bout_mean / bout_sd             #t score b/se
  bout_tPval<- (1 - pt( abs(bout_t), nobs ))*2 #two tailed test = zero probability = z-prob
  #name definition
  if( is.null(var_names)){
    bout_names<- as.matrix(object$names)
  }else{
    bout_names<- as.matrix(var_names)
  }
  
  if(is.null(file)){file <- ""}#output to the console
  #HEADER
  write(sprintf("--------MCMC spatial autoregressive probit--------"), file, append=T)
  #sprintf("Dependent Variable")
  write(sprintf("Execution time  = %6.3f %s", object$time, attr(object$time, "units"))  , file, append=T)
  write(sprintf("N steps for TMVN= %6d"  , object$nsteps), file, append=T)
  write(sprintf("N draws         = %6d, N omit (burn-in)= %6d", ndraw, nomit), file, append=T)
  write(sprintf("N observations  = %6d, K covariates    = %6d", nobs, nvar)  , file, append=T)
  write(sprintf("# fo 0 Y values = %6d, # of 1 Y values = %6d", object$zip, nobs - object$zip) , file, append=T)
  write(sprintf("Min rho         = % 6.3f, Max rho         = % 6.3f", object$rmin, object$rmax), file, append=T)
  write(sprintf("--------------------------------------------------"), file, append=T)
  write(sprintf(""), file, append=T)
  #ESTIMATION RESULTS
  coefficients <- cbind(bout_mean, bout_sd, bout_sig, bout_t, bout_tPval)
  dimnames(coefficients) <- list(bout_names, 
        c("Estimate", "Std. Dev", "Bayes p-level", "t-value", "Pr(>|z|)"))
  printCoefmat(coefficients, digits = digits,
    signif.stars = FALSE)
  if (getOption("show.signif.stars")) {               
    # The solution: using cat() instead of print() and use line breaks
    # cat(paste(strwrap(x, width = 70), collapse = "\\\\\n"), "\n")
    # http://r.789695.n4.nabble.com/Sweave-line-breaks-td2307755.html
    Signif <- symnum(1e-6, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))
    x <- paste("Signif. codes: ", attr(Signif, "legend"), "\n", sep="")
    cat(paste(strwrap(x, width = getOption("width")), collapse = "\\\n"), "\n")
  }
  return(coefficients)
}

# c.sarprobit works in the same way as boot:::c.boot().
c.sarprobit <- function(...) {
 args <- list(...)
 nm <- lapply(args, names)
 if (!all(sapply(nm, function(x) identical(x, nm[[1]]))))
   stop("arguments are not all the same type of \"sarprobit\" object")
 res <- args[[1]]
 res$time  <- as.difftime(max(as.numeric(sapply(args, "[[", "time"), units="secs")), units="secs") 
 res$ndraw <- sum(sapply(args, "[[", "ndraw")) 
 res$nomit <- sum(sapply(args, "[[", "nomit"))
 res$B     <- do.call(rbind, lapply(args, "[[", "B"))
 res$bdraw <- do.call(rbind, lapply(args, "[[", "bdraw"))
 res$pdraw <- do.call(rbind, lapply(args, "[[", "pdraw"))
 res$coefficients <- colMeans(res$B) 
 res$beta  <- res$coefficients[1:res$nvar]
 res$rho   <- res$coefficients[res$nvar+1]
 res
}


# extract the coefficients
coef.sarprobit <- function(object, ...) {
 if (!inherits(object, "sarprobit")) 
        stop("use only with \"sarprobit\" objects")
 return(object$coefficients)
}

# extract the coefficients
coefficients.sarprobit <- function(object, ...) {
 UseMethod("coef", object)
}

# plot MCMC results for class "sarprobit" (draft version);
# diagnostic plots for results (trace plots, ACF, posterior density function)
# method is very similar to plot.lm()
#
# @param x
# @param which
# @param ask
# @param trueparam a vector of "true" parameter values to be marked in posterior density plot
plot.sarprobit <- function(x, which=c(1, 2, 3), 
  ask = prod(par("mfcol")) < length(which) && dev.interactive(), ..., trueparam=NULL) {
 if (!inherits(x, "sarprobit")) 
        stop("use only with \"sarprobit\" objects")
 if (!is.numeric(which) || any(which < 1) || any(which > 3)) 
        stop("'which' must be in 1:3")
        
 names <- x$names
 B <- x$B
 k <- ncol(B)
  
 show <- rep(FALSE, 3)
 show[which] <- TRUE
 
 if (ask) {
   oask <- devAskNewPage(TRUE)
   on.exit(devAskNewPage(oask))
 }
 if (show[1L]) {
  # trace plots
  for (i in 1:k) {
    plot(1:nrow(B), B[,i], type="l", xlab="iteration", ylab=names[i], main=substitute("Trace plot of "*x, list(x=names[i])), ...)
    if (!is.null(trueparam)) abline(h=trueparam[i], col="red", lty=2)
  }
 }

 if (show[2L]) {
   # ACFs
   for (i in 1:k) {
     acf(B[,i], main=substitute("ACF of "*x, list(x=names[i])), ...)
   }
 }
 
 if (show[3L]) {
   # posterior distribution
   for (i in 1:k) {
     plot(density(B[,i]), main=substitute("Posterior distribution of "*x, list(x=names[i])), ...)
     if (!is.null(trueparam)) abline(v=trueparam[i], col="red", lty=2)
   }
 }
}

# return fitted values of SAR probit (on reponse scale vs. linear predictor scale)
fitted.sarprobit <- function(object, ...) {
  object$fitted.value
}

predict.glm <- function (object, newdata = NULL, 
    type = c("link", "response", "terms"), 
    se.fit = FALSE, dispersion = NULL, terms = NULL, 
    na.action = na.pass, ...) 
{
    type <- match.arg(type)
    na.act <- object$na.action
    object$na.action <- NULL
    if (!se.fit) {
        if (missing(newdata)) {
            pred <- switch(type, 
              link = object$linear.predictors, 
              response = object$fitted.values, 
              terms = predict.lm(object, se.fit = se.fit, scale = 1, type = "terms", terms = terms)
            )
            if (!is.null(na.act)) 
                pred <- napredict(na.act, pred)
        }
        else {
            pred <- predict.lm(object, newdata, se.fit, scale = 1, 
                type = ifelse(type == "link", "response", type), 
                terms = terms, na.action = na.action)
            switch(type, response = {
                pred <- family(object)$linkinv(pred)
            }, link = , terms = )
        }
    }
    else {
        if (inherits(object, "survreg")) 
            dispersion <- 1
        if (is.null(dispersion) || dispersion == 0) 
            dispersion <- summary(object, dispersion = dispersion)$dispersion
        residual.scale <- as.vector(sqrt(dispersion))
        pred <- predict.lm(object, newdata, se.fit, scale = residual.scale, 
            type = ifelse(type == "link", "response", type), 
            terms = terms, na.action = na.action)
        fit <- pred$fit
        se.fit <- pred$se.fit
        switch(type, response = {
            se.fit <- se.fit * abs(family(object)$mu.eta(fit))
            fit <- family(object)$linkinv(fit)
        }, link = , terms = )
        if (missing(newdata) && !is.null(na.act)) {
            fit <- napredict(na.act, fit)
            se.fit <- napredict(na.act, se.fit)
        }
        pred <- list(fit = fit, se.fit = se.fit, residual.scale = residual.scale)
    }
    pred
}


