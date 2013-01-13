# Estimating a Probit Model with Spatial Errors (SEM Probit)
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

if (FALSE) {
 library(tmvtnorm)
 library(Matrix)    # sparseMatrix

 setwd("G:/R/spatialprobit/R")
 source("sar_base.r")
 source("matrix_operations.r")
 source("stats_distributions.r")
 source("utility_functions.r")
 setwd("G:/R/spatialprobit/Misc")
 source("Metropolis-Hastings-rho.R")
  
}

# Bayesian estimation of the probit model with spatial errors (SEM probit)
#
# @param formula 
semprobit <- function(formula, W, data, subset, ...) {
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
  sem_probit_mcmc(y, X, W, ...)    
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
# @param showProgress
sem_probit_mcmc <- function(y, X, W, ndraw=1000, burn.in=100, thinning=1, 
  prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12), start=list(rho=0.75, beta=rep(0, ncol(X))),
  m=10, showProgress=FALSE){  

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
  mu <- X %*% beta
  
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
  
  # 3. sample from rho | beta, z using Metropolis-Hastings with burn.in=20
  rho  <- draw_rho_metropolis(type="SEM", n=1, z, I_n, W, X, burn.in=20, start.value=rho, c=1)$rho_t[20+1]
  
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
  fitted.values   <- X %*% beta                     # E[z | beta] = (X * beta)
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
  #results$total     <- total
  #results$direct    <- direct
  #results$indirect  <- indirect
  results$W <- W
  results$X <- X

  #results$predicted <- # prediction required. The default is on the scale of the linear predictors
  class(results)    <- "semprobit"
  return(results)
}

# extract the coefficients
coef.semprobit <- function(object, ...) {
 if (!inherits(object, "semprobit")) 
        stop("use only with \"semprobit\" objects")
 return(object$coefficients)
}

# extract the coefficients
coefficients.semprobit <- function(object, ...) {
 UseMethod("coef", object)
}


# define summary method
# summary method for class "semprobit"
summary.semprobit <- function(object, var_names=NULL, file=NULL, digits = max(3, getOption("digits")-3), ...){
  # check for class "sarprobit"
  if (!inherits(object, "semprobit")) 
        stop("use only with \"semprobit\" objects")
        
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
  write(sprintf("# of 0 Y values = %6d, # of 1 Y values = %6d", object$zip, nobs - object$zip) , file, append=T)
  write(sprintf("Min rho         = % 6.3f, Max rho         = % 6.3f", object$rmin, object$rmax), file, append=T)
  write(sprintf("--------------------------------------------------"), file, append=T)
  write(sprintf(""), file, append=T)
  #ESTIMATION RESULTS
  coefficients <- cbind(bout_mean, bout_sd, bout_sig, bout_t, bout_tPval)
  dimnames(coefficients) <- list(bout_names, 
        c("Estimate", "Std. Dev", "Bayes p-level", "t-value", "Pr(>|z|)"))
  printCoefmat(coefficients, digits = digits,
    signif.stars = getOption("show.signif.stars"))      
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
  return(invisible(coefficients))
} 

# plot MCMC results for class "semprobit" (draft version);
# diagnostic plots for results (trace plots, ACF, posterior density function)
# method is very similar to plot.lm()
#
# @param x
# @param which
# @param ask
# @param trueparam a vector of "true" parameter values to be marked in posterior density plot
plot.semprobit <- function(x, which=c(1, 2, 3), 
  ask = prod(par("mfcol")) < length(which) && dev.interactive(), ..., trueparam=NULL) {
 if (!inherits(x, "semprobit")) 
        stop("use only with \"semprobit\" objects")
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

if (FALSE) {
# example:
library(tmvtnorm)
library(Matrix)

n <- d <- 300
m <- 3
W <- sparseMatrix(i=rep(1:d, each=m), 
  j=replicate(d, sample(x=1:d, size=m, replace=FALSE)), x=1/m, dims=c(d, d))
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
eps <- rnorm(n=n, mean=0, sd=1)   # Normierung!!!
rho <- 0.75
X   <- cbind(x1=1, x2=runif(n=n, -2, 2))
beta <- c(-0.2, 0.5)

z <- as.vector(X %*% beta + solve(I_n - rho * W) %*% eps)     # SEM model
y <- as.numeric(z >= 0)
mu <- X %*% beta

# truncation points for z, depend only on y, can be precalculated
lower <- ifelse(y > 0, 0,  -Inf)
upper <- ifelse(y > 0, Inf,   0)

S <- I_n - rho * W
H <- t(S) %*% S            # precision matrix H for beta | rho, z, y

znew <- as.double(rtmvnorm.sparseMatrix(n=1, mean=mu, H=H, lower=lower, upper=upper, burn.in=10))

# stürzt in rtmvnorm.sparseMatrix ab, wenn I_n/W dense "matrix" sind.
# dann ist auch H eine dense matrix. TODO: Fehler abfangen!
Rprof("SEM-probit.out")
fit <- sem_probit_mcmc(y, X, W, ndraw=500, burn.in=100, thinning=2, showProgress=TRUE)
Rprof(NULL)
summaryRprof("SEM-probit.out")
summary(fit)
plot(fit)
}


