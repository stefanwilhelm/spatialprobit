require(pracma)
require(truncnorm)
require(msm)
require(Matrix)

#source("sar_base.r")

# PURPOSE: Bayesian estimates of the spatial autoregressive probit model
#          y = rho*W*y + XB + e, e = N(0,I_n)
#          y is a binary 0,1 nx1 vector
#          B = N(c,T), 
#          1/sige = Gamma(nu,d0), 
#          rho = Uniform(rmin,rmax), or rho = beta(a1,a2); 
#-------------------------------------------------------------
# USAGE: results = sarp_g(y,x,W,ndraw,nomit,prior)
# where: y = dependent variable vector (nobs x 1)
#        x = independent variables matrix (nobs x nvar), 
#            the intercept term (if present) must be in the first column of the matrix x
#        W = spatial weight matrix (standardized, row-sums = 1)
#    ndraw = # of draws
#    nomit = # of initial draws omitted for burn-in            
#    prior = a structure variable with:
#            prior.nsteps = # of samples used by truncated normal Gibbs sampler
#            prior.beta  = prior means for beta,   c above (default 0)
#            priov.bcov  = prior beta covariance , T above (default 1e+12)
#            prior.a1    = parameter for beta(a1,a2) prior on rho see: 'help beta_prior'
#            prior.a2    = (default = 1.0, a uniform prior on rmin,rmax) 
#            prior.eig   = 0 for default rmin = -1,rmax = +1, 1 for eigenvalue calculation of these
#            prior.rmin  = (optional) min rho used in sampling (default = -1)
#            prior.rmax  = (optional) max rho used in sampling (default = 1)  
#            prior.lflag = 0 for full lndet computation (default = 1, fastest)
#                        = 1 for MC approx (fast for large problems)
#                        = 2 for Spline approx (medium speed)
#            prior.order = order to use with prior.lflag = 1 option (default = 50)
#            prior.iter  = iters to use with prior.lflag = 1 option (default = 30) 
#            prior.lndet = a matrix returned by sar, sar_g, sarp_g, etc.
#                          containing log-determinant information to save time
#-------------------------------------------------------------
# RETURNS:  a structure:
#          results.meth     = 'sarp_g'
#          results.beta     = posterior mean of bhat based on draws
#          results.rho      = posterior mean of rho based on draws
#          results.sige     = posterior mean of sige based on draws
#          results.sigma    = posterior mean of sige based on (e'*e)/(n-k)
#          results.bdraw    = bhat draws (ndraw-nomit x nvar)
#          results.pdraw    = rho  draws (ndraw-nomit x 1)
#          results.sdraw    = sige draws (ndraw-nomit x 1)
#          results.total    = a matrix (ndraw,nvars-1) total x-impacts
#          results.direct   = a matrix (ndraw,nvars-1) direct x-impacts
#          results.indirect = a matrix (ndraw,nvars-1) indirect x-impacts
#          results.total_obs= a matrix (ndraw,nvars-1) observation-level total x-impacts
#          results.vmean  = mean of vi draws (nobs x 1) 
#          results.rdraw  = r draws (ndraw-nomit x 1) (if m,k input)
#          results.bmean  = b prior means, prior.beta from input
#          results.bstd   = b prior std deviations sqrt(diag(prior.bcov))
#          results.novi   = 1 for prior.novi = 1, 0 for prior.rval input
#          results.nobs   = # of observations
#          results.nvar   = # of variables in x-matrix
#          results.ndraw  = # of draws
#          results.nomit  = # of initial draws omitted
#          results.nsteps = # of samples used by Gibbs sampler for TMVN
#          results.y      = y-vector from input (nobs x 1)
#          results.zip    = # of zero y-values
#          results.yhat   = mean of posterior predicted (nobs x 1)
#          results.resid  = residuals, based on posterior means
#          results.rsqr   = r-squared based on posterior means
#          results.rbar   = adjusted r-squared
#          results.a1     = a1 parameter for beta prior on rho from input, or default value
#          results.a2     = a2 parameter for beta prior on rho from input, or default value
#          results.time1  = time for eigenvalue calculation
#          results.time2  = time for log determinant calcluation
#          results.time3  = time for sampling
#          results.time   = total time taken  
#          results.rmax   = 1/max eigenvalue of W (or rmax if input)
#          results.rmin   = 1/min eigenvalue of W (or rmin if input)          
#          results.tflag  = 'plevel' (default) for printing p-levels
#                         = 'tstat' for printing bogus t-statistics 
#          results.lflag  = lflag from input
#          results.cflag  = 1 for intercept term, 0 for no intercept term
#          results.iter   = prior.iter option from input
#          results.order  = prior.order option from input
#          results.limit  = matrix of [rho lower95,logdet approx, upper95] 
#                           intervals for the case of lflag = 1
#           results.lndet = a matrix containing log-determinant information
#                           (for use in later function calls to save time)
#           results.mlike = log marginal likelihood (a vector ranging over
#                           rho values that can be integrated for model comparison)
# --------------------------------------------------------------
# NOTES: - the intercept term (if you have one)
#          must be in the first column of the matrix x
# --------------------------------------------------------------
# SEE ALSO: (sarp_gd, sarp_gd2 demos) prt
# --------------------------------------------------------------
# REFERENCES: LeSage and Pace (2009) Chapter 10 on Bayesian estimation 
#             of spatial probit regression models.
# For lndet information see: Chapter 4 
#----------------------------------------------------------------

# written by:
# James P. LeSage, last updated 3/2010
# Dept of Finance & Economics
# Texas State University-San Marcos
# 601 University Drive
# San Marcos, TX 78666
# jlesage@spatial-econometrics.com

sar_probit_mcmc.miguel <- function(y,x,W,ndraw,nomit,prior){
  n          <- length( y )
  n1         <- nrow( x )
  n2         <- nrow( W )
  k          <- ncol( x )
  yin        <- y
  In         <- speye( n )
  
  eflag      <-  0   # default to not computing eigenvalues
  ldetflag   <-  0   # default to 1999 Pace and Barry MC determinant approx
  order      <- 50   # there are parameters used by the MC det approx
  iter       <- 30   # defaults based on Pace and Barry recommendation
  rmin       <- -1   # use -1,1 rho interval as default
  rmax       <-  1  
  detval     <-  0   # just a flag
  rho        <-  0.5 
  a1         <-  1.0 
  a2         <-  1.0 
  c          <-  zeros(k,1);   # diffuse prior for beta
  Ts         <-  speye(k)*1e+12
  prior_beta <- 0   # flag for diffuse prior on beta
  novi_flag  <- 0   # do vi-estimates
  inform_flag<- 0
  metflag    <- 0
  nsample    <- 5  #same as in the demo
  
  #validate inputs
  if( length(c(which(y == 0 ),which(y == 1))) != length( y ) ){
    print('sarp_g: not all y-values are 0 or 1')
    return(NULL)
  }
  if( n1 != n2 && n1 != n ){
    print('sarp_g: wrong size weight matrix W')
    return(NULL)
  }
  ind <- match( n, apply(x,2,sum))
  if( is.na(ind) ){
    cflag <- 0
    p     <- k
  }else if( ind == 1 ){
    cflag <- 1
    p     <- k - 1
  }else{
    print('sarp_g: intercept term must be in first column of the x-matrix')
    return(NULL)
  }
  
  results       <- NULL
  results$nobs  <- n
  results$nvar  <- k
  results$y     <- y 
  results$zip   <- n - sum(y) # number of zero values in the y-vector
  results$cflag <- cflag
  results$p     <- p
  results$order <- order  
  results$iter  <- iter

  #start timer
  timet <- Sys.time() 
  time1 <- 0; time2 <- 0; time3 <- 0

  #compute eigen values of spatial weight matrix W (rmin, rmax)
  tmp   <- sar_eigs(eflag,W)
  rmin  <- tmp$rmin
  rmax  <- tmp$rmax
  time1 <- tmp$time
  rm(tmp);gc()
  
  results$time1 <- time1
  
  tmp    <- sar_lndet(ldetflag, W, rmin, rmax)
  detval <- tmp$detval
  time2  <- tmp$time
  rm(tmp);gc()
  
  results$time2 <- time2

  # pre-calculate traces for the x-impacts calculations
  # TODO: consider converting this to a suport function
  iiter     <- 50
  o         <- 100
  diag_ests <- Matrix( data=0, nrow=n, ncol=o )
  for( iii in 1:iiter ){
    u        <- randn(n,1)       #u    -> nx1    
    umat     <- u[, ones(1,o)]   #umat -> nx100  
    wumat    <- zeros(n,o)       #wumat-> nx100
    wumat[,1]<-u[,1]             
    wu       <-u[,1]            #wu    -> nx1
    for(ii in 2:o ){
      wu         <- W %*% wu     #nx1 matrix
      wumat[ ,ii] <- wu[,1]
    }
    diag_estimates_iii <- umat * wumat
    diag_ests          <- diag_ests + diag_estimates_iii
  }
  estimated_diags <- diag_ests / iiter
  
  # storage for draws
  bsave        <- zeros(ndraw-nomit,k);
  psave        <- zeros(ndraw-nomit,1);
  ymean        <- zeros(n,1);
  acc_rate     <- zeros(ndraw,1);
  total        <- zeros(ndraw-nomit,p);
  total_obs    <- zeros(n,p);
  direct       <- zeros(ndraw-nomit,p);
  indirect     <- zeros(ndraw-nomit,p);
  avg_total    <- zeros(p,1);
  avg_direct   <- zeros(p,1);
  avg_indirect <- zeros(p,1);

  # ====== initializations
  print('sarp_g: MCMC estimation ...')
  pb<-txtProgressBar(min=0,max=ndraw, initial=0, style=3)
  # compute this stuff once to save time
  sige <- 1
  #IMMUTABLE MATRICES
  TI   <- qr.solve(Ts)
  TIc  <- TI %*% c
  Wadd <-  W + t(W)         # W + W'
  WpW  <- t(W) %*% W        # W'W
  xp   <- t(x)              # X'
  xpx  <- xp %*% x          # X'X
  xpxI <- qr.solve( xpx )   #(X'X)^{-1}
  #MUTABLE MATRICES
  xpy  <- xp %*% y          # X'y
  Wy   <-  W %*% y          # W'y
  xpWy <- xp %*% Wy         # X'Wy
  
  dd   <- Matrix( data=0, nrow=n,ncol=n,sparse=T)
  z    <- zeros(n,1)
  
  iter <- 1
  while( iter <= ndraw ){
    #---- DRAW BETA ----
    A    <- xpx + sige * TI
    AI   <- qr.solve( A )               #T* = ( X'X + T^(-1) )^-1
    ys   <- y - rho*Wy                  
    b    <- xp %*% ys + sige * TIc      #X'(y - rho*Wy) + T^(-1)
    b0   <- AI %*% b                    #C* = ( X'X + T^(-1) )^-1 * X'(y - rho*Wy) + T^(-1)
    bhat <- norm_rnd( sige * AI  ) + b0 #Draw from N(0,T*)
    xb   <- x %*% bhat
    
    #---- DRAW RHO ----
    #see LeSage 2009 chapter 5 - page 132 for the explanation of the
    #code below which is used for numerical integration of the rho prior.
    #I changed from the original code to match the notation of the book
    #using c0 and cd below instead of b0 and bd ....
    #c0   <- mldivide( as.matrix(xpx) , as.matrix(xpy) ) #(X'X)^-1 X'y
    #cd   <- mldivide( as.matrix(xpx) , as.matrix(xpWy)) #(X'X)^(-1) * X'Wy
    c0   <-  xpxI %*% xpy   #(X'X)^-1 X'y
    cd   <-  xpxI %*% xpWy  #(X'X)^(-1) * X'Wy
    e0   <-  y - x %*% c0
    ed   <- Wy - x %*% cd
    epe0 <- as.real(t(e0) %*% e0)
    eped <- as.real(t(ed) %*% ed)
    epe0d<- as.real(t(ed) %*% e0)    
    rho  <- draw_rho(detval,epe0,eped,epe0d,n,k,rho,a1,a2)

    # update z-values,  
    hh     <- In - rho*W
    hhI    <- qr.solve( hh )
    #mu     <- mldivide( as.matrix(hh), as.matrix(xb) )
    mu     <- hhI %*% xb
    # tauinv <- h'*h
    tauinv <- In - rho * Wadd + rho^2*WpW           #nxn
    aa     <- as.matrix(diag(tauinv))               #nx1
    h      <- ones(n,1)/sqrt(aa)
    c      <- matdiv(-tauinv,aa)
    ctilde <- c - diag(diag(c))
    
    for( initer in 1:nsample){
      for(i in 1:n){
        aa    <- ctilde[i,] %*% z
        muuse <- ( -mu[i, 1] - aa )/h[i,1]
        if( yin[i,1] == 0 ){
          #rigth truncated normal
          t1 <- rtnorm(    1, upper=muuse, mean=0, sd=1)
          #t1 <- rtruncnorm(1, b=muuse, mean=0,sd=1 )
        }else if( yin[i,1] == 1){
          #left truncated normal
          t1 <- rtnorm(    1, lower=muuse, mean=0, sd=1)
          #t1 <- rtruncnorm(1, a=muuse, mean=0,sd=1 )
        }
        z[i,1] <- aa + h[i,1]*t1
      }
    }

    #recalc variables
    y   <- mu + z 
    Wy  <- W %*% y
    xpy <- xp %*% y
    xpWy<- xp  %*% Wy
    
    #After the burn-in phase 
    #save all draws for posterior analysis
    if( iter > nomit ){
      bsave[iter-nomit, ] <- bhat[,1] 
      psave[iter-nomit,1] <- rho
      ymean               <- ymean + y
      
      # compute effects estimates
      rhovec <- as.matrix( rho^(0:(o-1)) ) #100x1
      if( cflag == 1 ){#has intercept
        beff <- as.matrix( bhat[2:length(bhat), 1] )
      }else if(cflag == 0){
        beff <- bhat
      }
      #s    <- mldivide( as.matrix(hh), as.matrix(In) )
      s         <- hhI %*% In
      pdfz      <- dnorm( as.matrix( mu ) ) #standard normal pdf
      diag( dd )<- pdfz
      for( kk in 1:p ){
        avg_direct[kk,1]   <- as.real(t(pdfz) %*% (estimated_diags %*% rhovec * beff[ kk,1 ] /n))
        tmp                <- apply( dd %*% s * beff[kk,1], 2, sum )
        avg_total[kk,1]    <- mean( tmp ) 
        total_obs[ ,kk]    <- total_obs[ ,kk] + tmp 
        avg_indirect[kk,1] <- avg_total[kk,1] - avg_direct[kk,1]
      }  
      total[iter-nomit, ]   <- avg_total[,1]    # an ndraw-nomit x p matrix
      direct[iter-nomit, ]  <- avg_direct[,1]   # an ndraw-nomit x p matrix
      indirect[iter-nomit, ]<- avg_indirect[,1] # an ndraw-nomit x p matrix 
    }#end if iter > nomit
    iter = iter + 1
    setTxtProgressBar(pb, iter)
  } #end sampling loop
  close(pb) #close progress bar
  
  time3         <- Sys.time() - timet
  results$time3 <- time3
  total_obs     <- total_obs/(ndraw-nomit)
  beta          <- t( apply(bsave,2,mean ) )
  rho           <- mean( psave )
  ymean         <- ymean/(ndraw-nomit)
  results$sige  <- 1
  sige          <- 1
  # compute log marginal likelihood based on ymean (just an approximation)
  Wy            <- W %*% ymean;
  AI            <- qr.solve( xp %*% x + sige * TI )
  c0            <- AI %*% ( xp %*% ymean + sige * TIc )
  cd            <- AI %*% ( xp %*% Wy    + sige * TIc)
  e0            <- ymean - x %*% c0
  ed            <-    Wy - x %*% cd
  epe0          <- as.real(t(e0) %*% e0)
  eped          <- as.real(t(ed) %*% ed)
  epe0d         <- as.real(t(ed) %*% e0)
  logdetx       <- log(det(t(x) %*% x + sige*TI))
                           
  if( inform_flag == 0 ){
    mlike = sar_marginal(detval,e0,ed,epe0,eped,epe0d,n,k,logdetx,a1,a2)
  }else if( inform_flag == 1 ) {
    print('sarp_g: informative priors not implemented')
    #mlike = sar_marginal2(detval,e0,ed,epe0,eped,epe0d,nobs,k,logdetx,a1,a2,c,TI,x,y,sige,W)
  }
  
  results$meth      <- 'sarp_g'
  results$ymean     <- ymean
  results$total     <- total
  results$direct    <- direct
  results$indirect  <- indirect
  results$total_obs <- total_obs
  results$beta      <- beta
  results$rho       <- rho
  results$bdraw     <- bsave
  results$pdraw     <- psave
  results$bmean     <- c
  results$bstd      <- sqrt(diag(Ts))
  results$ndraw     <- ndraw
  results$nomit     <- nomit
  results$nsteps    <- nsample
  results$time      <- Sys.time() - timet
  results$a1        <- a1
  results$a2        <- a2
  results$tflag     <- 'plevel'
  results$rmax      <- rmax 
  results$rmin      <- rmin
  results$lflag     <- ldetflag
  results$lndet     <- detval
  results$priorb    <- inform_flag
  results$mlike     <- mlike
  results$names     <- c(colnames(x), 'rho')
  return(results)
}


# =========================================================================
# support functions below
# =========================================================================

# PURPOSE:
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
draw_rho <- function(detval,epe0,eped,epe0d,n,k,rho,a1,a2){
  nmk     <- (n-k)/2
  nrho    <- nrow(detval)
  iota    <- ones(nrho,1)
  #LeSage 2009 page 132 justifies the code below
  z       <-  epe0 * iota - 2*detval[,1] * epe0d + (detval[,1] * detval[,1]) * eped #CHECK epe %*% vs *
  z       <- -nmk * log(z)
  #C = gammaln(nmk)*iota -nmk*log(2*pi)*iota - 0.5*logdetx*iota;
  den     <- detval[,2] + z
  bprior  <- beta_prior(detval[,1],a1,a2)
  den     <- den + log(bprior)
  n       <- length(den)
  y       <- as.matrix( detval[ ,1] )
  adj     <- max(den)
  den     <- den - adj
  x       <- exp(den)
  # trapezoid rule used for numerical integral approximation
  isum    <- sum((y[2:n,1] + y[1:(n-1),1]) * (x[2:n,1] - x[1:(n-1),1])/2) 
  z       <- abs(x/isum)
  den     <- cumsum(z)
  rnd     <- runif(1,min=0,max=1)*sum(z)
  ind     <- which( den <= rnd )
  idraw   <- max(ind) #returns -inf if the ind lis is empty
  if (  idraw > 0 && idraw < nrho ){
    results <- detval[idraw,1]
  }else{
    results <- rho
  }
  return( results )
}


                   
 # PURPOSE: returns a vector of the log-marginal over a grid of rho-values
# -------------------------------------------------------------------------
# USAGE: out = sar_marginal(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2)
# where:       detval = an ngrid x 2 matrix with rho-values and lndet values
#                  e0 = y - x*b0;
#                 ed = Wy - x*bd;
#               epe0 = e0'*e0;
#               eped = ed'*ed;
#              epe0d = ed'*e0;
#               nobs = # of observations
#               nvar = # of explanatory variables
#            logdetx = log(det(x'*x))
#                 a1 = parameter for beta prior on rho
#                 a2 = parameter for beta prior on rho
# -------------------------------------------------------------------------
# RETURNS: out = a structure variable
#        out = log marginal, a vector the length of detval
# -------------------------------------------------------------------------
# NOTES: -this does not take any prior on beta, sigma into account, uses diffuse priors
#         see sar_marginal2() 
# -------------------------------------------------------------------------

# written by:
# James P. LeSage, last updated 3/2010
# Dept of Finance & Economics
# Texas State University-San Marcos
# 601 University Drive
# San Marcos, TX 78666
# jlesage@spatial-econometrics.com                   
sar_marginal <- function(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2){
  n   <- nrow(detval)
  nmk <- (nobs-nvar)/2
  # C is a constant of integration that can vary with nvars, so for model
  # comparisions involving different nvars we need to include this
  bprior <-  beta_prior(detval[ ,1],a1,a2)
  C      <-  log(bprior) + lgamma(nmk) - nmk*log(2*pi) - 0.5*logdetx
  iota   <-  ones(n,1)
  #CHECK THIS %*% eped or * eped
  z      <-  epe0 * iota - 2*detval[ ,1] * epe0d + detval[ ,1] * detval[ ,1 ] * eped
  den    <-  C + detval[ ,2] - nmk * log(z)
  return( as.real(den) )
}
                   
# PURPOSE: returns a vector of the log-marginal over a grid of rho-values
#          for the case of an informative prior on beta
# -------------------------------------------------------------------------
# USAGE: out = sar_marginal(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2)
# where:       detval = an ngrid x 2 matrix with rho-values and lndet values
#                  e0 = y - x*b0;
#                 ed = Wy - x*bd;
#               epe0 = e0'*e0;
#               eped = ed'*ed;
#              epe0d = ed'*e0;
#               nobs = # of observations
#               nvar = # of explanatory variables
#            logdetx = log(det(x'*x))
#                 a1 = parameter for beta prior on rho
#                 a2 = parameter for beta prior on rho
#                 c = prior mean for beta
#                TI = prior var-cov for beta
#                xs = x*sqrt(V) or x if homoscedastic model
#                ys = y*sqrt(V) or y is homoscedastic model
# -------------------------------------------------------------------------
# RETURNS: out = a structure variable
#        out = log marginal, a vector the length of detval
# -------------------------------------------------------------------------
# NOTES: -this is only an approximation based on the posterior mean Vi-estimates
# -------------------------------------------------------------------------

# written by:
# James P. LeSage, last updated 3/2010
# Dept of Finance & Economics
# Texas State University-San Marcos
# 601 University Drive
# San Marcos, TX 78666
# jlesage@spatial-econometrics.com

# sar_marginal2 <- function(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,a1,a2,c,TI,xs,ys,sige,W){
#   n = length(detval);
#   nmk = (nobs-nvar)/2;
#   # C is a constant of integration that can vary with nvars, so for model
#   # comparisions involving different nvars we need to include this
#   bprior = beta_prior(detval(:,1),a1,a2);
#   C = log(bprior) + gammaln(nmk) - nmk*log(2*pi);
#   iota = ones(n,1);
#   z = epe0*iota - 2*detval(:,1)*epe0d + detval(:,1).*detval(:,1)*eped;
#   # add quadratic terms based on prior for beta
#   Q1 = zeros(n,1);
#   Q2 = zeros(n,1);
#   xpxi = inv(xs'*xs);
#   sTI = sige*TI;
#   xpxis = inv(xs'*xs + sTI);
# logdetx = log(det(xpxis));
# C = C - 0.5*logdetx;
#           for i=1:n;
#            rho = detval(i,1);
#            D = speye(nobs) - rho*W;
#            bhat = xpxi*(xs'*D*ys);
#            beta = xpxis*(xs'*D*ys + sTI*c); 
#            Q1(i,1) = (c - beta)'*sTI*(c - beta);
#            Q2(i,1) = (bhat - beta)'*(xs'*xs)*(bhat - beta);
#           end;
# 
# den = C + detval(:,2) - nmk*log(z + Q1 + Q2);
# out = real(den);
