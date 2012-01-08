\name{lndetfull}
\alias{lndetfull}
\title{
This function computes the Chebyshev log determinant approximation. 
}
\description{
This function computes the Chebyshev log determinant approximation.
}
\usage{
lndetfull(W, lmin, lmax)
}
\arguments{
  \item{W}{spatial weight matrix}
  \item{lmin}{
%%     ~~Describe \code{lmin} here~~
}
  \item{lmax}{
%%     ~~Describe \code{lmax} here~~
}
}

\references{
Pace and LeSage 2004, Chebyshev Approximation of log-determinants of spatial weight matrices, Computational Statistics and Data Analysis, vol 45
}
\author{
Miguel Godinho de Matos <miguelgodinhomatos@cmu.edu>
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function( W, lmin, lmax ){
  rvec    <- seq(lmin, lmax, 0.01)
  In      <- speye(nrow(W))
  niter   <- length(rvec)
  
  results       <- NULL
  results$rho   <- NULL
  results$lndet <- NULL
  
  for(i in 1:niter){
    rho             <- rvec[i]
    z               <- In - rho*W
    results$rho[i]  <- rho 
    results$lndet[i]<- as.real(determinant(z, log=TRUE)$modulus) 
  }  
  return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
