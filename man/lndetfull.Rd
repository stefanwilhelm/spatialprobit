\name{lndetfull}
\alias{lndetfull}
\title{
Computes Pace and Barry's grid for log det(I-rho*W) using sparse matrices.
Chebyshev log determinant approximation for a spatial weight matrix. 
}
\description{
This function computes the Chebyshev log determinant approximation for
\eqn{\ln{|I_n - \rho W|}}{ln|I_n - rho*W|}
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

\value{
\item{lndet}{a vector of log determinants for 0 < rho < 1}
\item{rho}{a vector of \eqn{rho}{rho} values associated with lndet values}
}

\references{
LeSage, J. and Pace, R. K. (2009), \emph{Introduction to Spatial Econometrics}, CRC Press 

LeSage, J. and Pace, R. K. (2004), Chebyshev Approximation of log-determinants of spatial weight matrices, \emph{Computational Statistics and Data Analysis}, \bold{45}

Pace, R. K. and Barry, R. (1997), Quick Computation of Spatial Autoregressive Estimators, \emph{Geographical Analysis}
}

\author{Miguel Godinho de Matos <miguelgodinhomatos@cmu.edu>}

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
