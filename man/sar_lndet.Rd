\name{sar_lndet}
\alias{sar_lndet}
\alias{lndetfull}
\alias{lndetChebyshev}

\title{
Approximation of the log determinant \eqn{\ln{|I_n - \rho W|}}{ln|I_n - rho*W|} of a spatial weight matrix 
}

\description{
Compute the log determinant \eqn{\ln{|I_n - \rho W|}}{ln|I_n - rho*W|} of a spatial weight matrix W
using either the exact approach, or using some approximations like 
the Chebyshev log determinant approximation or Pace and Barry approximation.
}
\usage{
sar_lndet(ldetflag, W, rmin, rmax)
lndetfull(W, rmin, rmax)
lndetChebyshev(W, rmin, rmax)
}
\arguments{
  \item{ldetflag}{flag to compute the exact or approximate log-determinant (Chebychev approximation, Pace and Barry approximation). See details.}
  \item{W}{spatial weight matrix}
  \item{rmin}{minimum eigen value}
  \item{rmax}{maximum eigen value}
}

\details{
\code{ldetflag=0} will compute the exact log-determinant at some gridpoints, whereas
\code{ldetflag=1} will compute the Chebyshev log-determinant approximation.

\bold{Exact log-determinant:}\cr
The exact log determinant \eqn{\ln|I_n - \rho W|}{ln|I_n - rho W|} 
is evaluated on a grid from \eqn{\rho=-1,...,+1}{\rho=-1,\ldots,+1}. The gridpoints
are then approximated by a spline function.

\bold{Chebychev approximation:}\cr
This option provides the Chebyshev log-determinant approximation as proposed 
by Pace and LeSage (2004). The implementation is faster than the full 
log-determinant method.

\bold{Pace and Barry approximation:}\cr
To be implemented.
}

\value{
\item{detval}{a 2-column \code{Matrix} with gridpoints for rho from \code{rmin},\ldots,\code{rmax} and corresponding log-determinant}
\item{time}{execution time}
}

\references{
Pace, R. K. and Barry, R. (1997), Quick Computation of Spatial Autoregressive Estimators, \emph{Geographical Analysis}, \bold{29}, 232--247

LeSage, J. and Pace, R. K. (2004), Chebyshev Approximation of log-determinants of spatial weight matrices, \emph{Computational Statistics and Data Analysis}, \bold{45}, 179--196.

LeSage, J. and Pace, R. K. (2009), \emph{Introduction to Spatial Econometrics}, CRC Press, chapter 4 
}

\author{
James P. LeSage, Adapted to R by Miguel Godinho de Matos <miguelgodinhomatos@cmu.edu>
}

\examples{
require(Matrix)

# sparse matrix representation for spatial weight matrix W (d x d) 
# and m nearest neighbors
d <- 10
m <- 3
W <- sparseMatrix(i=rep(1:d, each=m), 
  j=replicate(d, sample(x=1:d, size=m, replace=FALSE)), x=1/m, dims=c(d, d))

# exact log determinant
ldet1 <- sar_lndet(ldetflag=0, W, rmin=-1, rmax=1)

# Chebychev approximation of log determinant
ldet2 <- sar_lndet(ldetflag=1, W, rmin=-1, rmax=1)

plot(ldet1$detval[,1], ldet1$detval[,2], type="l", col="black", 
  xlab="rho", ylab="ln|I_n - rho W|",
  main="Log-determinant ln|I_n - rho W| Interpolations")
lines(ldet2$detval[,1], ldet2$detval[,2], type="l", col="red")
legend("bottomleft", legend=c("Exact log-determinant", "Chebychev approximation"), 
  lty=1, lwd=1, col=c("black","red"), bty="n")

}
