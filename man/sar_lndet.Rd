\name{sar_lndet}
\alias{sar_lndet}
\title{
compute the log determinant \eqn{\ln{|I_n - \rho W|}}{ln|I_n - rho*W|} using the user-selected (or default) method
}
\description{
compute the log determinant \eqn{\ln{|I_n - \rho W|}}{ln|I_n - rho*W|} using the user-selected (or default) method
}
\usage{
sar_lndet(ldetflag, W, rmin, rmax)
}
\arguments{
  \item{ldetflag}{flag to compute the exact or approximate log determinant (Chebychev approximation, Pace and Barry approximation). See details.}
  \item{W}{spatial weight matrix}
  \item{rmin}{minimum eigen value}
  \item{rmax}{maximum eigen value}
}

\details{
 Chebychev approximation, 
 Pace and Barry approximation
}

\value{
\item{detval}{a 2-column \code{Matrix} with gridpoints from \code{rmin},\ldots,\code{rmax}}
\item{time}{execution time}
}

\references{
John LeSage, R. Kelley Pace: Introduction to Spatial Econometrics. CRC Press. 2009
}

\author{
James P. LeSage, Adapted to R by Miguel Godinho de Matos <miguelgodinhomatos@cmu.edu>
}

\seealso{\code{\link{lndetfull}}}

\examples{
require(Matrix)

# sparse matrix representation for spatial weight matrix W (d x d) and m nearest neighbors
d <- 100
m <- 6
W <- sparseMatrix(i=rep(1:d, each=m), j=replicate(d, sample(x=1:d, size=m, replace=FALSE)), x=1/m, dims=c(d, d))

# exact log determinant
#det1 <- sar_lndet(ldetflag=0, W, rmin=-1, rmax=1)

# Chebychev approximation of log determinant
#det2 <- sar_lndet(ldetflag=1, W, rmin=-1, rmax=1)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
