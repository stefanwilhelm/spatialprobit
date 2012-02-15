\name{sar_eigs}
\alias{sar_eigs}
\title{
compute the eigenvalues for the spatial weight matrix W
}
\description{
compute the eigenvalues \eqn{\lambda}{lambda} for the spatial weight matrix \eqn{W} and lower and upper bound
for parameter \eqn{\rho}{rho}.
}
\usage{
sar_eigs(eflag, W)
}
\arguments{
  \item{eflag}{if eflag==1, then eigen values
}
  \item{W}{spatial weight matrix}
}
\value{
function returns a list of
\item{rmin }{minimum value for \eqn{\rho}{rho}. if eflag==1, then \eqn{1/\lambda_{min}}{1/min(lambda)} , else -1}
\item{rmax }{maximum value for \eqn{\rho}{rho}. Always 1.}
\item{time}{execution time}
}

\author{
James P. LeSage, Adapted to R by Miguel Godinho de Matos <miguelgodinhomatos@cmu.edu>
}

\examples{
set.seed(123)
# sparse matrix representation for spatial weight matrix W (d x d) 
# and m nearest neighbors
d <- 100
m <- 6
W <- sparseMatrix(i=rep(1:d, each=m), 
  j=replicate(d, sample(x=1:d, size=m, replace=FALSE)), x=1/m, dims=c(d, d))
sar_eigs(eflag=1, W)    
}
