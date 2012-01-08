\name{buildSpatialWeightMatrix}
\alias{buildSpatialWeightMatrix}
\title{Build Spatial Weight Matrix}
\description{
Build a spatial weight matrix using the m nearest neighbors of (x, y) coordinates 
}
\usage{
buildSpatialWeightMatrix(X, m = 6)
}

\arguments{
  \item{X}{number of locational coordinates}
  \item{m}{number of nearest neighbors}
}

\details{
Determine the m nearest neighbors for a set of points represented by (x, y) coordinates
and build a spatial weight matrix.
}

\value{
The method returns a sparse spatial weight matrix W with dimension 
\code{nrow(X)} \eqn{\times}{x} \code{nrow(X)} and 
m non-zero entries per row which represent the m nearest neighbors.
}

\author{
Stefan Wilhelm <stefan.wilhelm@financial.com>
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(X, m=6) {
  n <- nrow(X)
  # spatial weight matrix W based on the 6 nearest neighbors
  D <- matrix(NA, n, m)  # (n x m) index matrix to the 6 nearest neigbhors from point i
  for (i in 1:n) {
    p <- X[i,]
    # euclidean dist from all points to p
    d <- sqrt((X[,1] - p[1])^2 + (X[,2] - p[2])^2)
    # determine the m nearest neighbors (rank 1 is the point itself, ranks 2..(m+1) are the m nearest neighbors
    # TODO: if 2 points have the same distance, e.g. ranks=1, 2, 2.5, 2.5 are odd
    D[i, ] <- which(rank(d) \%in\% 2:(m+1))
  }
  # sparse matrix representation for spatial weight matrix W
  W <- sparseMatrix(i = rep(1:n, m), j=as.vector(D), x=1/m)
  #W <- sparseMatrix(i = 1:n, j=1:n, x=1) # identity matrix for testing
  return(W)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
