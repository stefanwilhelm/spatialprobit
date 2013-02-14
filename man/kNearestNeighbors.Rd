\name{kNearestNeighbors}
\alias{kNearestNeighbors}
\title{Build Spatial Weight Matrix from k Nearest Neighbors}
\description{
Build a spatial weight matrix W using the k nearest neighbors of (x, y) coordinates 
}
\usage{
kNearestNeighbors(x, y, k = 6)
}

\arguments{
  \item{x}{x coordinate}
  \item{y}{y coordinate}
  \item{k}{number of nearest neighbors}
}

\details{
Determine the k nearest neighbors for a set of n points represented by (x, y) 
coordinates and build a spatial weight matrix W (n \eqn{\times}{x} n).
W will be a sparse matrix representation and row-standardised.

This method is a convenience method for quickly creating a spatial weights matrix based
on planar coordinates. More ways to create W are 
available in \code{\link[spdep]{knearneigh}} of package \code{spdep}.
}

\value{
The method returns a sparse spatial weight matrix W with dimension 
(n \eqn{\times}{x} n) and \code{k} non-zero entries per row 
which represent the \code{k} nearest neighbors.
}

\seealso{
 \code{\link[spdep]{nb2listw}} and \code{\link[spdep]{knearneigh}}
 for computation of neighbors lists, spatial weights and standardisation.
}

\author{
Stefan Wilhelm <wilhelm@financial.com>
}
\examples{
require(Matrix)
# build spatial weight matrix W from random (x,y) coordinates
W <- kNearestNeighbors(x=rnorm(100), y=rnorm(100), k=6)
image(W, main="spatial weight matrix W")
}