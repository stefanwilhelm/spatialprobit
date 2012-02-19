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
W will be a sparse matrix representation.
}

\value{
The method returns a sparse spatial weight matrix W with dimension 
(n \eqn{\times}{x} n) and \code{k} non-zero entries per row 
which represent the \code{k} nearest neighbors.
}

\author{
Stefan Wilhelm <wilhelm@financial.com>
}
\examples{
require(Matrix)
# build spatial weight matrix W from (x,y) coordinates
W <- kNearestNeighbors(x=rnorm(100), y=rnorm(100), k=6)
image(W, main="spatial weight matrix W")
}