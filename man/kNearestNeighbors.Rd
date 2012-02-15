\name{kNearestNeighbors}
\alias{kNearestNeighbors}
\title{Build Spatial Weight Matrix from k nearest neighbors}
\description{
Build a spatial weight matrix using the k nearest neighbors of (x, y) coordinates 
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
Determine the k nearest neighbors for a set of points represented by (x, y) 
coordinates and build a spatial weight matrix.
}

\value{
The method returns a sparse spatial weight matrix W with dimension 
\code{nrow(X)} \eqn{\times}{x} \code{nrow(X)} and 
k non-zero entries per row which represent the k nearest neighbors.
}

\author{
Stefan Wilhelm <stefan.wilhelm@financial.com>
}
\examples{
# build spatial weight matrix W from (x,y) coordinates
W <- kNearestNeighbors(x=rnorm(100), y=rnorm(100), k=6)
}