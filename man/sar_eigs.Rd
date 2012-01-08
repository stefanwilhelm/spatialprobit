\name{sar_eigs}
\alias{sar_eigs}
\title{
compute the eigenvalues for the spatial weight matrix W
}
\description{
compute the eigenvalues for the spatial weight matrix W
}
\usage{
sar_eigs(eflag, W)
}
\arguments{
  \item{eflag}{
%%     ~~Describe \code{eflag} here~~
}
  \item{W}{spatial weight matrix}
}
\value{
function returns a list of
\item{rmin }{minimum eigen value}
\item{rmax }{maximum eigen value}
\item{time}{execution time}
}

\author{
James P. LeSage, Adapted to R by Miguel Godinho de Matos <miguelgodinhomatos@cmu.edu>
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (eflag, W) 
{
    results <- NULL
    results$rmin <- -1
    results$rmax <- 1
    results$time <- 0
    if (eflag == 1) {
        t0 <- Sys.time()
        lambda <- as.real(eigen(W, only.values = TRUE)$values)
        results$rmin <- 1/min(lambda)
        results$rmax <- 1
        results$time <- Sys.time() - t0
    }
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
