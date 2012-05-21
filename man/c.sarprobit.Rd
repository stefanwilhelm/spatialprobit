\name{c.sarprobit}
\alias{c.sarprobit}
\title{Combine different SAR probit estimates into one}
\description{
This method combines SAR probit estimates into one SAR probit object,
e.g. when collecting the results of a parallel MCMC.
}
\usage{
\method{c}{sarprobit}(...)
}
\arguments{
  \item{...}{A vector of \code{sarprobit} objects.}
}
\value{
This functions returns an object of class \code{sarprobit}.
}
\author{
Stefan Wilhelm <wilhelm@financial.com>
}

\seealso{
\code{\link{sarprobit}} for SAR probit model fitting
}

\examples{
\dontrun{
## parallel estimation using mclapply() under Unix (will not work on Windows)
library(parallel)
mc <- 2 # number of cores; set as appropriate to your hardware
run1 <- function(...) sar_probit_mcmc(y, X, W, ndraw=500, burn.in=200, thinning=1)
system.time( {
## To make this reproducible:
set.seed(123, "L'Ecuyer")
sarprobit.res <- do.call(c, mclapply(seq_len(mc), run1))
})
summary(sarprobit.res)


## parallel estimation using parLapply() under Windows
library(parallel)
ndraw <- 1000  # the number of MCMC draws
mc <- 4        # the number of cores; set as appropriate to your hardware

run1 <- function(...) {
  args <- list(...)
  library(spatialprobit)
  sar_probit_mcmc(y=args$y, X=args$X2, W=args$W, ndraw=args$ndraw, burn.in=100, thinning=1)
}

parallelEstimation <- function(mc, ndraw, y, X, W) {
  cl <- makeCluster(mc)
  ## To make this reproducible:
  clusterSetRNGStream(cl, 123)
  library(spatialprobit) # needed for c() method on master
  sarprobit.res <- do.call(c, parLapply(cl, seq_len(mc), run1, y=y, X2=X, W=W, ndraw=ndraw/mc))
  stopCluster(cl)
  return(sarprobit.res)
}

# parallel estimation using 1, 2, 4 and 8 cores
system.time(p1 <- parallelEstimation(mc=1, ndraw=5000, y=y, X=X, W=W))
system.time(p2 <- parallelEstimation(mc=2, ndraw=5000, y=y, X=X, W=W))
system.time(p4 <- parallelEstimation(mc=4, ndraw=5000, y=y, X=X, W=W))
system.time(p8 <- parallelEstimation(mc=8, ndraw=5000, y=y, X=X, W=W))
}
}