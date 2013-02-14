\name{plot.sarprobit}
\alias{plot.sarprobit}
\alias{plot.semprobit}
\alias{plot.sartobit}
\title{
Plot Diagnostics for \code{sarprobit}, \code{semprobit} or \code{sartobit} objects
}
\description{
Three plots (selectable by which) are currently available:
MCMC trace plots, autocorrelation plots and posterior density plots.
}
\usage{
\method{plot}{sarprobit}(x, 
which = c(1, 2, 3), 
ask = prod(par("mfcol")) < length(which) && dev.interactive(),
...,
trueparam = NULL)

\method{plot}{semprobit}(x, 
which = c(1, 2, 3), 
ask = prod(par("mfcol")) < length(which) && dev.interactive(),
...,
trueparam = NULL)

\method{plot}{sartobit}(x, 
which = c(1, 2, 3), 
ask = prod(par("mfcol")) < length(which) && dev.interactive(),
...,
trueparam = NULL)
}

\arguments{
  \item{x}{a \code{sarprobit} or \code{semprobit} object}
  \item{which}{if a subset of the plots is required, specify a subset of the numbers 1:3.}
  \item{ask}{logical; if TRUE, the user is asked before each plot, see \code{par(ask=.)}.}
  \item{...}{other parameters to be passed through to plotting functions.}
  \item{trueparam}{a vector of "true" parameter values to be marked as vertical lines in posterior density plot}
}

\author{
Stefan Wilhelm <wilhelm@financial.com>
}

\examples{
\dontrun{
library(Matrix)
# number of observations
n <- 200

# true parameters
beta <- c(0, 1, -1)
rho <- 0.75

# design matrix with two standard normal variates as "covariates"
X <- cbind(intercept=1, x=rnorm(n), y=rnorm(n))

# sparse identity matrix
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)

# number of nearest neighbors in spatial weight matrix W
m <- 6

# spatial weight matrix with m=6 nearest neighbors
W <- sparseMatrix(i=rep(1:n, each=m), 
  j=replicate(n, sample(x=1:n, size=m, replace=FALSE)), x=1/m, dims=c(n, n))

# innovations
eps <- rnorm(n=n, mean=0, sd=1)

# generate data from model 
S <- I_n - rho * W
z <- solve(qr(S), X \%*\% beta + eps)
y <- as.vector(z >= 0)  # 0 or 1, FALSE or TRUE

# estimate SAR probit model
fit1 <- sar_probit_mcmc(y, X, W, ndraw=500, burn.in=100, thinning=1, prior=NULL)
plot(fit1, which=c(1,3))
}
}