\name{marginal.effects}
\alias{marginal.effects}
\alias{marginal.effects.sarprobit}
\alias{marginal.effects.sartobit}
\alias{impacts.sarprobit}
\alias{impacts.sartobit}

\title{Marginal effects for spatial probit and Tobit models (SAR probit, SAR Tobit)}
\description{
Estimate marginal effects (average direct, indirect and total impacts)
for the SAR probit and SAR Tobit model.
}
\usage{
\method{marginal.effects}{sarprobit}(object, o = 100, ...)
\method{marginal.effects}{sartobit}(object, o = 100, ...)
\method{impacts}{sarprobit}(obj, file=NULL, 
  digits = max(3, getOption("digits")-3), ...)
\method{impacts}{sartobit}(obj, file=NULL, 
  digits = max(3, getOption("digits")-3), ...)
}

\arguments{
  \item{object}{Estimated model of class \code{sarprobit} or \code{sartobit}}
  \item{obj}{Estimated model of class \code{sarprobit} or \code{sartobit}}
  \item{o}{maximum value for the power \eqn{tr(W^i), i=1,...,o} to be estimated}
  \item{digits}{number of digits for printing}
  \item{file}{Output to file or console}
  \item{...}{additional parameters}
}
\details{
\code{impacts()} will extract and print the marginal effects from a fitted model,
while \code{marginal.effects(x)} will estimate the marginal effects anew for a fitted
model.

In spatial models, a change in some explanatory variable \eqn{x_{ir}} for observation \eqn{i} 
will not only affect the observations \eqn{y_i} directly (direct impact), 
but also affect neighboring observations \eqn{y_j} (indirect impact). 
These impacts potentially also include feedback loops from observation 
\eqn{i} to observation \eqn{j} and back to \eqn{i}. 
(see LeSage (2009), section 2.7 for interpreting parameter estimates in spatial models).

For the \eqn{r}-th non-constant explanatory variable, 
let \eqn{S_r(W)} be the \eqn{n \times n}{n x n} matrix 
that captures the impacts from observation \eqn{i} to \eqn{j}.

The direct impact of a change in \eqn{x_{ir}} on its own observation \eqn{y_i} can be written as
\deqn{
  \frac{\partial y_i}{\partial x_{ir}} = S_r(W)_{ii}
}{
  d y_i / d x_{ir} = S_r(W)_{ii}
} 
and the indirect impact from observation \eqn{j} to observation \eqn{i} as
\deqn{
  \frac{\partial y_i}{\partial x_{jr}} = S_r(W)_{ij}
}{
  d y_i / d x_{jr} = S_r(W)_{ij}
}.

LeSage(2009) proposed summary measures for direct, indirect and total effects, 
e.g. averaged direct impacts across all \eqn{n} observations.
See LeSage(2009), section 5.6.2., p.149/150 for marginal effects estimation in general spatial models
and section 10.1.6, p.293 for marginal effects in SAR probit models.

We implement these three summary measures:
\enumerate{

\item average direct impacts:   
\deqn{M_r(D) = \bar{S_r(W)_{ii}} = n^{-1} tr(S_r(W))}{M_r(D) = average(S_r(W)_{ii}) = n^{-1} tr(S_r(W))}

\item average total impacts:    
\deqn{M_r(T) = n^{-1} 1'_n S_r(W) 1_n}

\item average indirect impacts: 
\deqn{M_r(I) = M_r(T) - M_r(D)}
}

The average direct impact is the average of the diagonal elements, the average total
impacts is the mean of the row (column) sums.

For the average direct impacts \eqn{M_r(D)}, 
there are efficient approaches available, see LeSage (2009), chapter 4, pp.114/115.

The computation of the average total effects \eqn{M_r(T)} and hence 
also the average indirect effects \eqn{M_r(I)} are more subtle, 
as \eqn{S_r(W)} is a dense \eqn{n \times n}{n x n} matrix.
In the LeSage Spatial Econometrics Toolbox for MATLAB (March 2010), 
the implementation in \code{sarp_g} computes the matrix inverse of \eqn{S= (I_n -  \rho W)}
which all the negative consequences for large n.
We implemented \eqn{n^{-1} 1'_n S_r(W) 1_n} via a QR decomposition of 
\eqn{S = (I_n -  \rho W)} (already available from a previous step) and solving a linear equation, 
which is less costly and will work better for large \eqn{n}.

\bold{SAR probit model}\cr

Specifically, for the SAR probit model the \eqn{n \times n}{n x n} matrix of marginal effects is
\deqn{
  S_r(W) = \frac{\partial E[y | x_r]}{\partial x_{r}'} = \phi\left((I_n -  \rho W)^{-1} \bar{x}_r \beta_r \right) \odot (I_n -  \rho W)^{-1} I_n \beta_r 
}{
  S_r(W) = d E[y | x_r] / d x_r' = phi((I_n - rho W)^{-1} I_n mean(x_r) beta_r) * (I_n -  rho W)^{-1} I_n beta_r
}

\bold{SAR Tobit model}\cr

Specifically, for the SAR Tobit model the \eqn{n \times n}{n x n} matrix of marginal effects is
\deqn{
  S_r(W) = \frac{\partial E[y | x_r]}{\partial x_{r}'} = \Phi\left((I_n -  \rho W)^{-1} \bar{x}_r \beta_r / \sigma \right) \odot (I_n -  \rho W)^{-1} I_n \beta_r 
}{
  S_r(W) = d E[y | x_r] / d x_r' = Phi((I_n - rho W)^{-1} I_n mean(x_r) beta_r / sigma) * (I_n -  rho W)^{-1} I_n beta_r
}

}

\section{Warning}{
1. Although the direct impacts can be efficiently estimated, the computation
of the indirect effects require the inversion of a \eqn{n \times n}{n x n} matrix
and will break down for large \eqn{n}.

2. \eqn{tr(W^i)} is determined with simulation, 
so different calls to this method will produce different estimates.
}

\references{
LeSage, J. and Pace, R. K. (2009), \emph{Introduction to Spatial Econometrics}, CRC Press 
}
\author{
Stefan Wilhelm <wilhelm@financial.com>
}

\seealso{
 \code{\link{marginal.effects.sartobit}}
}

\examples{
\dontrun{
require(spatialprobit)

# number of observations
n <- 10

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
set.seed(12345)
sarprobit.fit1 <- sar_probit_mcmc(y, X, W, ndraw=500, burn.in=100, 
  thinning=1, prior=NULL, computeMarginalEffects=TRUE)
summary(sarprobit.fit1)

# print impacts
impacts(sarprobit.fit1)

################################################################################
#
# Example from LeSage/Pace (2009), section 10.3.1, p. 302-304
#
################################################################################

# Value of "a" is not stated in book! 
# Assuming a=-1 which gives approx. 50% censoring
 library(spatialprobit)

a <- -1   # control degree of censored observation
n <- 1000
rho <- 0.7
beta <- c(0, 2)
sige <- 0.5
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
x <- runif(n, a, 1)
X <- cbind(1, x)
eps <- rnorm(n, sd=sqrt(sige))
param <- c(beta, sige, rho)

# random locational coordinates and 6 nearest neighbors
lat <- rnorm(n)
long <- rnorm(n)
W <- kNearestNeighbors(lat, long, k=6)

y <- as.double(solve(I_n - rho * W) \%*\% (X \%*\% beta + eps))
table(y > 0)

# set negative values to zero to reflect sample truncation
ind <- which(y <=0)
y[ind] <- 0

# Fit SAR Tobit (with approx. 50\% censored observations)
fit_sartobit <- sartobit(y ~ x,W,ndraw=1000,burn.in=200, showProgress=TRUE)

# print impacts
impacts(fit_sartobit)
}
}
