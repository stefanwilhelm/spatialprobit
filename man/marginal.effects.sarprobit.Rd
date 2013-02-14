\name{marginal.effects}
\alias{marginal.effects}
\alias{marginal.effects.sarprobit}
\alias{impacts.sarprobit}
\title{Compute marginal effects for the SAR probit model}
\description{
Estimate marginal effects (average direct, indirect and total impacts)
for the SAR probit model.
}
\usage{
\method{marginal.effects}{sarprobit}(object, o = 100, ...)
\method{impacts}{sarprobit}(obj, file=NULL, 
  digits = max(3, getOption("digits")-3), ...)
}

\arguments{
  \item{object}{Estimated SAR probit model of class \code{sarprobit}}
  \item{obj}{Estimated SAR probit model of class \code{sarprobit}}
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

For the \eqn{r}-th non-constant explanatory variable, \eqn{S_r(W)} is defined as
\deqn{S_r(W) = \phi\left((I_n -  \rho W)^{-1} \bar{x}_r \beta_r \right) \odot (I_n -  \rho W)^{-1} I_n \beta_r}

The direct impact of a change in \eqn{x_{ir}} on its 
own observation \eqn{y_i} can be written as
\deqn{
  \frac{\partial y_i}{\partial x_{ir}} = S_r(W)_{ii}
}, 
the indirect impact from observation \eqn{j} to observation \eqn{i} as
\deqn{
  \frac{\partial y_i}{\partial x_{jr}} = S_r(W)_{ij}
}.

LeSage(2009) proposed summary measures for direct, indirect and total effects, 
e.g. averaged direct impacts across all \eqn{n} observations.
See LeSage(2009), section 5.6.2., p.149/150 for marginal effects estimation in general spatial models
and section 10.1.6, p.293 for marginal effects in SAR probit models.

We implement these summary measures: \cr
average direct impacts:   \eqn{M_r(D) = average S_r(W)_{ii} = n^{-1} tr(S_r(W))} \cr           
average total impacts:    \eqn{M_r(T) = n^{-1} 1'_n S_r(W) 1_n} \cr     
average indirect impacts: \eqn{M_r(I) = M_r(T) - M_r(D)} \cr

Specifically, for the SAR probit model the \eqn{n \times n}{n x n} matrix of marginal effects is
\deqn{
  \frac{\partial E[y | x_r]}{\partial x_{r}'} = \phi\left((I_n -  \rho W)^{-1} \bar{x}_r \beta_r \right) \odot (I_n -  \rho W)^{-1} I_n \beta_r 
}{
 d E[y | x_r] / d x_r' = phi((In -  rho W)^{-1} I_n mean(x_r) beta_r) * (I_n -  \rho W)^{-1} I_n \beta_r
}
The average direct impact is the average of the diagonal elements, the average total
impacts is the mean of the row (column) sums.

For the average direct impacts \eqn{M_r(D)}, 
there are efficient approaches available, see LeSage (2009), chapter 4, pp.114/115.

The computation of the average total effects \eqn{M_r(T)} and hence 
also the average indirect effects \eqn{M_r(I)} are more subtle, 
as \eqn{S_r(W)} is a dense \eqn{n \times n}{n x n} matrix.
In the LeSage Spatial Econometrics Toolbox for MATLAB (March 2010), 
the implementation in \code{sarp_g} computes the matrix inverse of \eqn{S=I_n -  \rho W}
which all the negative consequences for large n.
We implemented \eqn{n^{-1} 1'_n S_r(W) 1_n} via a QR decomposition of 
\eqn{S=I_n -  \rho W} (already available from a previous step) and solving a linear equation, 
which is less costly and will work better for large \eqn{n}.
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

\examples{
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
sarprobit.fit1 <- sar_probit_mcmc(y, X, W, ndraw=1000, burn.in=50, 
  thinning=1, prior=NULL, computeMarginalEffects=TRUE)
summary(sarprobit.fit1)

# average direct effects
colMeans(sarprobit.fit1$direct)

me <- marginal.effects(sarprobit.fit1)
colMeans(me$direct)

plot(density(me$direct[,1]))
lines(density(sarprobit.fit1$direct[,1]), col="red")
plot(me$direct[,1],sarprobit.fit1$direct[,1]) 
abline(a=0, b=1, lty=1, col="red")

# problem: both should be the same, wrong: 
# tr(W^i) is determined with simulation, so there will be different realisations in sarprobit.fit1
# and marginal.effects.sarprobit()
quantile(me$direct[,1], prob=c(0.025, 0.975)) # 95% confidence interval
quantile(sarprobit.fit1$direct[,1], prob=c(0.025, 0.975))
}
