\name{sartobit}
\alias{sartobit}
\alias{sar_tobit_mcmc}
\title{Bayesian estimation of the SAR Tobit model}

\description{
Bayesian estimation of the spatial autoregressive Tobit model (SAR Tobit model).
}

\usage{
sartobit(formula, W, data, ...)

sar_tobit_mcmc(y, X, W, ndraw = 1000, burn.in = 100, thinning = 1, 
  prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0), 
  start = list(rho = 0.75, beta = rep(0, ncol(X)), sige = 1),
  m=10, computeMarginalEffects=FALSE, showProgress=FALSE)
}
\arguments{
  \item{y}{dependent variables. vector of zeros and ones}
  \item{X}{design matrix}
  \item{W}{spatial weight matrix}
  \item{ndraw}{number of MCMC iterations}
  \item{burn.in}{number of MCMC burn-in to be discarded}
  \item{thinning}{MCMC thinning factor, defaults to 1.}
  \item{prior}{A list of prior settings for 
               \eqn{\rho \sim Beta(a1,a2)}{or rho ~ Beta(a1,a2)} 
               and \eqn{\beta \sim N(c,T)}{beta ~ N(c,T)}. Defaults to diffuse prior for beta.}
  \item{start}{list of start values}
  \item{m}{Number of burn-in samples in innermost Gibbs sampler. Defaults to 10.}
  \item{computeMarginalEffects}{Flag if marginal effects are calculated. Defaults to FALSE. 
   We recommend to enable it only when sample size is small.}
  \item{showProgress}{Flag if progress bar should be shown. Defaults to FALSE.}
  \item{formula}{an object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{data}{an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. 
              If not found in data, the variables are taken from environment(formula), typically the environment from which sarprobit is called.}
  \item{...}{additional arguments to be passed}
}
\details{
Bayesian estimates of the spatial autoregressive Tobit model (SAR Tobit model)
\deqn{
  z = \rho W y + X \beta + \epsilon, \epsilon \sim N(0, \sigma^2_{e} I_n)
}{
  z = rho*W*y + X*beta + e, e ~ N(0,sige*I_n)
}
\deqn{
  z = (I_n - \rho W)^{-1} X \beta + (I_n - \rho W)^{-1} \epsilon
}{
  z = (I_n - rho W)^{-1} X beta + (I_n - rho W)^{-1} e
}
\deqn{
  y = max(z, 0)
}
where \eqn{y} \eqn{(n \times 1)}{(n x 1)} is only observed for \eqn{z \ge 0}{z >= 0} and censored to 0 otherwise. 
\eqn{\beta}{beta} is a \eqn{(k \times 1)}{(k x 1)} 
vector of parameters associated with the \eqn{(n \times k)}{(n x k)} data matrix X.

The prior distributions are 
 \eqn{\beta \sim N(c,T)}{beta ~ N(c,T)} 
 and  
 \eqn{\rho \sim Uni(rmin,rmax)}{rho ~ Uni(rmin,rmax)} 
 or 
 \eqn{\rho \sim Beta(a1,a2)}{or rho ~ Beta(a1,a2)}.
}
\value{
Returns a structure of class \code{sartobit}:
%\item{meth }{ 'sartobit' }
\item{beta }{ posterior mean of bhat based on draws }
\item{rho }{ posterior mean of rho based on draws }
\item{bdraw }{ beta draws (ndraw-nomit x nvar) }
\item{pdraw }{ rho  draws (ndraw-nomit x 1) }
\item{sdraw }{ sige draws (ndraw-nomit x 1) }
\item{total }{ a matrix (ndraw,nvars-1) total x-impacts }
\item{direct }{ a matrix (ndraw,nvars-1) direct x-impacts }
\item{indirect }{ a matrix (ndraw,nvars-1) indirect x-impacts }
%\item{total_obs }{ a matrix (ndraw,nvars-1) observation-level total x-impacts }
\item{rdraw }{ r draws (ndraw-nomit x 1) (if m,k input) }
%\item{bmean }{ b prior means, prior.beta from input }
%\item{bstd }{ b prior std deviations sqrt(diag(prior.bcov)) }
\item{nobs }{ # of observations }
\item{nvar }{ # of variables in x-matrix }
\item{ndraw }{ # of draws }
\item{nomit }{ # of initial draws omitted }
\item{nsteps }{ # of samples used by Gibbs sampler for TMVN }
\item{y }{ y-vector from input (nobs x 1) }
\item{zip }{ # of zero y-values }
%\item{yhat }{ mean of posterior predicted (nobs x 1) }
%\item{resid }{ residuals, based on posterior means }
%\item{rsqr }{ r-squared based on posterior means }
%\item{rbar }{ adjusted r-squared }
\item{a1 }{ a1 parameter for beta prior on rho from input, or default value }
\item{a2 }{ a2 parameter for beta prior on rho from input, or default value }
\item{time }{ total time taken }  
\item{rmax }{ 1/max eigenvalue of W (or rmax if input) }
\item{rmin }{ 1/min eigenvalue of W (or rmin if input) }          
\item{tflag }{ 'plevel' (default) for printing p-levels; 
               'tstat' for printing bogus t-statistics }
\item{lflag }{ lflag from input }
\item{cflag }{ 1 for intercept term, 0 for no intercept term }
%\item{iter }{  prior.iter option from input }
\item{lndet }{ a matrix containing log-determinant information
               (for use in later function calls to save time) }
%\item{mlike }{ log marginal likelihood (a vector ranging over
%               rho values that can be integrated for model comparison) }
}
\references{
LeSage, J. and Pace, R. K. (2009), \emph{Introduction to Spatial Econometrics}, CRC Press, chapter 10,
 section 10.3, 299--304 
}

\author{
adapted to and optimized for R by Stefan Wilhelm <wilhelm@financial.com> 
based on Matlab code from James P. LeSage}

\seealso{
\code{\link{sarprobit}}, \code{\link{sarorderedprobit}} or \code{\link{semprobit}} 
for SAR probit/SAR Ordered Probit/ SEM probit model fitting
}

\examples{
\dontrun{
# Example from LeSage/Pace (2009), section 10.3.1, p. 302-304
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

# full information
yfull <- y

# set negative values to zero to reflect sample truncation
ind <- which(y <=0)
y[ind] <- 0

# Fit SAR (with complete information)
fit_sar <- sartobit(yfull ~ X-1, W,ndraw=1000,burn.in=200, showProgress=FALSE)
summary(fit_sar)

# Fit SAR Tobit (with approx. 50\% censored observations)
fit_sartobit <- sartobit(y ~ x,W,ndraw=1000,burn.in=200, showProgress=TRUE)

par(mfrow=c(2,2))
for (i in 1:4) {
 ylim1 <- range(fit_sar$B[,i], fit_sartobit$B[,i])
 plot(fit_sar$B[,i], type="l", ylim=ylim1, main=fit_sartobit$names[i], col="red")
 lines(fit_sartobit$B[,i], col="green")
 legend("topleft", legend=c("SAR", "SAR Tobit"), col=c("red", "green"), 
   lty=1, bty="n")
}

# Fit SAR Tobit (with approx. 50\% censored observations)
fit_sartobit <- sartobit(y ~ x,W,ndraw=1000,burn.in=0, showProgress=TRUE, 
    computeMarginalEffects=TRUE)
# Print SAR Tobit marginal effects
impacts(fit_sartobit)
#--------Marginal Effects--------
#
#(a) Direct effects
#  lower_005 posterior_mean upper_095
#x     1.013          1.092     1.176
#
#(b) Indirect effects
#  lower_005 posterior_mean upper_095
#x     2.583          2.800     3.011
#
#(c) Total effects
#  lower_005 posterior_mean upper_095
#x     3.597          3.892     4.183
#mfx <- marginal.effects(fit_sartobit)
}
}
