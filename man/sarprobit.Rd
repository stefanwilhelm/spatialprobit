\name{sarprobit}
\alias{sar_probit_mcmc}
\alias{sarprobit}
\title{Bayesian estimation of the SAR probit model}

\description{
Bayesian estimation of the spatial autoregressive probit model (SAR probit model).
}

\usage{
sarprobit(formula, W, data, subset, ...)

sar_probit_mcmc(y, X, W, ndraw = 1000, burn.in = 100, thinning = 1, 
  prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0), 
  start = list(rho = 0.75, beta = rep(0, ncol(X))),
  m=10, computeMarginalEffects=TRUE, showProgress=FALSE)
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
  \item{computeMarginalEffects}{Flag if marginal effects are calculated. Defaults to TRUE}
  \item{showProgress}{Flag if progress bar should be shown. Defaults to FALSE.}
  \item{formula}{an object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{data}{an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. 
              If not found in data, the variables are taken from environment(formula), typically the environment from which sarprobit is called.}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{...}{additional arguments to be passed}
}
\details{
Bayesian estimates of the spatial autoregressive probit model (SAR probit model)
\deqn{
  z = \rho W z + X \beta + \epsilon, \epsilon \sim N(0, I_n)
}{
  z = rho*W*z + X*beta + e, e ~ N(0,I_n)
}
\deqn{
  z = (I_n - \rho W)^{-1} X \beta + (I_n - \rho W)^{-1} \epsilon
}{
  z = (I_n - rho W)^{-1} X beta + (I_n - rho W)^{-1} e
}    
where y is a binary 0,1 \eqn{(n \times 1)}{(n x 1)} vector of observations for z < 0 and z >= 0. 
\eqn{\beta}{beta} is a \eqn{(k \times 1)}{(k x 1)} 
vector of parameters associated with the \eqn{(n \times k)}{(n x k)} data matrix X.
The error variance \eqn{\sigma_e}{sigma_e} is set to 1 for identification.

The prior distributions are 
 \eqn{\beta \sim N(c,T)}{beta ~ N(c,T)} 
 and  
 \eqn{\rho \sim Uni(rmin,rmax)}{rho ~ Uni(rmin,rmax)} 
 or 
 \eqn{\rho \sim Beta(a1,a2)}{or rho ~ Beta(a1,a2)}.
}
\value{
Returns a structure of class \code{sarprobit}:
%\item{meth }{ 'sarprobit' }
\item{beta }{ posterior mean of bhat based on draws }
\item{rho }{ posterior mean of rho based on draws }
\item{bdraw }{ beta draws (ndraw-nomit x nvar) }
\item{pdraw }{ rho  draws (ndraw-nomit x 1) }
%\item{sdraw }{ sige draws (ndraw-nomit x 1) }
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
LeSage, J. and Pace, R. K. (2009), \emph{Introduction to Spatial Econometrics}, CRC Press, chapter 10 
}

\author{
adapted to and optimized for R by Stefan Wilhelm <wilhelm@financial.com> 
and Miguel Godinho de Matos <miguelgodinhomatos@cmu.edu>
based on code from James P. LeSage}

\seealso{
\code{\link{sar_lndet}} for computing log-determinants
}

\examples{
library(Matrix)
set.seed(2)

# number of observations
n <- 100

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
# W must not have non-zeros in the main diagonal!
i <- rep(1:n, each=m)
j <- rep(NA, n * m)
for (k in 1:n) {
  j[(((k-1)*m)+1):(k*m)] <- sample(x=(1:n)[-k], size=m, replace=FALSE)
}
W <- sparseMatrix(i, j, x=1/m, dims=c(n, n))


# innovations
eps <- rnorm(n=n, mean=0, sd=1)

# generate data from model 
S <- I_n - rho * W
z <- solve(qr(S), X \%*\% beta + eps)
y <- as.vector(z >= 0)  # 0 or 1, FALSE or TRUE

# estimate SAR probit model
sarprobit.fit1 <- sar_probit_mcmc(y, X, W, ndraw=100, thinning=1, prior=NULL)
summary(sarprobit.fit1)
}
