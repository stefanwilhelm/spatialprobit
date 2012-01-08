\name{sar_probit_mcmc}
\alias{sar_probit_mcmc}
\title{Bayesian estimation of the spatial probit model}

\description{
Bayesian estimation of the spatial autoregressive probit model.
}

\usage{
sar_probit_mcmc(y, X, W, ndraw = 1000, burn.in = 100, thinning = 1, prior, start = list(rho = 0.75, beta = rep(0, ncol(X))))
}
\arguments{
  \item{y}{dependent variables. vector of zeros and ones}
  \item{X}{design matrix}
  \item{W}{spatial weight matrix}
  \item{ndraw}{number of MCMC iterations}
  \item{burn.in}{number of MCMC burn-in to be discarded}
  \item{thinning}{MCMC thinning factor, defaults to 1.}
  \item{prior}{
%%     ~~Describe \code{prior} here~~
}
  \item{start}{list of start values}
}
\details{
Bayesian estimates of the spatial autoregressive probit model
\deqn{
  y = \rho W y + X B + e, e \sim N(0,I_n)
  y is a binary 0,1 n x 1 vector
  B = N(c,T), 
  1/sige = Gamma(nu,d0), 
  \rho \sim Uniform(rmin,rmax), or \rho \sim beta(a1,a2)
}{
  y = rho*W*y + XB + e, e = N(0,I_n)
  y is a binary 0,1 n x 1 vector
  B = N(c,T), 
  1/sige = Gamma(nu,d0), 
  rho = Uniform(rmin,rmax), or rho = beta(a1,a2)
}
}
\value{
Returns a structure of class \code{sarprobit}:
\item{meth }{ 'sarp_g' }
\item{beta }{ posterior mean of bhat based on draws }
\item{rho }{ posterior mean of rho based on draws }
\item{sige }{ posterior mean of sige based on draws }
\item{sigma }{ posterior mean of sige based on (e'*e)/(n-k) }
\item{bdraw }{ bhat draws (ndraw-nomit x nvar) }
\item{pdraw }{ rho  draws (ndraw-nomit x 1) }
\item{sdraw }{ sige draws (ndraw-nomit x 1) }
\item{total }{ a matrix (ndraw,nvars-1) total x-impacts }
\item{direct }{ a matrix (ndraw,nvars-1) direct x-impacts }
\item{indirect }{ a matrix (ndraw,nvars-1) indirect x-impacts }
\item{total_obs }{ a matrix (ndraw,nvars-1) observation-level total x-impacts }
\item{vmean }{ mean of vi draws (nobs x 1) } 
\item{rdraw }{ r draws (ndraw-nomit x 1) (if m,k input) }
\item{bmean }{ b prior means, prior.beta from input }
\item{bstd }{ b prior std deviations sqrt(diag(prior.bcov)) }
\item{novi }{ 1 for prior.novi = 1, 0 for prior.rval input }
\item{nobs }{ # of observations }
\item{nvar }{ # of variables in x-matrix }
\item{ndraw }{ # of draws }
\item{nomit }{ # of initial draws omitted }
\item{nsteps }{ # of samples used by Gibbs sampler for TMVN }
\item{y }{ y-vector from input (nobs x 1) }
\item{zip }{ # of zero y-values }
\item{yhat }{ mean of posterior predicted (nobs x 1) }
\item{resid }{ residuals, based on posterior means }
\item{rsqr }{ r-squared based on posterior means }
\item{rbar }{ adjusted r-squared }
\item{a1 }{ a1 parameter for beta prior on rho from input, or default value }
\item{a2 }{ a2 parameter for beta prior on rho from input, or default value }
\item{time1 }{ time for eigenvalue calculation }
\item{time2 }{ time for log determinant calcluation }
\item{time3 }{ time for sampling }
\item{time }{ total time taken }  
\item{rmax }{ 1/max eigenvalue of W (or rmax if input) }
\item{rmin }{ 1/min eigenvalue of W (or rmin if input) }          
\item{tflag }{ 'plevel' (default) for printing p-levels; 
               'tstat' for printing bogus t-statistics }
\item{lflag }{ lflag from input }
\item{cflag }{  1 for intercept term, 0 for no intercept term }
\item{iter }{  prior.iter option from input }
\item{order }{  prior.order option from input }
\item{limit }{  matrix of [rho lower95,logdet approx, upper95],  
                intervals for the case of lflag = 1 }
\item{lndet }{ a matrix containing log-determinant information
                       (for use in later function calls to save time) }
\item{mlike }{ log marginal likelihood (a vector ranging over
                       rho values that can be integrated for model comparison) }
}
\references{
LeSage, J. and Pace, R. K. (2009), \emph{Introduction to Spatial Econometrics}, CRC Press, chapter 10 
}

\author{
Stefan Wilhelm <stefan.wilhelm@financial.com>, Miguel Godinho de Matos <miguelgodinhomatos@cmu.edu>
}

\seealso{
\code{\link{sar_lndet}}
}

\examples{
\dontrun{
library(Matrix)
n <- 200
beta <- c(0, 1, -1)
rho <- 0.75
# design matrix with two standard normal variates as "coordinates"
X <- cbind(intercept=1, x=rnorm(n), y=rnorm(n))
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
m <- 6
W <- sparseMatrix(i=rep(1:n, each=m), j=replicate(n, sample(x=1:n, size=m, replace=FALSE)), x=1/m, dims=c(n, n))
# identity matrix as spatial weight matrix for testing
#W <- I_n
  
# create samples from epsilon using independence of distributions (rnorm()) to avoid dense matrix I_n
eps <- rnorm(n=n, mean=0, sd=1)
 
S <- solve(I_n - rho * W)
z <- S \%*\% X \%*\% beta + S \%*\% eps
y <- as.vector(z >= 0)  # 0 or 1, FALSE or TRUE


# estimate model
fit1 <- sar_probit_mcmc(y, X, W, ndraw=500, burn.in=100, thinning=1, prior=NULL)
summary(fit1)
}
}
