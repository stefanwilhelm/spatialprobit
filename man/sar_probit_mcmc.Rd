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
  y = rho*W*y + XB + e, e = N(0,I_n)
  y is a binary 0,1 n x 1 vector
  B = N(c,T), 
  1/sige = Gamma(nu,d0), 
  rho = Uniform(rmin,rmax), or rho = beta(a1,a2);
}
\value{
Returns a structure of class \code{sarprobit}:
\item{results.meth }{ 'sarp_g' }
\item{results.beta }{ posterior mean of bhat based on draws }
\item{results.rho }{ posterior mean of rho based on draws }
\item{results.sige }{ posterior mean of sige based on draws }
\item{results.sigma }{ posterior mean of sige based on (e'*e)/(n-k) }
\item{results.bdraw }{ bhat draws (ndraw-nomit x nvar) }
\item{results.pdraw }{ rho  draws (ndraw-nomit x 1) }
\item{results.sdraw }{ sige draws (ndraw-nomit x 1) }
\item{results.total }{ a matrix (ndraw,nvars-1) total x-impacts }
\item{results.direct }{ a matrix (ndraw,nvars-1) direct x-impacts }
\item{results.indirect }{ a matrix (ndraw,nvars-1) indirect x-impacts }
\item{results.total_obs }{ a matrix (ndraw,nvars-1) observation-level total x-impacts }
\item{results.vmean }{ mean of vi draws (nobs x 1) } 
\item{results.rdraw }{ r draws (ndraw-nomit x 1) (if m,k input) }
\item{results.bmean }{ b prior means, prior.beta from input }
\item{results.bstd }{ b prior std deviations sqrt(diag(prior.bcov)) }
\item{results.novi }{ 1 for prior.novi = 1, 0 for prior.rval input }
\item{results.nobs }{ # of observations }
\item{results.nvar }{ # of variables in x-matrix }
\item{results.ndraw }{ # of draws }
\item{results.nomit }{ # of initial draws omitted }
\item{results.nsteps }{ # of samples used by Gibbs sampler for TMVN }
\item{results.y }{ y-vector from input (nobs x 1) }
\item{results.zip }{ # of zero y-values }
\item{results.yhat }{ mean of posterior predicted (nobs x 1) }
\item{results.resid }{ residuals, based on posterior means }
\item{results.rsqr }{ r-squared based on posterior means }
\item{results.rbar }{ adjusted r-squared }
\item{results.a1 }{ a1 parameter for beta prior on rho from input, or default value }
\item{results.a2 }{ a2 parameter for beta prior on rho from input, or default value }
\item{results.time1 }{ time for eigenvalue calculation }
\item{results.time2 }{ time for log determinant calcluation }
\item{results.time3 }{ time for sampling }
\item{results.time }{ total time taken }  
\item{results.rmax }{ 1/max eigenvalue of W (or rmax if input) }
\item{results.rmin }{ 1/min eigenvalue of W (or rmin if input) }          
\item{results.tflag }{ 'plevel' (default) for printing p-levels; 
                       'tstat' for printing bogus t-statistics }
\item{results.lflag }{ lflag from input }
\item{results.cflag }{  1 for intercept term, 0 for no intercept term }
\item{results.iter }{  prior.iter option from input }
\item{results.order }{  prior.order option from input }
\item{results.limit }{  matrix of [rho lower95,logdet approx, upper95],  
                        intervals for the case of lflag = 1 }
\item{results.lndet }{ a matrix containing log-determinant information
                       (for use in later function calls to save time) }
\item{results.mlike }{ log marginal likelihood (a vector ranging over
                       rho values that can be integrated for model comparison) }
}
\references{
John LeSage. Introduction to Spatial Econometrics. chapter 10. 2009
}
\author{
Stefan Wilhelm <stefan.wilhelm@financial.com>, Miguel Godinho de Matos
}

\seealso{
\code{\link{sar_lndet}}
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
