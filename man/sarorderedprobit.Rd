\name{sarorderedprobit}
\alias{sarorderedprobit}
\alias{sar_ordered_probit_mcmc}

\title{Bayesian estimation of the SAR ordered probit model}

\description{
Bayesian estimation of the spatial autoregressive ordered probit model (SAR ordered probit model).
}

\usage{
sarorderedprobit(formula, W, data, subset, ...)

sar_ordered_probit_mcmc(y, X, W, ndraw = 1000, burn.in = 100, thinning = 1, 
  prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0), 
  start = list(rho = 0.75, beta = rep(0, ncol(X)), 
  phi = c(-Inf, 0:(max(y)-1), Inf)),
  m=10, computeMarginalEffects=TRUE, showProgress=FALSE)

}
\arguments{
  \item{y}{dependent variables. vector of discrete choices from 1 to J (\{1,2,...,J\})}
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
  \item{computeMarginalEffects}{Flag if marginal effects are calculated. Defaults to TRUE. Currently without effect.}
  \item{showProgress}{Flag if progress bar should be shown. Defaults to FALSE.}
  \item{formula}{an object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{data}{an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. 
              If not found in data, the variables are taken from environment(formula), typically the environment from which sarprobit is called.}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{...}{additional arguments to be passed}
}
\details{
Bayesian estimates of the spatial autoregressive ordered probit model (SAR ordered probit model)
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
where y is a \eqn{(n \times 1)}{(n x 1)} vector of 
discrete choices from 1 to J, \eqn{y \in \{1,2,...,J\}}, where

  \eqn{y = 1}{y = 1}   for \eqn{-\infty  \le z \le \phi_1 = 0}{-Inf <= z <= phi_1 = 0} \cr
  \eqn{y = 2}{y = 2}   for \eqn{\phi_1   \le z \le \phi_2}{phi_1 <= z < phi_2} \cr
  ...     \cr
  \eqn{y = j}{y = j}   for \eqn{\phi_{j-1} \le z \le \phi_j}{phi_(j-1) <= z < phi_j} \cr
  ...     \cr
  \eqn{y = J}{y = J}   for \eqn{\phi_{J-1} \le z \le \infty}{phi_(J-1) <= z < Inf}

The vector \eqn{\phi=(\phi_1,...,\phi_{J-1})'}{phi} \eqn{(J-1 \times 1)}{(J - 1 x 1)} represents the cut points
(threshold values) that need to be estimated. \eqn{\phi_1 = 0}{phi_1 = 0} is set to zero by default.
  
\eqn{\beta}{beta} is a \eqn{(k \times 1)}{(k x 1)} 
vector of parameters associated with the \eqn{(n \times k)}{(n x k)} data matrix X.

\eqn{\rho}{rho} is the spatial dependence parameter.

The error variance \eqn{\sigma_e}{sigma_e} is set to 1 for identification.

Computation of marginal effects is currently not implemented.

}
\value{
Returns a structure of class \code{sarprobit}:
\item{beta }{ posterior mean of bhat based on draws }
\item{rho }{ posterior mean of rho based on draws }
\item{phi }{ posterior mean of phi based on draws }
\item{coefficients }{ posterior mean of coefficients }
\item{fitted.values }{ fitted values }
\item{fitted.reponse }{ fitted reponse }
\item{ndraw }{ \# of draws }
\item{bdraw }{ beta draws (ndraw-nomit x nvar) }
\item{pdraw }{ rho  draws (ndraw-nomit x 1) }
\item{phidraw }{ phi  draws (ndraw-nomit x 1) }
\item{a1 }{ a1 parameter for beta prior on rho from input, or default value }
\item{a2 }{ a2 parameter for beta prior on rho from input, or default value }
\item{time }{ total time taken }  
\item{rmax }{ 1/max eigenvalue of W (or rmax if input) }
\item{rmin }{ 1/min eigenvalue of W (or rmin if input) }
\item{tflag }{ 'plevel' (default) for printing p-levels; 
               'tstat' for printing bogus t-statistics }
\item{lflag }{ lflag from input }
\item{cflag }{ 1 for intercept term, 0 for no intercept term }
\item{lndet }{ a matrix containing log-determinant information
               (for use in later function calls to save time) }          
\item{W }{spatial weights matrix }
\item{X }{regressor matrix }
}

\references{
LeSage, J. and Pace, R. K. (2009), \emph{Introduction to Spatial Econometrics}, CRC Press, chapter 10, section 10.2 
}

\author{
Stefan Wilhelm <wilhelm@financial.com> 
}

\seealso{
  \code{\link{sarprobit}} for the SAR binary probit model
}
