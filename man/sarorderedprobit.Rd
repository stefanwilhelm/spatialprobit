\name{sarorderedprobit}
\alias{sarorderedprobit}
\alias{sar_ordered_probit_mcmc}

\title{Bayesian estimation of the SAR probit model}

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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
