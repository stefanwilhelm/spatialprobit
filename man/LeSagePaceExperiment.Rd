\name{LeSagePaceExperiment}
\alias{LeSagePaceExperiment}
\title{Replicate the LeSage and Pace (2009), section 10.1.5 experiment}
\description{
 This method replicates the experiment from LeSage and Pace (2009), section 10.1.5.
 It first generates data from a SAR probit model and then estimates the model with our implementation. 
}
\usage{
LeSagePaceExperiment(n = 400, beta = c(0, 1, -1), rho = 0.75, ndraw = 1000, 
 burn.in = 200, thinning = 1, m = 10, computeMarginalEffects=TRUE, ...)
}
\arguments{
  \item{n}{sample size}
  \item{beta}{parameter vector}
  \item{rho}{spatial dependence parameter}
  \item{ndraw}{number of draws}
  \item{burn.in}{number of burn-in samples}
  \item{thinning}{thinning parameter}
  \item{m}{Gibbs sampler burn-in size for drawing from the truncated multinormal distribution}
  \item{computeMarginalEffects}{Should marginal effects be computed?}
  \item{...}{Additional parameters to be passed to \code{\link{sar_probit_mcmc}}}
}
\value{
 Returns a structure of class \code{sarprobit}
}
\references{
 LeSage, J. and Pace, R. K. (2009), \emph{Introduction to Spatial Econometrics}, CRC Press, section 10.1.5
}

\author{
Stefan Wilhelm <wilhelm@financial.com>
}

\examples{
\dontrun{
# LeSage/Pace(2009), Table 10.1, p.291: n=400, m=10
res1 <- LeSagePaceExperiment(n=400, beta=c(0,1,-1), rho=0.75, 
 ndraw=1000, burn.in=200, thinning=1, m=10)
res1$time
res1$coefficients
summary(res1)

# LeSage/Pace(2009), Table 10.1, p.291: n=1000, m=1
res2 <- LeSagePaceExperiment(n=1000, beta=c(0,1,-1), rho=0.75, 
  ndraw=1000, burn.in=200, thinning=1, m=1)
res2$time
res2$coefficients
summary(res2)

# LeSage/Pace(2009), Table 10.2, p.291: n=400, m=1
res400.1 <- LeSagePaceExperiment(n=400, beta=c(0,1,-1), rho=0.75, 
  ndraw=1000, burn.in=200, thinning=1, m=1)
summary(res400.1)

# LeSage/Pace(2009), Table 10.2, p.291: n=400, m=2
res400.2 <- LeSagePaceExperiment(n=400, beta=c(0,1,-1), rho=0.75, 
  ndraw=1000, burn.in=200, thinning=1, m=2)
summary(res400.2)

# LeSage/Pace(2009), Table 10.2, p.291: n=400, m=10
res400.10 <- LeSagePaceExperiment(n=400, beta=c(0,1,-1), rho=0.75, 
  ndraw=1000, burn.in=200, thinning=1, m=10)
summary(res400.10)
}
}
