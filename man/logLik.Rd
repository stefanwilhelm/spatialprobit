\name{logLik}
\alias{logLik.sarprobit}
\alias{logLik.semprobit}
\title{Log Likelihood for spatial probit models (SAR probit, SEM probit)}

\description{
The functions return the log likelihood for the spatial autoregressive probit model
(SAR probit, spatial lag model) and the spatial error model probit (SEM probit).
}

\usage{
\method{logLik}{sarprobit}(object, ...)
\method{logLik}{semprobit}(object, ...)
}

\arguments{
  \item{object}{a fitted \code{sarprobit} or \code{semprobit} object}
  \item{\dots}{further arguments passed to or from other methods}
}

\value{
 returns an object of class \code{logLik}
}

\author{
Stefan Wilhelm <wilhelm@financial.com>
}

\seealso{
  \code{\link[spdep]{logLik.sarlm}}
}