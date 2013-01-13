\name{coef.sarprobit}
\alias{coef.sarprobit}
\alias{coefficients.sarprobit}
\alias{coef.semprobit}
\alias{coefficients.semprobit}

\title{Extract Model Coefficients}
\description{
coef is a generic function which extracts model coefficients 
from objects returned by modeling functions. coefficients is an alias for it.
}
\usage{
\method{coef}{sarprobit}(object, ...)
\method{coefficients}{sarprobit}(object, ...)
\method{coef}{semprobit}(object, ...)
\method{coefficients}{semprobit}(object, ...)
}
\arguments{
  \item{object}{class \code{sarprobit} or \code{semprobit} with the results from 
    \code{sarprobit()}/\code{sar_probit_mcmc()} or \code{semprobit()}/\code{sem_probit_mcmc()}}
  \item{...}{Additional arguments}
}
\value{
vector of named model coefficients
}
\author{
Stefan Wilhelm <wilhelm@financial.com>
}