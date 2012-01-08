\name{coef.sarprobit}
\alias{coef.sarprobit}
\alias{coefficients.sarprobit}

\title{Extract Model Coefficients}
\description{
coef is a generic function which extracts model coefficients from objects returned by modeling functions. coefficients is an alias for it.
}
\usage{
\method{coef}{sarprobit}(object, ...)
\method{coefficients}{sarprobit}(object, ...)
}
\arguments{
  \item{object}{class \code{sarprobit} with the results from \code{sar_probit_mcmc}}
  \item{...}{Additional arguments}
}
\value{
This functions does not return any values.
}
\author{
Stefan Wilhelm
}