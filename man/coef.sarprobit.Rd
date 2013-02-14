\name{coef.sarprobit}
\alias{coef.sarprobit}
\alias{coefficients.sarprobit}
\alias{coef.semprobit}
\alias{coefficients.semprobit}
\alias{coef.sartobit}
\alias{coefficients.sartobit}

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
\method{coef}{sartobit}(object, ...)
\method{coefficients}{sartobit}(object, ...)
}
\arguments{
  \item{object}{class \code{sarprobit}, \code{semprobit} or \code{sartobit} with model fits}
  \item{...}{Additional arguments}
}
\value{
vector of named model coefficients
}
\author{
Stefan Wilhelm <wilhelm@financial.com>
}