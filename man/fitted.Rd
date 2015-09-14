\name{fitted}
\alias{fitted.sarprobit}
\alias{fitted.sarorderedprobit}
\alias{fitted.semprobit}
\alias{fitted.sartobit}
\title{Fitted values of spatial probit/Tobit models}

\description{
Calculate fitted values of spatial probit/Tobit models. 
}

\usage{
\method{fitted}{sarprobit}(object, ...)
\method{fitted}{sarorderedprobit}(object, ...)
\method{fitted}{semprobit}(object, ...)
\method{fitted}{sartobit}(object, ...)
}

\arguments{
  \item{object}{a fitted model of class \code{sarprobit},\code{sarorderedprobit},\code{semprobit} or \code{sartobit}}
  \item{\dots}{further arguments passed to or from other methods}
}

\value{
 A numeric vector of the fitted values. 
}

\author{
Stefan Wilhelm <wilhelm@financial.com>
}