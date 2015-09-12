\name{summary}
\alias{summary.sarprobit}
\alias{summary.sarorderedprobit}
\alias{summary.semprobit}
\alias{summary.sartobit}

\title{Print the results of the spatial probit/Tobit estimation via MCMC}
\description{
Print the results of the spatial probit/Tobit estimation via MCMC
}
\usage{
\method{summary}{sarprobit}(object, var_names = NULL, file = NULL, 
  digits = max(3, getOption("digits")-3), ...)
  
\method{summary}{sarprobit}(object, var_names = NULL, file = NULL, 
  digits = max(3, getOption("digits")-3), ...)  

\method{summary}{semprobit}(object, var_names = NULL, file = NULL, 
  digits = max(3, getOption("digits")-3), ...)

\method{summary}{sartobit}(object, var_names = NULL, file = NULL, 
  digits = max(3, getOption("digits")-3), ...)  
}
\arguments{
  \item{object}{class \code{sarprobit}, \code{semprobit} or \code{sartobit} or with model fits from 
    \code{\link{sarprobit}}, \code{\link{semprobit}} or \code{\link{sartobit}}}
  \item{var_names}{vector with names for the parameters under analysis}
  \item{file}{file name to be printed. If NULL or "" then print to console.}
  \item{digits}{integer, used for number formatting with signif() (for summary.default) or format() (for summary.data.frame).}
  \item{...}{Additional arguments}
}
\value{
This functions does not return any values.
}
\author{
Miguel Godinho de Matos <miguelgodinhomatos@cmu.edu>, Stefan Wilhelm <wilhelm@financial.com>
}

\seealso{
\code{\link{sarprobit}}, \code{\link{sarorderedprobit}}, \code{\link{semprobit}} or \code{\link{sartobit}} 
for SAR probit/SAR Ordered Probit/ SEM probit/ SAR Tobit model fitting
}