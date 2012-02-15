\name{summary.sarprobit}
\alias{summary.sarprobit}
\title{Print the results of the SAR probit estimation via MCMC}
\description{
Print the results of the SAR probit estimation via MCMC
}
\usage{
\method{summary}{sarprobit}(object, var_names = NULL, file = NULL, 
  digits = max(3, getOption("digits")-3), ...)
}
\arguments{
  \item{object}{class \code{sarprobit} with the results from \code{sar_probit_mcmc}}
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

\examples{
# TO BE DONE
}