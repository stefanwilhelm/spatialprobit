\name{summary.sarprobit}
\alias{summary.sarprobit}
\title{Print the results of the SAR probit estimation via MCMC}
\description{
Print the results of the SAR probit estimation via MCMC
}
\usage{
\method{summary}{sarprobit}(object, var_names = NULL, file = NULL, ...)
}
\arguments{
  \item{object}{class \code{sarprobit} with the results from \code{sar_probit_mcmc}}
  \item{var_names}{vector with names for the parameters under analysis}
  \item{file}{file name to be printed. If NULL or "" then print to console.}
  \item{...}{Additional arguments}
}
\value{
This functions does not return any values.
}
\author{
Miguel Godinho de Matos <miguelgodinhomatos@cmu.edu>
}

\examples{
# TO BE DONE
}