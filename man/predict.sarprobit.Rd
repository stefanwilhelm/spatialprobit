

\name{predict.sarprobit}
\alias{predict.sarprobit}
\title{Predict method for SAR Probit Model Fits}
\description{
Predicted values based on SAR probit model object.
}
\usage{
\method{predict}{sarprobit}(object, newdata, na.action, ...)
}
\arguments{
  \item{object}{Object class inheriting from \code{sarprobit}}
  \item{newdata}{An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.}
  \item{na.action}{function determining what should be done with missing values in \code{newdata}. The default is to predict \code{NA}.}
  \item{...}{Additional arguments}
}
\value{
This functions does not return any values.
}
\author{
Stefan Wilhelm
}