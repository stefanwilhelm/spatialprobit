\name{sar_lndet}
\alias{sar_lndet}
\title{
compute the log determinant |I_n - rho*W| using the user-selected (or default) method
}
\description{
compute the log determinant |I_n - rho*W| using the user-selected (or default) method
}
\usage{
sar_lndet(ldetflag, W, rmin, rmax)
}
\arguments{
  \item{ldetflag}{
%%     ~~Describe \code{ldetflag} here~~
}
  \item{W}{spatial weight matrix}
  \item{rmin}{minimum eigen value}
  \item{rmax}{maximum eigen value}
}

\value{
\item{detval}{}
\item{time}{execution time}
}

\references{
John LeSage. Introduction to Spatial Econometrics. 2009
}

\author{
James P. LeSage, Adapted to R by Miguel Godinho de Matos <miguelgodinhomatos@cmu.edu>
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
