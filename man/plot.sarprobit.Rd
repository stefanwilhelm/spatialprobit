\name{plot.sarprobit}
\alias{plot.sarprobit}
\title{
Plot Diagnostics for an \code{sarprobit} Object
}
\description{
Three plots (selectable by which) are currently available:
MCMC trace plots, autocorrelation plots and posterior density plots
}
\usage{
\method{plot}{sarprobit}(x, 
which = c(1, 2, 3), 
ask = prod(par("mfcol")) < length(which) && dev.interactive(),
...,
trueparam = NULL)
}

\arguments{
  \item{x}{a \code{sarprobit} Object}
  \item{which}{if a subset of the plots is required, specify a subset of the numbers 1:3.}
  \item{ask}{logical; if TRUE, the user is asked before each plot, see \code{par(ask=.)}.}
  \item{...}{other parameters to be passed through to plotting functions.}
  \item{trueparam}{...}
}

\author{
Stefan Wilhelm <stefan.wilhelm@financial.com>
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, which = c(1, 2, 3), trueparam = NULL) 
{
    names <- x$names
    B <- x$B
    k <- ncol(B)
    show <- rep(FALSE, 3)
    show[which] <- TRUE
    if (show[1L]) {
        for (i in 1:k) {
            plot(1:nrow(B), B[, i], type = "l", main = substitute("Trace plot of " * 
                x, list(x = names[i])))
            if (!is.null(trueparam)) 
                abline(h = trueparam[i], col = "red", lty = 2)
        }
    }
    if (show[2L]) {
        for (i in 1:k) {
            acf(B[, i], main = substitute("ACF of " * x, list(x = names[i])))
        }
    }
    if (show[3L]) {
        for (i in 1:k) {
            plot(density(B[, i]), main = substitute("Posterior distribution of " * 
                x, list(x = names[i])))
            if (!is.null(trueparam)) 
                abline(v = trueparam[i], col = "red", lty = 2)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
