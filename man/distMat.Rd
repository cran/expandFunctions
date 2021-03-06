% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expandFunctions.R
\name{distMat}
\alias{distMat}
\title{Make a matrix with coefficients distributed as dist}
\usage{
distMat(p, q, dist = rnorm, distOpt = NULL)
}
\arguments{
\item{p}{Number of rows}

\item{q}{Number of columns}

\item{dist}{distribution of coefficients to draw from;
default is rnorm.}

\item{distOpt}{Named list of additional parameters for dist.
\emph{Always omit the first parameter,n, of the
distribution sampling function}. Defaults may
be omitted if desired (see examples).}
}
\value{
A pXq matrix with coefficients distributed
        as dist with parameters defined in `...`
}
\description{
Generate a pXq matrix with coefficients drawn
             from the univariate distribution dist with
             options distOpt.
}
\details{
The user may provide their own distribution function,
         but note that the number of values to return, n,
         \emph{must} be the first argument, just as with
         the built-in distributions.  The first argument does
         not have to be named.
}
\examples{
X <- distMat(10,2)
X <- distMat(10,2,distOpt=list(mean=1,sd=2))
X <- distMat(5,3,rnorm,list(mean=1,sd=2))
X <- distMat(5,3,rnorm,list(sd=2))
X <- distMat(50,3,rt,list(df=3))
}

