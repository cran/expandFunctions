% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expandFunctions.R
\name{eReplace}
\alias{eReplace}
\title{Replace values in an R object coerible to a matrix}
\usage{
eReplace(X, a, b)
}
\arguments{
\item{X}{R object coercible to a matrix}

\item{a}{Value to be replaced}

\item{b}{Value to replace}
}
\value{
X with all a's replaced with b's.  a may be NA
}
\description{
Replace values in an R object coerible to a matrix.
  It is useful for replacing NA with other values, etc.,
  such as with padding values.
}
\examples{
A <- matrix(1:6,ncol=2)
A <- eReplace(A,1,NA)
A <- eReplace(A,NA,-9999)
A <- eReplace(A,-9999,0)
}
\seealso{
\code{\link[base]{replace}}, which performs the same
operation on vectors, and on which this operation is based.
}

