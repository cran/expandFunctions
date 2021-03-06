% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expandFunctions.R
\name{lagshift}
\alias{lagshift}
\title{Helper function for eLag.}
\usage{
lagshift(x, i, lagMax, pad)
}
\arguments{
\item{x}{Input vector}

\item{i}{Shift (integer)}

\item{lagMax}{Maximum lag that will be needed}

\item{pad}{Scalar used for padding.}
}
\value{
vector padded front and back with padding appropriate
for generating lag.
}
\description{
Generates shifted columns.
}
\examples{
lagshift(1:3,0,1,NA)
}

