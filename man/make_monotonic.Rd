% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_monotonic.R
\name{make_monotonic}
\alias{make_monotonic}
\title{Make a vector strictly monotonic}
\usage{
make_monotonic(vect, decreasing = T)
}
\arguments{
\item{vect}{A numeric vector.}

\item{decreasing}{Should the output indicate monotonically decreasing?}
}
\value{
A logical vector.
}
\description{
Given a vector \code{v} of length greater than 1, returns a logical vector
the same length as \code{v} indicating whether each element of \code{v} is
greater (or less-than, if \code{decreasing==TRUE}) than all previous elements.
}
