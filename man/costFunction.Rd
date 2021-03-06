% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/assumptionChecking.R
\name{costFunction}
\alias{costFunction}
\title{Evaluate fit of a field}
\usage{
costFunction(field, bandwidths, rawPairs, approximate = FALSE,
  verbose = FALSE)
}
\arguments{
\item{field}{The \code{probField} of interest}

\item{bandwidths}{A vector or list containing the new kernel and distance parameters to asses the cost at.
See \link{setBandwidths} and \link{getBandwidths} for more details}

\item{rawPairs}{Pairs of interest to be passed down to \link{plausibleIndependentPairs}}

\item{approximate}{Logical whether or not to approximate the cost}

\item{verbose}{Logical whether or not to print the bandwidth and the cost to the console}
}
\value{
A numeric value equal to negative log the value from \link{plausibleUniformIndependent} or
and error code of 1000. Or in if using the approximation, return two times the product of the
KS statistics.
}
\description{
Evaluate the fit of a field given a set of new bandwidths.
}

