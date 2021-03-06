% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/assumptionChecking.R
\name{plausibleUniformProbabilities}
\alias{plausibleUniformProbabilities}
\title{Assess potential uniformity}
\usage{
plausibleUniformProbabilities(field)
}
\arguments{
\item{field}{The \code{probField} of interest}
}
\value{
An \code{h-test} object. See return of \link{ks.test}
}
\description{
Perform a Kolmogorov-Smirnov test to see if the jackknife probabilities
of survey nodes in a probability field could have been generated from a
uniform distribution
}

