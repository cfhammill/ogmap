% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/assumptionChecking.R
\name{plausibleUniformIndependent}
\alias{plausibleUniformIndependent}
\title{Assess field assumptions}
\usage{
plausibleUniformIndependent(field, rawPairs)
}
\arguments{
\item{field}{The \code{probField} of interest}

\item{rawPairs}{The pairs of interest to be passed to \link{plausibleIndependentPairs}}
}
\value{
a numeric value equal to the product of the two p-values (see \link{plausibleUniformProbabilities}
 and \link{plausibleIndependentPairs}) divided by 1 minus the log that product.
}
\description{
Assess the two major assumptions made for field construction - pair independence and uniformity
}

