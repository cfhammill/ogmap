% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/statistic.R
\name{probFieldStatistic}
\alias{probFieldStatistic}
\title{Statistic}
\usage{
probFieldStatistic(field, nodes, designFunction, resultFunction, ...)
}
\arguments{
\item{field}{The \code{probField} of interest}

\item{nodes}{The \code{probFieldNodeList} denoting where to evaluate the field}

\item{designFunction}{A special function for generating the statistical requirements for generating inference.
See details for requirements of the function}

\item{resultFunction}{A special function used on the results of the \code{designFunction} to generate
inference. See details for requirements of the function}

\item{...}{Additional named parameters necessary for the \code{designFunction}}
}
\value{
An object of class \code{probFieldStatistic} containing:
\itemize{
\item{\code{observed}:}{ The observed values from the field}
\item{\code{design}:}{ The results of the design function}
\item{\code{result}:}{ The result function such that it can be called by the user}
}
}
\description{
Create a statistic object for a \code{probField} and a \code{probFieldNodeList} denoting where to evaluate the field.
The statistic can be used in resampling procedures for putting confidence limits on parameter estimates for a
field.
}
\details{
The \code{probFieldStatistic} class is concerned with generating inference from a \code{probField} at some
internal portion of its domain. Statistics of these type are frequently used in resampling routines so the
class is designed to facillitate this. The \code{designFunction} is a special function used to calculate the
statistical machinery necessary to calculate the result such that it needs only be calculated once yet provide
the capability of rapidly resampling many times. The \code{resultFunction} is required to take the machinery
produced by the \code{designFunction} and produce a value for the statistic. The only major constraint on these two
functions is that they interoperate seamlessly. The \code{resultFunction} must take as its first argument the
result of the corresponding \code{designFunction}. Additionally, it would be ideal to allow an indices argument
within the \code{resultFunction} for compatibility with resampling routines.
}

