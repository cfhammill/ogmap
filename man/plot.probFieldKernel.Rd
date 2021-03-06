% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/distanceAndKernel.R
\name{plot.probFieldKernel}
\alias{plot.probFieldKernel}
\title{Visualize a Kernel}
\usage{
\method{plot}{probFieldKernel}(x, domain, domainMax = NULL, nvals = 1000,
  type = "l", ...)
}
\arguments{
\item{x}{A \code{probFieldKernel} of interest}

\item{domain}{Either a vector of values to evaluate the kernel at, or the minimum value for
a sequence}

\item{domainMax}{Optional maximum value of the sequence}

\item{nvals}{The number of sequence values to generate (if xmax is supplied)}

\item{type}{A plot parameter indicating which type of plot to draw, defaults to line
see \link{par} for details}

\item{...}{Additional plot parameters to be passed to plot}
}
\description{
Create a simple plot showing the value of a kernel evaluated over a range of distances
}
\details{
If xmax is supplied, a sequence of 1000 values
}

