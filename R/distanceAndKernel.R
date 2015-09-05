########################################################################
# distanceAndKernel.R
# Copyright 2015 Geoff Evans
#
# This file is part of the Ogmap R Package.
#
# Ogmap is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Ogmap is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Ogmap.  If not, see <http://www.gnu.org/licenses/>.
########################################################################


#' @title Create an inter-node distance function
#' @description Create a distance function to be used in calculating the distance between two
#' \code{probFieldNodes}
#' @param fun A function of two nodes
#' @param ... Addition scale parameters used in the calculation of distance
#' @param guides An optional named list of references (either named or positional) to the guide members
#' that are used in the distance calculation
#' @return A \code{probFieldDistance} object, which is a wrapped function with attributes
#' representing scale parameters
#' @details
#' This function constructs a distance metric of the users choosing to be used in calculations
#' within and involving probability fields. Run \code{linearDistance} (without brackets) and
#' arguments to see a simple sample implementation. For compatibility the function must take
#' two arguments (each a \code{probFieldNode}) extract the guide variables and calculate a
#' meaningful distance metric to be used by \code{probFieldKernel} objects. For a complex example
#' see the code for \link{trawlDistance}.
#' @aliases probFieldDistance
#' @export
distance <- function(fun, ..., guides = NULL){
  distanceFull <- structure(fun,
                            distanceParameters = list(...),
                            guides = guides,
                            call = match.call(),
                            class = c("probFieldDistance","list"))
  distanceFull
}

#' @title Retrieve Distance Parameters
#' @description Extract the additional parameters used by a distance function
#' @param distance A \code{probFieldDistance}
#' @return A named vector containing additional parameters and their values
#' @details getScale is an alias to getDistanceParameters, the name of this function
#' derives from the scale parameter used in \link{linearDistance}
#' @export
getDistanceParameters <- function(distance){
  unlist(attr(distance, "distanceParameters"))
}

#'@rdname getDistanceParameters
#'@export
getScale <- getDistanceParameters


#' @title Linear Distance Between Nodes
#' @description Generates a \code{probFieldDistance} object that calculates the absolute linear
#' distance between two nodes divided by a given scale
#' @param scale The scaling parameter
#' @param guides A vector containing the names of the guide variables used in the distance calculation,
#' must be provided for multidimensional guide variables
#' @return A \code{probFieldDistance} object which wraps a function that extracts the guide
#' values from two given nodes and returns the absolute linear distance
#' @details
#' This code is designed to make a \code{probFieldDistance} object that will calculate a
#' linear distance between the guides of two nodes. It performs pair-wise subtractions
#' and element-wise division, the distances in each guide variable are summed
#' @export
linearDistance <- function(scale, guides = NULL) {

  distance(function(node1, node2){

    #ensure scale is non-zero
    if(scale == 0){
      scale <- 1
    } else {
      scale <- abs(scale)
    }

    g1 <- unlist(getGuide(node1))
    g2 <- unlist(getGuide(node2))

    if(!is.null(guides)){
      g1 <- g1[unlist(guides)]
      g2 <- g2[unlist(guides)]
    }

    sum(abs(g1 - g2) / scale)
    },
    scale = scale, guides = guides)
}

#' @title Distance Between Trawl Surveys
#' @description
#' Generates a \code{probFieldDistance} object that calculates the scaled distance between two trawl survey nodes
#' given latitude, longitude, and depth, with respect to independent horizontal and vertical scales.
#' @param horizontalScale The scaling to be applied to the between lat/long distance
#' @param verticalScale The scaling to be applied to the depth distance
#' @param latGuide The name of the guide parameter containing the latitudinal coordinate
#' within the \code{probFieldNode}s
#' @param longGuide The name of the guide parameter containing the longitudinal coordinate
#' within the \code{probFieldNode}s
#' @param verticalGuide The name of the guide parameter containing the depth coordinate
#' within the \code{probFieldNode}s
#' @param earthRadius The radius of the earth in kilometers at the approximate latitude of interest
#' defaults to the radius of the earth at 55N
#' @details
#' This function was included for historical reasons and as a sample of a more complex distance function.
#' Internally, the function generated creates and evaluates both an \link{earthDistance} and a \link{linearDistance}
#' The \code{earthDistance} evaluates the scaled earth-surface inter-lat/long coordinate distance given the radius
#' of the earth at your latitude of interest, the default is the earth's radius at 55N. The \link{linearDistance}
#' evaluates the scaled distance between node depths. The results of evaluating these two functions is summed to
#' give the overall scaled distance between two nodes. None of the actual computation is performed when this function is
#' called, instead a function of two nodes is created that will perform all of the computation on demand. See the documentation
#' for \link{probFieldDistance} for more details.
#' @export
trawlDistance <- function(horizontalScale, verticalScale,
                          latGuide, longGuide, verticalGuide,
                          earthRadius = 6367){
  distance(
    function(node1, node2){
      g1 <- getGuide(node1)
      g2 <- getGuide(node2)
      hdist <- earthDistance(horizontalScale, latGuide, longGuide)(node1, node2)
      vdist <- linearDistance(verticalScale, verticalGuide)(node1, node2)
      as.numeric(hdist + vdist)
    }, horizontalScale = horizontalScale, verticalScale = verticalScale,
    guides = list(latGuide = latGuide, longGuide = longGuide, verticalGuide = verticalGuide)
    )
}

#' @title Earth Surface Distance
#' @description Create \code{probFieldDistance} object to calculate the distance between two nodes
#' along the surface of the earth, adjusting for curvature, given lat and long coordinates.
#' @param horizontalScale The scaling factor to apply to the distance
#' @param latGuide The name of the guide parameter containing the latitudinal coordinate
#' within the \code{probFieldNode}s
#' @param longGuide The name of the guide parameter containing the longitudinal coordinate
#' within the \code{probFieldNode}s
#' @param earthRadius The radius of the earth in kilometers at the approximate latitude of interest
#' defaults to the radius of the earth at 55N
#' @export
earthDistance <- function(horizontalScale, latGuide, longGuide, earthRadius = 6367){
  degreesPerRad = pi/180
  distance(
    function(node1, node2){

      g1 <- getGuide(node1)
      g1Lat <- g1[latGuide]
      g1Long <- g1[longGuide]

      g2 <- getGuide(node2)
      g2Lat <- g2[latGuide]
      g2Long <- g2[longGuide]

      latDist <- g1Lat - g2Lat
      longDist <- (g1Long - g2Long) * cos(degreesPerRad * (g1Lat + g2Lat)/2)

      dist <- earthRadius * degreesPerRad * sqrt(latDist^2 + longDist^2) / horizontalScale

      as.numeric(dist)
    },
    horizontalScale = horizontalScale,
    guides = list(latGuide = latGuide, longGuide = longGuide))
}

#' @title Create a kernel weighting function
#' @description Create a custom \code{probFieldKernel} object to evaluate the relevance of an
#' observation at a given distance
#' @param fun A function of distance
#' @param ... Extra shape parameters
#' @return A \code{probFieldKernel} object that wraps a kernel function with an additional
#' attribute containing the shaping parameters, and the call used to produce the kernel function
#' @details
#' The function allows a user to create a custom kernel function to assign a weighting to an
#' observation at a given distance. The provided function needs only be a function of a number
#' representing a distance.
#' @aliases probFieldKernel
#' @export
kernel <- function(fun, ...){
  kernelFull <- structure(fun,
                          kernelParameters = list(...),
                          call = match.call(),
                          class = c("probFieldKernel", "list"))
  kernelFull
}

#' @title Kernel Parameters
#' @description Return a vector of the additional shape parameters used in a kernel function
#' @param kernel The \code{probFieldKernel} of interest
#' @return A named vector containing the additional shape parameters used in the kernel function
#' @details \code{getShape} is an alias to \code{getKernelParameters} this function's name derives
#' from the single shape parameter necessary for many kernel functions. See the implementation of
#' \link{subcauchy} for an example
#' @export
getKernelParameters <- function(kernel){
  unlist(attr(kernel, "kernelParameters"))
}

#' @rdname getKernelParameters
#' @export
getShape <- getKernelParameters

#' @title Subcauchy Kernel
#' @description Generate a \code{probFieldKernel} object used to evaluate the relevance of an observation at
#' a given distance using a subcauchy kernel governed by the parameter shape
#' @param shape The shape parameter governing the subcauchy kernel
#' @return A \code{probFieldKernel} object which wraps a subcauchy kernel function with the shape
#' exponent \code{shape} and an additional attribute containing the value of \code{shape}
#' @export
subcauchy <- function(shape){
  kernel(function(d) 1 / (1 + d^shape),
         shape = shape)
}

#' @title Logistic Kernel
#' @description
#' Generate a \code{probFieldKernel} object to evaluate the relevance of an observation at a give distance
#' using a logistic kernel governed by a single shape parameter
#' @param shape The shape parameter governing the logistic kernel
#' @return  A \code{probFieldKernel} object which wraps a logistic kernel function with the shape
#' parameter \code{shape} and an additional attribute containing the value of \code{shape}
logist <- function(shape){
  kernel(function(d){
    (shape + 1) / (shape + (shape + 2)^d)
  },
  shape = shape)
}

#' Debug Kernel
#' A simple kernel for debugging
#' @return A \code{probFieldKernel} wrapping a function of distance which returns 0.5 if distance is less than
#' 1, 0.5 if distance is between .5 and 1.5 exclusive, and .25 if distance is less than 2.5 and zero for any other
#' input
#' @export
debugKernel <- function(){
  kernel(function(d){
    if(d < 0.5) return(1)
    if(d < 1.5) return(0.5)
    if(d < 2.5) return(.25)
    return(0)
  })
}

#' @title Visualize a Kernel
#' @description
#' Create a simple plot showing the value of a kernel evaluated over a range of distances
#' @param x A \code{probFieldKernel} of interest
#' @param domain Either a vector of values to evaluate the kernel at, or the minimum value for
#' a sequence
#' @param domainMax Optional maximum value of the sequence
#' @param nvals The number of sequence values to generate (if xmax is supplied)
#' @param type A plot parameter indicating which type of plot to draw, defaults to line
#' see \link{par} for details
#' @param ... Additional plot parameters to be passed to plot
#' @return NULL
#' @details
#' If xmax is supplied, a sequence of 1000 values
#' @export
plot.probFieldKernel <- function(x, domain, domainMax=NULL, nvals = 1000, type = "l", ...){
  if(!is.null(domainMax)) domain <- seq(domain, domainMax, length.out = nvals)
  plot(domain, x(domain), type = type, ...)
}

