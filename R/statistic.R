########################################################################
# statistic.R
# Copyright Her Majesty the Queen in right of Canada, 2015
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


#' @title Statistic
#' @description
#' Create a statistic object for a \code{probField} and a \code{probFieldNodeList} denoting where to evaluate the field.
#' The statistic can be used in resampling procedures for putting confidence limits on parameter estimates for a
#' field.
#' @param field The \code{probField} of interest
#' @param nodes The \code{probFieldNodeList} denoting where to evaluate the field
#' @param designFunction A special function for generating the statistical requirements for generating inference.
#' See details for requirements of the function
#' @param resultFunction A special function used on the results of the \code{designFunction} to generate
#' inference. See details for requirements of the function
#' @param ... Additional named parameters necessary for the \code{designFunction}
#' @return An object of class \code{probFieldStatistic} containing:
#' \itemize{
#' \item{\code{observed}:}{ The observed values from the field}
#' \item{\code{design}:}{ The results of the design function}
#' \item{\code{result}:}{ The result function such that it can be called by the user}
#' }
#' @details
#' The \code{probFieldStatistic} class is concerned with generating inference from a \code{probField} at some
#' internal portion of its domain. Statistics of these type are frequently used in resampling routines so the
#' class is designed to facillitate this. The \code{designFunction} is a special function used to calculate the
#' statistical machinery necessary to calculate the result such that it needs only be calculated once yet provide
#' the capability of rapidly resampling many times. The \code{resultFunction} is required to take the machinery
#' produced by the \code{designFunction} and produce a value for the statistic. The only major constraint on these two
#' functions is that they interoperate seamlessly. The \code{resultFunction} must take as its first argument the
#' result of the corresponding \code{designFunction}. Additionally, it would be ideal to allow an indices argument
#' within the \code{resultFunction} for compatibility with resampling routines.
#' @export
probFieldStatistic <- function(field, nodes, designFunction, resultFunction, ...){
  design <- designFunction(field, nodes, ...)

  structure(list(
    observed = field$modifiedData[,field$response],
    design = design,
    result = function(...){resultFunction(designResult = design, ...)}
  ), class = c("probFieldStatistic", "list"))
}

#' @title Calculate Expected Values
#' @description
#' Generate a \code{probFieldStatistic} to provide easy resampling of expected values for one or many nodes
#' @param field The \code{probField} of interest
#' @param nodes Either a \code{probFieldNode} or \code{probFieldNodeList} of nodes for which to calculate the
#' expectation
#' @return
#' A \code{probFieldStatistic} object containing the observed values from the survey nodes in the field,
#' a callable result function to allow statistic resampling, and the results of applying the design function
#' to the field and nodes provided, in this case a matrix of steps probabilities for the nodes provided and
#' the values from the field.
#' @export
expectationStat <- function(field, nodes){

  designFunction <- function(field, nodes){
    nodes <- evaluateNodes(nodes, field)
    values <- getVariable(field$survey)
    steps <- getSteps(nodes)
    list(steps = steps, values = values)
  }

  resultFunction <- function(designResult, indices = NULL){
    if(is.null(indices)) indices <- 1:length(designResult$values)
    resultMatrix <- as.matrix(designResult$values[indices] * designResult$steps)
    colSums(resultMatrix)
  }

  probFieldStatistic(field, nodes, designFunction, resultFunction)
}

#' @title Weighted Sum of Means
#' @description
#' Calculated a weighted sum of means for a group of nodes within the domain of a probability field
#' @param field The \code{probField} of a interest
#' @param nodes Ideally a \code{probFieldNodeList}, but generally any list of \code{probFieldNode} objects indicating
#' where in the domain of \code{field} to generate inference
#' @param weights A vector of weights (0-1) for each node in \code{nodes} to weight the statistic
#' @return A \link{probFieldStatistic} object containing the observations from the field, the design results (in this case
#' a list containing the variance matrix and the weighted step probabilities from the field), and a function obtain and
#' resample the statistic. See \code{probFieldStatistic} for details of the result function.
#' @export
weightedSumOfMeans <- function(field, nodes, weights){
  designFun <- function(field, nodes, weights){
    nodesWithCDFs <- evaluateNodes(nodes, field)
    steps <- as.matrix(getSteps(nodesWithCDFs))
    survWeights <- (steps %*% weights)[,1]

    varmat <- matrix(0, nrow = length(field$survey), ncol = length(field$survey))

    for(i in 1:length(field$survey)){
      n <- field$survey[[i]]
      steps <- n$steps
      var <- diag(steps) - outer(steps,steps)
      varmat <- varmat + survWeights[i]^2 * var
    }

    list(varianceMatrix = varmat, surveyWeights = survWeights, observations = getVariable(field))
  }

  resultFun <- function(designResult, indices = NULL, verbose = FALSE){
    if(is.null(indices)) indices <- 1:length(designResult$surveyWeights)
    stat <- sum(designResult$observations[indices] * designResult$surveyWeights)
    if(verbose) print(stat)
    stat
  }

  probFieldStatistic(field, nodes, designFun, resultFun, weights = weights)
}

#' @title Resample a Statistic
#' @description
#' Generate a set of bootstrap resampled statistics from a given \code{probField} and an
#' associated \code{probFieldStatistic}
#' @param x A \code{probFieldStatistic}
#' @param field The association \code{probField}
#' @param bootSamples The number of bootstrap resamples from the field's survey indices
#' @param doubleBootSamples The number of bootstrap resamples to take from each of the sets of
#' resampled indices above
#' @param ... Additional parameters (currently only included for the reference in the generic)
#' @return A list containing
#' \itemize{
#' \item{T1}{ The (numeric) statistic calculated from the field}
#' \item{T2}{ A vector of length \code{bootSamples} containing resampled statistics}
#' \item{T3}{ A list of length \code{bootSamples} elements, each element consisting of a vector
#' of length \code{doubleBootSamples} containing resampled statistics (see details)}
#' }
#' @details
#' The resample method relies on a \code{probField}'s ability to generate random survey indices from
#' the CDFs contained in its \code{probFieldNode}s (see \link{resampledSurveyIndices}).
#' Since the result function of a \code{probFieldStatistic} is designed to accept survey indices for
#' which to re-calculate the statistic, a distribution of potential statistics can be created.
#' The second bootstrap resampling does...
#' @export
resample <- function(x, ...) UseMethod("resample")

#' @describeIn resample
#' @export
resample.probFieldStatistic <- function(x, field, bootSamples, doubleBootSamples, ...){
  samples <- replicate(bootSamples, resampledSurveyIndices(field), simplify = FALSE)
  resamples <- lapply(bootSamples, function(s){
    replicate(doubleBootSamples, resampledSurveyIndices(field), simplify = FALSE)
  })

  T2s <- sapply(samples, x$result)
  T3s <- lapply(resamples, function(rs){
    sapply(rs, x$result)
  })

  list(T1 = x$result, T2 = T2s, T3 = T3s)
}

matchQuantile <- function(x, y, stat){
  rightXIndex <- min(which(x >= stat))

  xSurroundingStat <- data.frame(xi = c(rightXIndex, rightXIndex - 1),
                                 xv = x[c(rightXIndex, rightXIndex - 1)])
  xLine <- lm(xv ~ xi, data = xSurroundingStat)
  xInterpolated <- predict(xLine, list(xi = stat))
  xQuantile <- xInterpolated / length(x)

  quantile(y, xQuantile, names = FALSE, type = 4)
}
