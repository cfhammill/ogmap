########################################################################
# probField.R
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

#' @title Generate a probability field
#' @description Generate a probability field from a data.frame of survey observations
#' @param x Either a \link{formula} or a data.frame of the form described below
#' @param data A data frame containing one column of an observed random variable
#' and an arbitrary number of predictor variable
#' @param kernel A kernel function that generates a weighting given a distance
#' see \link{kernel} for more information on \code{probFieldKernel}s
#' @param distance A function to be applied to two \code{probFieldNode} objects
#' (produced by \link{node}) to generate the distance between them
#' @param response The name of the column containing the random variable, defaults to "y"
#' if specified name is not provided will use the first column with a warning
#' @param redundancyCorrected Logical Whether or not to apply a redundancy correction
#' to the survey nodes, see \link{calcCDF}
#' @param calcCDF Logical whether or not to calculate the step probabilities, CDFs, and redundancy
#' for each node and attach them to the node objects
#' @param breakTies Logical whether or not to add a small number to each random variable
#' to break ties, see \link{breakTies}
#' @param ... Additional parameters (currently only included for the reference in the generic)
#' @return An object of class \code{probField} containing
#' \itemize{
#' \item{\code{survey}}{ An object of class \code{probFieldNodeList} containing the survey nodes
#' in order of increasing values of their random variable}
#' \item{\code{kernel}}{ The kernel function used to construct the field}
#' \item{\code{distance}}{ The distance function used to construct the field}
#' \item{\code{formula}}{ The formula representing the dependencies amongst variables in the field,
#' namely, response as a function of the guide variables}
#' \item{\code{sortOrder}}{ A sorting index to order the original data}
#' \item{\code{originalOrder}}{ A sorting index to unsort the modified data}
#' \item{\code{calcCDF}}{ Whether CDFs, steps, and redundancy were calculated and appended to nodes}
#' \item{\code{redundancyCorrected}}{ Whether the CDFs and steps were corrected for redundancy}
#' \item{\code{breakTies}}{ Whether ties were broken}
#' }
#' @details
#' This is the gatekeeper function of the ogmap package, taking a data.frame containing
#' random variable observations and one or many predictor variables. First the frame is sorted
#' and converted to a series of nodes. If requested the function will calculate the CDFs and step
#' probabilities for all the other nodes using kernel weightings and a distance function supplied
#' by the user. The field represents a base model from which inference regarding the
#' distribution of potential outcomes given predictor variables, instead of just the expectation
#' of the distribution
#' @import parallel
#' @export
probField <- function(x, ...) UseMethod("probField")

#' @describeIn probField
#' @export
probField.formula <- function(x, data, kernel, distance,
                              redundancyCorrected = TRUE, calcCDF = TRUE,
                              breakTies = FALSE, ...){

  originalData <- data
  modifiedData <- model.frame(x, data)
  response = names(modifiedData)[1]

  #Hand off work to createField
  createField(originalData = originalData, modifiedData = modifiedData, kernel = kernel,
              distance = distance, response = response, redundancyCorrected = redundancyCorrected,
              calcCDF = calcCDF, breakTies = breakTies, formula = x)

}

#' @describeIn probField
#' @export
probField.data.frame <- function(x, kernel, distance, response = "y",
                      redundancyCorrected = TRUE, calcCDF = TRUE,
                      breakTies = FALSE, ...){
  originalData <- x
  modifiedData <- x

  if(is.numeric(response) && response <= nrow(originalData)){
    response <- names(originalData)[response]
  } else if(!response %in% names(originalData)) {
    response <- names(originalData)[1]
    warning(paste("No response variable present, using", response, "instead"))
  } else if(!is.character(response)){
    stop("Invalid response index")
  }

  #reorder the data frame such that the response variable is column one, to make the conversion to formula
  #work
  modifiedData <- modifiedData[ ,c(response, names(modifiedData)[! names(modifiedData) %in% response])]
  modelFormula <- formula(modifiedData)

  #Hand off work to createField
  createField(originalData = originalData, modifiedData = modifiedData, kernel = kernel,
              distance = distance, response = response, redundancyCorrected = redundancyCorrected,
              calcCDF = calcCDF, breakTies = breakTies, formula = modelFormula)
}

createField <- function(originalData, modifiedData, kernel, distance, response,
                        redundancyCorrected = TRUE, calcCDF = TRUE,
                        breakTies = FALSE, formula){

  #Check if the user wants to break ties, if so break them
  if(breakTies){
    modifiedData[,response] <-
      breakTies(modifiedData[,response])
  }

  sortedOrder <- order(modifiedData[,response])
  originalOrder <- order(sortedOrder)

  modifiedData <- modifiedData[sortedOrder,]
  survey <- as.nodes(modifiedData, response)

  if(calcCDF) survey <- calcCDF(survey, kernel, distance, redundancyCorrected)

  field <- structure(list(survey = survey, kernel = kernel, distance = distance, formula = formula,
                          originalData = originalData, modifiedData = modifiedData,
                          sortOrder = sortedOrder, originalOrder = originalOrder,
                          calcCDF = calcCDF, redundancyCorrected = redundancyCorrected,
                          breakTies = breakTies),
                     class = "probField")

  field
}

#' @title Subset a Field
#' @description Create a subset of the field
#' @param field The \code{probField} to subset
#' @param indices Numeric or logical index to the survey nodes of the probability field
#' @details Note that survey nodes are in order of increasing random variable, and not
#' the original order of the data, indices need to reflect this order
#' @return A \code{probField}
#' @export
subsetField <- function(field, indices){
  newField <- field

  newField$survey <- calcCDF(as.nodes(field$survey[indices]), field$kernel,
                       field$distance, field$redundancyCorrected)

  newField$originalData <- field$originalData[field$originalOrder[indices],]
  newField$modifiedData <- field$modifiedData[indices,]
  newField$sortOrder <- field$sortOrder[field$originalOrder[indices]]
  newField$originalOrder <- field$originalOrder[indices]

  newField
}



#' @title Tie breaking
#' @description The functions takes a vector of observations, generates a sample from 1 to the number of
#' observations without replacement, and multiplies that by the smallest positive double value
#' \code{.Machine$double.xmin} and adds that to the original observations. The procedure should
#' work in almost all cases except when values are measured with precision of 10e-308, in
#' which case ties should be unlikely at the outset.
#' @param x A vector of observations
#' @return A vector of observations with ties broken
#' @export
breakTies <- function(x){
  sample(1:length(x)) * .Machine$double.xmin + x
}

intraSurveyRelevance <- function(survey, kernel, distance){
  relevance <- matrix(1, nrow = length(survey), ncol = length(survey))
  for(i in 1:(length(survey) - 1)){
    for(j in (i+1):length(survey)){
      relevance[i, j] <- relevance[j, i] <-
        nodeRelevance(survey[[i]], survey[[j]], kernel, distance)
    }
  }
  relevance
}

intraSurveyRelevanceV <- function(survey, kernel, distance){
  relevance <- matrix(1, nrow = length(survey), ncol = length(survey))

  lapply(1:(length(survey) - 1), function(i){
    lapply((i+1):length(survey), function(j){
      rel <- nodeRelevance(survey[[i]], survey[[j]], kernel, distance)
      relevance[i,j] <<- rel
      relevance[j,i] <<- rel
      invisible(NULL)
    })
  })

  relevance
}


intraSurveyRelevanceP <- function(survey, kernel, distance){
  relevance <- matrix(1, nrow = length(survey), ncol = length(survey))

  mclapply(1:(length(survey) - 1), function(i){
    mclapply((i+1):length(survey), function(j){
      rel <- nodeRelevance(survey[[i]], survey[[j]], kernel, distance)
      relevance[i,j] <<- rel
      relevance[j,i] <<- rel
      invisible(NULL)
    })
  })

  relevance
}


intraSurveyRelevanceP2 <- function(survey, kernel, distance){
  nnodes <- length(survey)

  samplePoints <- expand.grid(1:nnodes,
                              1:nnodes)
  indicesFlat <- matrix(1:(nnodes^2), nrow = nnodes, ncol = nnodes)
  uniqueSamplePoints <- indicesFlat[upper.tri(indicesFlat)]

  indexI <- samplePoints[uniqueSamplePoints,1]
  indexJ <- samplePoints[uniqueSamplePoints,2]

  if(is.null(options('ogmapMaxCores')[[1]])){
    ncores <- detectCores() - 1
  } else {
    ncores <- min(options('ogmapMaxCores')[[1]], detectCores())
  }

  rel <- t(mcmapply(function(i, j){
    c(i, j, nodeRelevance(survey[[i]], survey[[j]], kernel, distance))
  }, i = indexI, j = indexJ, mc.cores = ncores))

  gc()
  relevance <- matrix(1, nrow = nnodes, ncol = nnodes)

  apply(rel, 1, function(row){
    relevance[row[1], row[2]] <<- row[3]
    relevance[row[2], row[1]] <<- row[3]
  })

  relevance
}

#' @title Inter Node Relevance
#' @description Calculation to determine the inferential relevance of one node to another.
#' Simply evaluates the kernel function on the distance between two nodes given
#' a distance function.
#' @param node1 a \code{probFieldNode} object
#' @param node2 another \code{probFieldNode}
#' @param kernel the kernel function of interest
#' @param distance the distance function of interest
#' @export
nodeRelevance <- function(node1, node2, kernel, distance){
  kernel(distance(node1, node2))
}

#' @title Generate the CDF, Step probabilities, and Redundancy
#' @description Generates the CDF, step probabilities, and redundancy for nodes either for a
#' list of nodes, or for one node relative to another list of nodes
#' @param x Either A \code{probFieldNode} or a \code{probFieldNodeList}
#' @param nodeList A \code{probFieldNodeList} with nodes containing precalculated redundancy
#' @param kernel The kernel function of interest
#' @param distance The distance function of interest
#' @param redundancyCorrected Whether or not to adjust for inter node redundancy
#' @param ... Additional parameters (currently only included for the reference in the generic)
#' @return An object with the same type as x with step probabilities, CDF(s), and redundancy
#' attached
#' @export
calcCDF <- function(x, ...){ UseMethod("calcCDF")}

#' @describeIn calcCDF Generate CDFs, step probabilities, and redundancy for a
#' \code{probFieldNodeList}
#' @export
calcCDF.probFieldNodeList <-
  function(x, kernel, distance, redundancyCorrected = TRUE, ...){
    survey <- x
    weights <- calcWeights(survey, kernel, distance)
    redundancy <- 1
    if(redundancyCorrected){
      redundancy <- colSums(weights)
      weights <- weights / redundancy
    }

    ogives <- apply(weights, 2, calcOgive)
    steps <- apply(weights, 2, calcSteps)

    for(i in 1:length(survey)){
      survey[[i]]$redundancy <- redundancy[i]
      survey[[i]]$steps <- steps[,i]
      survey[[i]]$ogive <- ogives[,i]
    }

    return(survey)
  }

#' @describeIn calcCDF
#' @export
calcCDF.list <- function(x, ...){
  class(x) <- c("probFieldNodeList", "list")
  calcCDF(x, ...)
}

#' @describeIn calcCDF Calculate the CDF, step probabilities, and redundancy for a \code{probFieldNode}
#' given a \code{probFieldNodeList}
#' @export
calcCDF.probFieldNode <-
  function(x, nodeList, kernel, distance, redundancyCorrected = TRUE, ...){
    weights <- calcWeights(x, nodeList, kernel, distance)
    redundancy <- 1
    if(redundancyCorrected) redundancy <- getRedundancy(nodeList)
    x$ogive <- calcOgive(weights/redundancy)
    x$steps <- calcSteps(weights/redundancy)

    return(x)
}

#' @title Retrieve a node
#' @description
#' Function to retrieve a specific node in a field from either it's guide, value, or index
#' @param field The probability field of interest
#' @param lookup The value be found
#' @param type A choice of index, value, or guide, (defaulting to index) for how to find
#' the node of interest
#' @param tolerance If matching by guide or value, what is the maximum difference between
#' the lookup and the node value to be considered a match
#' @param allowMultiples when matching by value or guide should multiple values potentially
#' be returned, if not the first match is returned
#' @return A \code{probFieldNode} or NULL if no match is found when allowMultiples is set
#' to FALSE, else, a \code{probFieldNodeList} containing one or more matches
#' @export
getNode <- function(field, lookup, type = c("index", "value", "guide"), tolerance = 10^-5,
                    allowMultiples = FALSE){
  if(! "probField" %in% class(field)) stop("Please supply a valid probField object")
  type <- match.arg(type)
  index <- switch(type,
         "index" = lookup,

         "value" = which(sapply(getVariable(field$survey),
                                              fuzzyEqual, lookup, tolerance = tolerance,
                                              global = TRUE)),

         "guide" = which(sapply(getGuide(field$survey),
                                              fuzzyEqual, lookup, tolerance = tolerance,
                                              global = TRUE))
         )
  if(length(index) == 0) return(NULL)
  if(allowMultiples) return(structure(field$survey[index], class = c("probFieldNodeList", "list")))
  return(field$survey[[index[1]]])
}

#' @title Retrieve bandwidths
#' @description
#' Return all additional parameters present in the kernel and distance components of a
#' probability field
#' @param field The \code{probField} of interest
#' @return A named character vector of additional parameters with an additional \code{from}
#' attribute specifying whether the parameters belong to the kernel or distance function
#' @export
getBandwidths <- function(field){
  shape <- unlist(attr(field$kernel, "kernelParameters"))
  scale <- unlist(attr(field$distance, "distanceParameters"))
  bandwidths <- c(shape, scale)
  attr(bandwidths, "from") <- rep(c("kernel", "distance"), c(length(shape), length(scale)))
  bandwidths
}

#' @title Set new bandwidths
#' @description
#' Set new bandwidths for a probability field and recalculate CDFs, steps, and redundancy
#' for survey nodes
#' @param field The \code{probField} of interest
#' @param bandwidths A vector, named vector, or list specifying the new bandwidths
#' @export
setBandwidths <- function(field, bandwidths){
  oldBandwidths <- getBandwidths(field)
  newField <- field

  paramNames <- names(oldBandwidths)
  kernelParams <- attr(oldBandwidths, "from") == "kernel"

  kernelArgs <- as.list(bandwidths[kernelParams])
  names(kernelArgs) <- paramNames[kernelParams]

  distanceArgs <- as.list(bandwidths[!kernelParams])
  names(distanceArgs) <- paramNames[!kernelParams]

  guides <- attr(newField$distance, "guides")

  list2env(kernelArgs, as.environment(-1))
  list2env(distanceArgs, as.environment(-1))

  if(!is.null(guides)) list2env(guides, as.environment(-1))

  newField$kernel <- eval(attr(newField$kernel, "call"))
  newField$distance <- eval(attr(newField$distance, "call"))

  newField$survey <- calcCDF(newField$survey, newField$kernel, newField$distance)

  newField
}

#' @title Assessing Variable Probability
#' @description
#' Function to calculate the probability of a random variable being observed at a
#' given node
#' @param field The \code{probField} to be assessed
#' @param randomVariable The random variable value of interest
#' @param node The \code{probFieldNode} of interest
#' @return The value of the CDF of the node of interest evaluated at the point where
#' randomVariable exceeds the observed value in the field
#' @description
#' This function assess the probability of an observed value, its inverse is the process
#' of assessing the value at a given probability see \link{variableAtProbability}
#' @export
probabilityOfVariable <- function(field, randomVariable, node){
  observedVar <- getVariable(field$survey)
  if(randomVariable <= observedVar[1]) stop("Variable out of range", call. = FALSE)

  index <- min(which(observedVar > randomVariable))
  node$ogive[index]
}

#' @title Finding probable value
#' @description
#' Assess the maximum value expected at a given probability from a field at a specific
#' node
#' @param field The \code{probField} of interest
#' @param probability The probability the value must not exceed
#' @param node The \code{probFieldNode} who's CDF is to be examined
#' @return A numeric value of the random variable
#' @description
#' This function is concerned with calculating the value indicated by a probability
#' for its inverse see \link(probabilityOfVariable)
#' @export
variableAtProbability <- function(field, probability, node){
  og <- node$ogive

  if(probability < og[2]){
    stop("Probability not in range", call. = FALSE)
  }

  index <- max(which(og[-1] < probability))
  getVariable(field$survey)[index + 1]
}

#' @title Node Summary Statistics
#' @description
#' The expected value, variance, or diversity of a given node
#' @param node The \code{probFieldNode} of interest
#' @param field The \code{probField} of interest
#' @param confidence The level of confidence for the interval
#' @param resamples The number of bootstrap resamples used to calculate confidence
#' limits
#' @return The numeric value of the statistic for the node for expectation, variance
#' and diversity, and a length-3 named vector with lower-limit, expectation, upper-limit
#' for confidence intervals
#' @export
nodeExpectation <- function(node, field){
  sum(getVariable(field$survey) * node$steps)
}

#' @describeIn nodeExpectation
#' @export
nodeVariance <- function(node, field){
  m <- nodeExpectation(node, field)
  msq <- sum(node$steps * getVariable(field$survey)^2)
  msq - m^2
}

#' @describeIn nodeExpectation
#' @export
nodeDiversity <- function(node){
1 / sum(node$steps^2)
}

#' @describeIn nodeExpectation
#' @export
nodeConfidenceIntervals <- function(node, field, confidence = .95, resamples = 1000){
  expectation <- nodeExpectation(node, field)

  resampledExpectation <- replicate(resamples, {
    resampledSurvey <- field$survey[resampledSurveyIndices(field)]
    class(resampledSurvey) <- c("probFieldNodeList", "list")

    resampledField <- field
    resampledField$survey <- resampledSurvey

    nodeExpectation(node, resampledField)
  }, simplify = TRUE)

  resampledExpectation <- sort(resampledExpectation)
  lowerIndex <- floor((1 - confidence)/2 * resamples)
  upperIndex <- resamples - lowerIndex

  confidenceIntervalVector <- c(resampledExpectation[lowerIndex],
                                expectation,
                                resampledExpectation[upperIndex])

  names(confidenceIntervalVector) <- c(paste0(lowerIndex / resamples * 100,"%"),
                                       "expectation",
                                       paste0(upperIndex / resamples * 100, "%"))
  confidenceIntervalVector
}

#' @title Random sample of survey indices
#' @description
#' Generate a random sample of survey indices for randomization and monte carlo methods
#' @param x The \code{probField} or \code{probFieldNodeList} of interest
#' @return A vector of the same length as the \code{probFieldNodeList} (or the survey within the \code{probField})
#' with random indices generated from the probability distribution at each node
#' @export
resampledSurveyIndices <- function(x) UseMethod("resampledSurveyIndices")

#' @describeIn resampledSurveyIndices
#' @export
resampledSurveyIndices.probField <- function(x) resampledSurveyIndices(x$survey)

#' @describeIn resampledSurveyIndices
#' @export
resampledSurveyIndices.probFieldNodeList <- function(x){
  probs <- runif(length(x))
  indices <- mapply(function(node, p){
    max( max(which(node$ogive < p)) - 1, 1)
  }, x, probs)

  indices
}

#' @title Get Survey Nodes
#' @description Return a \code{probFieldNodeList} either identical to its argument
#' if a \code{probFieldNodeList} is passed, or the contained list of survey nodes is a
#' \code{probField} is passed
#' @param surveySource The source of survey nodes
#' @return A \code{probFieldNodeList} or NULL if the argument is neither a \code{probField} or
#' \code{probFieldNodeList}
#' @export
getSurvey <- function(surveySource){
  survey <- NULL

  if("probField" %in% class(surveySource)){
    survey <- surveySource$survey
  } else if("probFieldNodeList" %in% class(surveySource)){
    survey <- surveySource
  }

  survey
}

#############################################
#Ancillary Functions
#############################################

calcWeights <- function(x, ...) UseMethod("calcWeights")

calcWeights.probFieldNodeList <-
  function(x, kernel, distance, ...){
    nodeList <- x

    if(.Platform$OS.type != "unix"){
      relevance <- intraSurveyRelevance(nodeList, kernel, distance)
    } else {
      relevance <- intraSurveyRelevanceP2(nodeList, kernel, distance)
    }

    relevance
  }

calcWeights.probFieldNode <-
  function(x, nodeList, kernel, distance, ...){
    relevance <- sapply(nodeList, nodeRelevance, node1 = x, kernel = kernel, distance = distance)

    relevance
  }

calcOgive <- function(weights){
  cWeights <- cumsum(c(0, weights))
  norm <- tail(cWeights, 1)
  if(norm == 0) norm <- 1

  cWeights / norm
}

calcSteps <- function(weights){
  cWeights <- cumsum(c(0, weights))
  norm <- tail(cWeights, 1)
  if(norm == 0) norm <- 1

  weights / norm
}


