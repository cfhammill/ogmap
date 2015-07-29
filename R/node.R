#' @title Create a node for a probability field
#' @description Create a node for a probability field, or to generate inference from a probability field.
#' A node contains a guide - a vector of predictor variables - and an optional random variable
#' @param guide An arbitrary length vector of predictor variables
#' (although for almost all purposes its length should match its associated field)
#' @param ranVar An option observation of a random variable of interest at the given guide
#' @param steps The step probabilities for ordered nodes in an associated field
#' @param ogive The cdf of step probabilities for ordered nodes in an associated field
#' @param redundancy The inferential redundancy of this node relative to the other nodes in
#' an associated field
#' @return A \code{probFieldNode} object containing all the parameter values
#' @details Users should strongly consider avoiding manually supplying steps,
#' ogives, and redundancy values. These values are typically calculated by
#' \link{calcCDF} and only make sense when in context of another set of nodes
#' @export
node <- function(guide, ranVar = NULL, steps = NULL, ogive = NULL, redundancy = 1){
  #Add weights as additional argument
  node <- structure(list(guide = guide, ranVar = ranVar,
                         steps = steps, ogive = ogive,
                         redundancy = redundancy),
                    class = c("probFieldNode", "list"))

}

nodeDistanceMatrix <- function(nodeSource, distance){
  nodes <- as.nodes(nodeSource)
  nNodes <- length(nodes)

  dists <- matrix(0, ncol = nNodes, nrow = nNodes)

  for(i in 1:(nNodes - 1)){
    for(j in (i + 1):nNodes){
      dists[i,j] <- dists[j,i] <- distance(nodes[[i]], nodes[[j]])
    }
  }

  dists
}

#' @title Add Information to a nodes from a field
#' @description
#' For a node or nodes within the domain of a probability field, calculate the CDF for the new
#' node relative to the survey nodes of the field, as well as the field's kernel and
#' distance functions.
#' @param nodes The new \code{probFieldNode} or \code{probFieldNodeList}
#' @param field The \code{probField} against which to evaluate the node's CDF
#' @export
evaluateNodes <- function(nodes, field) UseMethod("evaluateNodes")

#' @describeIn evaluateNodes
#' @export
evaluateNodes.probFieldNode <- function(nodes, field){
  calcCDF(nodes, field$survey, field$kernel, field$distance)
}

#' @describeIn evaluateNodes
#' @export
evaluateNodes.probFieldNodeList <- function(nodes, field){
  evaluatedNodes <- as.nodes(lapply(nodes, evaluateNodes, field))
  evaluatedNodes
}

#' @title getGuide
#' @description Get the guide variable(s) for a given node or collection thereof
#' @param x A \code{probFieldNode}, \code{probFieldNodeList}, or \code{probField}
#' @param returnVector Force the guide(s) to be returned as a vector
#' @param ... Additional parameters (currently only included for the reference in the generic)
#' @return For a \code{probFieldNode}, a vector containing the guide values for a node,
#' or if returnVector is set to false, potentially a one row data.frame depending on how
#' the node was constructed. For a collection of nodes returns a data.frame of guide variable values
#' @export
getGuide <- function(x, ...) UseMethod("getGuide")

#' @describeIn getGuide
#' @export
getGuide.probFieldNode <-  function(x, returnVector = TRUE, ...){
  if(! "probFieldNode" %in% class(x)) stop("Please supply a valid probFieldNode")
  if(returnVector) return(unlist(x$guide))
  return(x$guide)
}

#' @describeIn getGuide
#' @export
getGuide.probFieldNodeList <- function(x, ...){
  lapply(x, getGuide, returnVector = FALSE)
}

#' @describeIn getGuide
#' @export
getGuide.probField <- function(x, ...) getGuide(x$survey)

#' @title getVariable
#' @description Return the random variable observation of a \code{probFieldNode} or collection thereof
#' @param x A \code{probFieldNode}, \code{probFieldNodeList}, or \code{probField}
#' @param ... Additional parameters (currently only included for the reference in the generic)
#' @return For a node, either a number or NA depending on whether the node contains an
#' observation. For a nodeList, a vector of observations
#' @export
getVariable <- function(x) UseMethod("getVariable")

#' @describeIn getVariable
#' @export
getVariable.probFieldNode <- function(x){
  if(! "probFieldNode" %in% class(x)) stop("Please supply a valid probFieldNode")
  x$ranVar
}

#' @describeIn getVariable
#' @export
getVariable.probFieldNodeList <- function(x) sapply(x, getVariable)

#' @describeIn getVariable
#' @export
getVariable.probField <- function(x) getVariable(x$survey)

nodeList <- function(){

}

#' @title Convert data to nodes
#' @description Generic function to convert arbitrary data types into \code{probFieldNodeList} objects
#' @param x A data.frame or vector
#' @param response Which column of a data frame to be considered the random variable
#' @return A \code{probFieldNodeList}
#' @param ... Additional parameters (currently only included for the reference in the generic)
#' @details
#' If data frame is supplied each row is converted into a node with one column optionally representing
#' the response variable, if a numeric vector is supplied each element is considered a guide variable
#' for a new node. Lists are cast as class \code{probFieldNodeList}, \code{ProbField}s have their
#' survey nodes returned, and \code{probFieldNodeLists} are returned unmodified.
#' @export
as.nodes <- function(x, ...) UseMethod("as.nodes")

#' @describeIn as.nodes Convert a data.frame to a \code{probFieldNodeList}
#' @export
as.nodes.data.frame <- function(x, response = NULL, ...){
  nodeList <- structure(list(), class = c("probFieldNodeList", "list"))
  for(i in 1:nrow(x)){
    if(!is.null(response)){
      responseIndex <- which(names(x) == response)
      nodeList[[i]] <- node(guide = x[i,-responseIndex], ranVar = x[i, responseIndex])
    } else {
      nodeList[[i]] <- node(x[i,])
    }
  }

  return(nodeList)
}

#' @describeIn as.nodes
#' @export
as.nodes.numeric <- function(x, ...){
  structure(lapply(x, node),
            class = c("probFieldNodeList", "list"))
}

#' @describeIn as.nodes
#' @export
as.nodes.list <- function(x, ...){
  if(! "probFieldNode" %in% class(x[[1]])){
    stop("Cannot coerce list not containing probFieldNodes to type probFieldNodeList")
  }

  class(x) <- c("probFieldNodeList", "list")
  x
}

#' @describeIn as.nodes
#' @export
as.nodes.probFieldNodeList <- function(x, ...){
  x
}

#' @describeIn as.nodes
#' @export
as.nodes.probField <- function(x, ...){
  x$survey
}

#' @title getIndex
#' @description
#' Return the index identifying the position of a given node in a list of nodes
#' or a probability field
#' @param node The \code{probFieldNode} of interest
#' @param nodeListSource Either a \code{probFieldNodeList} or a \code{probField}
#' from which to retreave a nodeList of survey nodes.
#' @return A numeric index specifying the position of the given node in the list
#' @export
getIndex <- function(node, nodeListSource){
  nodeList <- as.nodes(nodeListSource)

  which(sapply(nodeList, identical, node, ignore.environment = TRUE))
}

#' @title Retrieve Step Probabilities
#' @description
#' Return the step probabilities for a \code{probFieldNode} or a group thereof
#' @param x A \code{probFieldNode}, \code{probFieldNodeList}, or \code{probField} from which to obtain step probabilities
#' @param asMatrix Logical whether the results should be returned as matrix or list when retrieving steps from a
#' collection of nodes
#' @param ... Additional parameters (currently only included for the reference in the generic)
#' @return A vector in the case x is a \code{probFieldNode} otherwise a matrix or list of such vectors at the users
#' discretion
#' @export
getSteps <- function(x, ...) UseMethod("getSteps")

#' @describeIn getSteps
#' @export
getSteps.probFieldNode <- function(x, ...) x$steps

#' @describeIn getSteps
#' @export
getSteps.probFieldNodeList <- function(x, asMatrix = TRUE, ...){
  steps <- lapply(x, getSteps)
  if(asMatrix) steps <- matrix(unlist(steps), ncol = length(steps))
  steps
}

#' @describeIn getSteps
#' @export
getSteps.probField <- function(x, asMatrix = TRUE, ...){
  nodeList <- x$survey
  getSteps(nodeList, asMatrix = asMatrix)
}

#' @title Retrieve Node CDFs
#' @description
#' Return the CDF for a \code{probFieldNode} or a group thereof
#' @param x A \code{probFieldNode}, \code{probFieldNodeList}, or \code{probField} from which to obtain the CDF(s)
#' @param asMatrix Logical whether the results should be returned as matrix or list when retrieving CDFs from a
#' collection of nodes
#' @param ... Additional parameters (currently only included for the reference in the generic)
#' @return A vector in the case x is a \code{probFieldNode} otherwise a matrix or list of such vectors at the users
#' discretion
#' @export
getCDF <- function(x, ...) UseMethod("getCDF")

#' @describeIn getCDF
#' @export
getCDF.probFieldNode <- function(x, ...) x$ogive

#' @describeIn getCDF
#' @export
getCDF.probFieldNodeList <- function(x, asMatrix = TRUE, ...){
  CDF <- lapply(x, getCDF)
  if(asMatrix) CDF <- matrix(unlist(CDF), ncol = length(CDF))
  CDF
}

#' @describeIn getCDF
#' @export
getCDF.probField <- function(x, asMatrix = TRUE, ...){
  nodeList <- x$survey
  getCDF(nodeList, asMatrix = asMatrix)
}

#' @title Retrieve Node Redundancy
#' @description
#' Return the redundancy of a survey node or a collection thereof
#' @param x A \code{probFieldNode}, \code{probFieldNodeList}, or \code{probField}
#' @return A numeric vector containing the redundancy of all nodes of interest
#' @export
getRedundancy <- function(x) UseMethod("getRedundancy")

#' @describeIn getRedundancy
#' @export
getRedundancy.probFieldNode <- function(x){x$redundancy}

#' @describeIn getRedundancy
#' @export
getRedundancy.probFieldNodeList <- function(x){sapply(x, getRedundancy)}

#' @describeIn getRedundancy
#' @export
getRedundancy.probField <- function(x){getRedundancy(x$survey)}

#' @title Find Relevant Indices
#' @description Find the indices corresponding to a subset of nodes from a group of nodes
#' (either a \code{probField} or \code{probFieldNodeList}) that has a given amount of
#' relevance to another group of nodes. The influence may be weighted by the importance
#' of nodes within the other group.
#' @param nodes The \code{probFieldNode} or \code{probFieldNodeList} to for which to find
#' their influential counterparts
#' @param weights A vector of weights representing the relative importance of influence on
#' each of the nodes given to the \code{nodes} parameter
#' @param surveySource Either a \code{probField} or \code{probFieldNodeList} from which to find
#' the indices of relevant member nodes
#' @param residue The percent of the right tail of the weighted relevance distribution to discard
#' @details surveySource should be a \code{probField} or \code{probFieldNodeList}, although
#' any object that can be coerced to a \code{probFieldNodeList} works, however doing so
#' is ill-advised
#' @return A vector of indices corresponding to the positions within \code{surveySource} that
#' are relevant to members of \code{nodes}
#' @export
findRelevantIndices <- function(nodes, surveySource, residue = 0.02, weights = NULL){

  survey <- as.nodes(surveySource)

  weightedRel <- weightedRelevance(nodes, surveySource, weights)

  cummulativeWeights <- cumsum(sort(weightedRel))
  aboveCutoff <- which(cummulativeWeights <= tail(cummulativeWeights, 1) * residue)

  if(length(aboveCutoff) == 0) return(numeric())

  breakPoint <- max(aboveCutoff)
  cutoffWeight <- cummulativeWeights[breakPoint + 1] - cummulativeWeights[breakPoint]
  relevantIndices <- which(weightedRel > cutoffWeight)

  relevantIndices
}

#' @title Weighted Relevance
#' @description Calculate the weighted relevance of each node in a \code{probField} or
#' \code{probFieldNodeList} to another group of nodes. The weighted relevance is the sum
#' of the step probability for each survey node summed accross the other nodes.
#' @inheritParams findRelevantIndices
#' @details surveySource should be a \code{probField} or \code{probFieldNodeList}, although
#' any object that can be coerced to a \code{probFieldNodeList} works, however doing so
#' is ill-advised
#' @return A numeric vector with length equal to the number of nodes in \code{surveySource}
#' with each element correspondign to the weighted sum of relevance to each of the members of
#' \code{nodes}
#' @export
weightedRelevance <- function(nodes, surveySource, weights = NULL){

  survey <- as.nodes(surveySource)
  if("probFieldNode" %in% class(nodes)) nodes <- as.nodes(list(nodes))

  if(is.null(weights)) weights <- rep(1, length(nodes))
  weights <- abs(weights)

  steps <- as.matrix(getSteps(nodes))
  weightedRelevance <- steps %*% weights

  as.numeric(weightedRelevance)
}

#' @title Plot a node
#' @description Plot the steps and CDF for a \code{probFieldNode}
#' @param x The \code{probFieldNode} of interest
#' @param ... Additional plot parameters to be added to both the steps and
#' CDF plots
#' @export
plot.probFieldNode <- function(x, ...){
  parO <- par(no.readonly = TRUE)
  par(mfrow = c(1,2))

  plot(getCDF(x), xlab = "", ylab = "CDF", type = "l", ...)
  plot(getSteps(x), xlab = "", ylab = "Step Probability", type = "l", ...)

  par(par)
  invisible(NULL)
}
