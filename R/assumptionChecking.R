########################################################################
# assumptionChecking.R
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

#' Jackknife probability of a node
#' @description
#' Calculate the jackknife probability of a node
#' @param x Either a \code{probFieldNode} or a numeric index to a node of interest
#' within a survey
#' @param field The probability field from which to obtain the survey nodes
#' @return A numeric probability
#' @export
jackknifeProb <- function(x, field){

  if("probFieldNode" %in% class(x)){
    node <- x
    index <- getIndex(node, field)
  }

  if(is.numeric(x)){
    node <- getNode(field, x)
    index <- x
  }

  node$ogive[index] / (1 - node$steps[index])
}

#' @title Assess potential uniformity
#' @description
#' Perform a Kolmogorov-Smirnov test to see if the jackknife probabilities
#' of survey nodes in a probability field could have been generated from a
#' uniform distribution
#' @param field The \code{probField} of interest
#' @return An \code{h-test} object. See return of \link{ks.test}
#' @export
plausibleUniformProbabilities <- function(field){
  nodeIndices <- seq(along = field$survey)
  probs <- sapply(nodeIndices, jackknifeProb, field)

  tryCatch(
    ks.test(abs(probs[-c(1, length(probs))] - 0.5), function(x){2 * x}),
    error = function(e){
      stop(paste0("ks.test error due on data = ", paste0(probs, collapse = " ")))
    },
    warning = function(w){
      stop(paste0("ks.test error due on data = ", paste0(abs(probs[-c(1, length(probs))] - 0.5), collapse = " ")))
    })
}




#' @title Jackknife pair probability difference
#' @description
#' Calculate the jackknife probability difference between a pair
#' of nodes
#' @param i Index of the first \code{probFieldNode} to be considered
#' @param j Index of the second \code{probFieldNode} to be considered
#' @param field the \code{probField} to be considered
#' @export
jackknifePairProbDiff <- function(i, j, field){

  if("probFieldNode" %in% class(i) && "probFieldNode" %in% class(j)){
    nodeI <- i
    nodeJ <- j
    indexI <- getIndex(i, field)
    indexJ <- getIndex(j, field)
  } else {
    nodeI <- getNode(field, i)
    nodeJ <- getNode(field, j)
    indexI <- i
    indexJ <- j
  }

  probI <- nodeI$ogive[i] / (1 - nodeI$steps[i] - nodeI$steps[j])
  probJ <- (nodeJ$ogive[j] - nodeJ$steps[i]) / (1 - nodeJ$steps[j] - nodeJ$steps[i])

  abs(probI - probJ)
}

#' @title Assess Paired Node Independence
#' @description
#' Assess the plausibility that pairs of nodes are independent by determining if the difference in probability
#' between the members of the pair are triangularly distributed
#' @param field The \code{probField} containing the nodes of interest
#' @param rawPairs A 2-column matrix or data.frame denoting the index of the members of the pairs with respect
#' to the original survey data
#' @return An \code{htest} object containing the results of a Kolmogorov-Smirnov test to determine if the paired
#' probability differences are triangularly distributed. See \link{ks.test} for more details
#' @export
plausibleIndependentPairs <- function(field, rawPairs){
  pairIndices <- rawPairs

  probDiffs <- apply(pairIndices, 1, function(row){
    jackknifePairProbDiff(row[1], row[2], field)
  })

  tryCatch(ks.test(probDiffs, function(x){x * (2-x)}),
           error = function(e){
             stop(paste0("ks.test error due on data = ", paste0(probDiffs, collapse = " ")))
             },
           warning = function(w){
             stop(stop(paste0("ks.test error due on data = ", paste0(probDiffs, collapse = " "))))
             })
}

#' @title Assess field assumptions
#' @description
#' Assess the two major assumptions made for field construction - pair independence and uniformity
#' @param field The \code{probField} of interest
#' @param rawPairs The pairs of interest to be passed to \link{plausibleIndependentPairs}
#' @return a numeric value equal to the product of the two p-values (see \link{plausibleUniformProbabilities}
#'  and \link{plausibleIndependentPairs}) divided by 1 minus the log that product.
#' @export
plausibleUniformIndependent <- function(field, rawPairs){
  pI <- plausibleIndependentPairs(field, rawPairs)$p.value
  pU <- plausibleUniformProbabilities(field)$p.value
  product <- pI * pU
  product * (1 - log(product))
}

#' @title Evaluate fit of a field
#' @description
#' Evaluate the fit of a field given a set of new bandwidths.
#' @param field The \code{probField} of interest
#' @param bandwidths A vector or list containing the new kernel and distance parameters to asses the cost at.
#' See \link{setBandwidths} and \link{getBandwidths} for more details
#' @param rawPairs Pairs of interest to be passed down to \link{plausibleIndependentPairs}
#' @param approximate Logical whether or not to approximate the cost
#' @param verbose Logical whether or not to print the bandwidth and the cost to the console
#' @return A numeric value equal to negative log the value from \link{plausibleUniformIndependent} or
#' and error code of 1000. Or in if using the approximation, return two times the product of the
#' KS statistics.
#' @export
costFunction <- function(field, bandwidths, rawPairs, approximate = FALSE, verbose = FALSE){
  updatedField <- setBandwidths(field, bandwidths)

  if(approximate){
    cost <- 2 * plausibleUniformProbabilities(field)$statistic +
      plausibleIndependentPairs(updatedField, rawPairs)$statistic

    cost <- as.numeric(cost) #Removes vector naming

  } else {
    pUI <- plausibleUniformIndependent(updatedField, rawPairs)
    if(!is.finite(pUI) || pUI <= 0) return(1000)
    cost <- -log(pUI)

  }


  if(verbose){print(paste0("Bandwidths ",
                           paste0(bandwidths, collapse = ", "),
                           " cost = ", cost))
              return(invisible(cost))
  }

  return(cost)
}

#' @title Find Close Pairs of Nodes
#' @description
#' Find a set of node pairing that describe a reasonable short path through the data.
#' Uses the nearest-neighbour solution to the travelling salesman problem using
#' a variety of start positions to find the shortest such path
#' @param nodeSource An object that can be coerced to a \code{probFieldNodeList} (through
#' \link{as.nodes})
#' @param distance The distance function used to evaluate inter-node distances
#' @param nStarts How many randomly selected nodes to try, defaults to all possible nodes,
#' alternatively a vector of indices to the starting nodes.
#' @return A 2-column matrix with the indices to each pair in the best path
#' @export
findNearbyPairs <- function(nodeSource, distance, nStarts = length(nodeSource)){
  nodes <- as.nodes(nodeSource)
  nNodes <- length(nodes)

  dists <- matrix(Inf, ncol = nNodes, nrow = nNodes)

  for(i in 1:(nNodes - 1)){
    for(j in (i + 1):nNodes){
      dists[i,j] <- dists[j,i] <- distance(nodes[[i]], nodes[[j]])
    }
  }

  saResults <- simulatedAnnealing(dists)

  pairs <- saResults$pairs
  pairs <- as.data.frame(pairs)
  pairs
}
