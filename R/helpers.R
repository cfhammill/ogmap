fuzzyEqual <- function(x, y, tolerance, global = FALSE){
  fe <- abs(x - y) < tolerance
  if(global) fe <- all(fe)
  fe
}

#' @title Set maximum number of cores
#' @description Set maximum number of cores to use when fitting \code{probfield}s
#' @param n The maximum number of cores
#' @return invisible null
#' @export
setMaxCores <- function(n){
  options("ogmapMaxCores" = n)
  invisible(NULL)
}

findNNPath <- function(startIndex, distanceMatrix){
  dm <- distanceMatrix
  currentNode <- startIndex
  path <- currentNode

  #Find the path
  for(i in 2:ncol(dm)){
    nextNode <- which(min(dm[-path,currentNode]) == dm[,currentNode])

    if(length(nextNode) > 1){
      nextNode <- nextNode[!nextNode %in% path]
      nextNode <- nextNode[1]
    }

    path <- as.numeric(c(path, nextNode))
    currentNode <- nextNode
  }

  #Determine path length
  n1 <- path[-length(path)]
  n2 <- path[-1]

  pathLength <- sum(mapply(
    function(indexI, indexJ){
      dm[indexI, indexJ]
    },
    n1[c(TRUE,FALSE)], n2[c(TRUE,FALSE)]))

  #Return node pairs, path, and pathlength of a solution
  pairs <- data.frame(node1 = n1[c(TRUE,FALSE)], node2 = n2[c(TRUE,FALSE)])
  list(pl = pathLength, path = path, pairs = pairs)
}

nnTravellingSalesman <- function(distanceMatrix, starts = ncol(distanceMatrix), returnBest = TRUE){
  if(length(starts) == 1) starts <- sample(1:ncol(distanceMatrix), starts)

  nnResults <- lapply(starts, findNNPath, distanceMatrix)

  if(returnBest){
    pathLengths <- sapply(nnResults, "[[", "pl")
    bestSolution <- which(pathLengths == min(pathLengths))[1]
    return(nnResults[[bestSolution]])
  }

  nnResults
}

simulatedAnnealing <- function(distanceMatrix, extraPairs = 0, temperature = 1){
  # This version of the algorithm takes the precomputed distance matrix from another function

  nNodes <- nrow(distanceMatrix)
  nOver <- 100 * nNodes
  nLimit <- 10 * nNodes
  nIndices <- nNodes  + nNodes %% 2 + extraPairs * 2

  #If there are an even number of nodes, add dummy rows and columns to the distance matrix
  if(nIndices > nNodes){
    nExtraNodes <- nIndices - nNodes
    distanceMatrix <- cbind(distanceMatrix,
                            matrix(rep(0, nrow(distanceMatrix) * nExtraNodes), ncol = nExtraNodes))
    distanceMatrix <- rbind(distanceMatrix,
                            matrix(rep(0, ncol(distanceMatrix) * nExtraNodes), nrow = nExtraNodes))
  }

  currentTotalDistance <- 0

  #Create the initial pairs
  oldPairs <- rbind(seq(1, nIndices, 2), seq(2, nIndices, 2))

  #Calculate original distance
  oldDistance <- sum(apply(oldPairs, 2, function(col) distanceMatrix[col[1], col[2]]))

  # Tally of changes and shuffles
  totalChanges <- 0
  totalSteps <- 0

  # Perform 100 temperature steps
  for(tempStep in 1:100){

    # Keep a tally of changes in each temp step
    nChanges <- 0

    # Perform nOver shuffles
    for(shuffles in 1:nOver){

      #Increment the shuffle counter
      totalSteps <- totalSteps + 1

      #Sample two random indices to swap, discard them if they are one apart (just swapping pair order)
      swapIndices <- sample(1:nIndices, 2)
      if(abs(swapIndices[1] - swapIndices[2]) == 1) next

      #Create an identical matrix of new pairs
      newPairs <- oldPairs

      #Exchange indices 1 and 2 in the new pairs
      newPairs[swapIndices[1]] <- oldPairs[swapIndices[2]]
      newPairs[swapIndices[2]] <- oldPairs[swapIndices[1]]

      #R is column major, so the cieling of the index divided by two yeilds the column
      #This code extracts the 4 pairs, the original 2 and the new 2
      op1 <- oldPairs[, ceiling(swapIndices[1] / 2)]
      op2 <- oldPairs[, ceiling(swapIndices[2] / 2)]
      np1 <- newPairs[, ceiling(swapIndices[1] / 2)]
      np2 <- newPairs[, ceiling(swapIndices[2] / 2)]

      #This calculates the distance between each pair
      od1 <- distanceMatrix[op1[1], op1[2]]
      od2 <- distanceMatrix[op2[1], op2[2]]
      nd1 <- distanceMatrix[np1[1], np1[2]]
      nd2 <- distanceMatrix[np2[1], np2[2]]

      #The distance change is the sum of new distances minus the old ones
      distanceChange <- nd1 + nd2 - od1 - od2

      #If the distance change is less than 0, or a random number on [0,1) is less than
      #e to the negative distance change divided by the temperature
      if(distanceChange < 0 || runif(1) < exp(-distanceChange/temperature)){
        oldPairs <- newPairs
        nChanges <- nChanges + 1
      }

      #If there have been more than nLimit swapped pairs move to the next temp step
      if(nChanges > nLimit) break
    }

    #Add the number of changes to the total changes tally
    totalChanges <- totalChanges + nChanges

    #If you didn't make any changes this temp step leave the temperature loop
    if(nChanges == 0) break

    #Decrease the temperature
    temperature <- temperature * .9
  }

  #Calculate the final distance between pairs
  finalDistance <- sum(apply(oldPairs, 2, function(col) distanceMatrix[col[1], col[2]]))

  #Remove any pairs with an imaginary node member
  oldPairs <- oldPairs[, oldPairs[1,] <= nNodes & oldPairs[2,] <= nNodes]

  #Sort the pairs such that they are organized low to high
  finalPairs <- t(apply(oldPairs, 2, sort))

  #Return some information about the solution and algorithm
  list(pairs = finalPairs,
       length = finalDistance,
       temperatureSteps = tempStep,
       numberChanges = totalChanges,
       totalSteps = totalSteps)
}
