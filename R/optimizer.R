#' @title Bandwidth Optimizer
#' @description determine the optimum kernel and distance parameters for a \code{probField} given paired observations
#' using either Nelder-Mead optimization via \link{optim} or covariance matrix adapting evolutionary strategy
#' optimization via \link{cma_es}
#' @param field The \code{probField} for which to optimize kernel and distance function parameters
#' @param pairs A 2 column matrix or \code{data.frame} containing integer indices to pairs of survey nodes to be
#' used in the optimization function
#' @param method  The optimization algorithm to use, currently either optim (Nelder-Mead; default) or cmaes (Covariance
#' matrix adapting evolutionary strategy)
#' @param ... Additional parameters to be passed to the optimization algorithm
#' @param returnField Whether or not to return the new field with optimized bandwidths or the bandwidth values
#' @param verbose Whether or not to print the optimized bandwidths
#' @return If returnField is TRUE returns a \code{probField} with the optimized bandwidths, if FALSE the results of
#' the optimizer function (see \link{optim} and \link{cma_es} for details)
#' @import cmaes
#' @export
optimizeBandwidths <- function(field, pairs,
                               method = c("optim", "cmaes"), ...,
                               verbose = FALSE,
                               returnField = TRUE){
  method <- match.arg(method)
  optimFun <- switch(method,
                     #optim = optim
                     optim = function(...) optim(..., method = "L-BFGS-B", lower = 0),
                     cmaes = cma_es)

  evaledPairs <- pairs
  startingBandwidths <- getBandwidths(field)

  res <- optimFun(startingBandwidths, costFunction, field = field, rawPairs = evaledPairs, ...)
  if(verbose) print(res$par)
  if(returnField) return(setBandwidths(field, res$par))
  return(res)
}
