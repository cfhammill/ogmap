#' @title Probability distribution estimation dependent on covariates
#'
#' @description
#' The ogmap package is designed to overcome the limitation of using only the expected value of a distribution for
#' inference. The package uses non-parametric weighting of survey observations to allow the entire distribution of
#' potential values given predictor variables to be examined.
#'
#'This package was designed with funding from the Canadian government's Department of Fisheries and Oceans. This
#'is a port of a python package designed, written, and conceived by Geoff Evans, with the port created by Chris Hammill
#'
#' @author Chris Hammill
#' @author Geoff Evans
#' @docType package
#' @name ogmap
NULL

#' @title Mock Stock Recruitment Data
#' @description A dataset containing the adult stock and juvenile recruitment for
#' \emph{Imaginarifish examplicus}
#' @format A data frame with 49 rows of 2 variables:
#' \describe{
#'   \item{stock}{Adult stock, in 10,000s}
#'   \item{recruitment}{Juvenile's recruited, in 10,000s}
#' }
#' @docType data
#' @name stockRecruitment
NULL
