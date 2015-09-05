########################################################################
# Ogmap.R
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
