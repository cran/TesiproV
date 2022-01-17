#' TesiproV: A package for the calculation of reliability and failure probability in civil engineering
#'
#' The Package provides three main types of objects:
#' \enumerate{
#' \item Objects for modeling base variables
#' \item Objects for modeling limit state functions and systems of them
#' \item Objects for modeling solving algorithms
#' }
#'
#'
#' By creating and combining those objects, one is able to model quite complex problems in terms of structural reliablity calculation.
#' For normally distributed variables there might be an workflow to calculate correlated problems (but no systems then).
#' There is also implemented a new distribution (logStudentT, often used for conrete compression strength) to show how one can implement
#' your very own or maybe combined multi modal distribution and use it with TesiproV.
#'
#'
#' @section Objects for base variables:
#' \code{PROB_BASEVAR}, \code{PROB_DETVAR}, \code{PARAM_BASEVAR}, \code{PARAM_DETVAR}
#'
#' @section Limit state functions:
#' \code{SYS_LSF}, \code{PROB_SYS}, \code{PARAM_SYS}
#'
#' @section Solving algorithms:
#' \code{PROB_MACHINE}
#'
#' @author (C) 2021 -  K. Nille-Hauf, T. Feiri, M. Ricker - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#'
#' @docType package
#' @name TesiproV

NULL
