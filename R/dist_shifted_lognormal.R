#' @name shifted_lnorm
#' @title Shifted lognormal distribution
#' @description
#' Density, distribution function, quantile function and random generation of the
#' shifted lognormal distribution. The random variable is \code{shift + Y} with
#' \eqn{\log(Y) \sim N(meanlog, sdlog^2)}. The parameter names and the order
#' (\code{meanlog}, \code{sdlog}, \code{shift}) follow the functions of the same
#' name in package \pkg{brms}, so that existing parameter vectors keep working.
#'
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of random values.
#' @param meanlog Mean of the underlying normal distribution on the log scale.
#' @param sdlog Standard deviation of the underlying normal distribution on the log scale.
#' @param shift Lower bound of the distribution.
#'
#' @return \code{dshifted_lnorm} gives the density, \code{pshifted_lnorm} the
#'   distribution function, \code{qshifted_lnorm} the quantile function and
#'   \code{rshifted_lnorm} generates random values.
#'
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#' @keywords internal
NULL

#' @rdname shifted_lnorm
#' @export
dshifted_lnorm <- function(x, meanlog = 0, sdlog = 1, shift = 0, log = FALSE) {
  stats::dlnorm(x - shift, meanlog = meanlog, sdlog = sdlog, log = log)
}

#' @rdname shifted_lnorm
#' @export
pshifted_lnorm <- function(q, meanlog = 0, sdlog = 1, shift = 0,
                           lower.tail = TRUE, log.p = FALSE) {
  stats::plnorm(q - shift, meanlog = meanlog, sdlog = sdlog,
                lower.tail = lower.tail, log.p = log.p)
}

#' @rdname shifted_lnorm
#' @export
qshifted_lnorm <- function(p, meanlog = 0, sdlog = 1, shift = 0,
                           lower.tail = TRUE, log.p = FALSE) {
  shift + stats::qlnorm(p, meanlog = meanlog, sdlog = sdlog,
                        lower.tail = lower.tail, log.p = log.p)
}

#' @rdname shifted_lnorm
#' @export
rshifted_lnorm <- function(n, meanlog = 0, sdlog = 1, shift = 0) {
  shift + stats::rlnorm(n, meanlog = meanlog, sdlog = sdlog)
}
