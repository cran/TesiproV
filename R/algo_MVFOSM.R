#' @name MVFOSM
#' @title Mean-Value First-Order Second-Moment (MVFOSM)
#' @description
#' Mean-Value First-Order Second-Moment (MVFOSM) method for the
#' approximation of failure probabilities in structural reliability analysis.
#'
#' The limit-state function is linearised at the mean values of the basic
#' variables. The reliability index is
#' \deqn{\beta = \frac{g(\mu)}{\sqrt{\sum_i \left(\frac{\partial g}{\partial x_i}\Big|_{\mu} \sigma_i\right)^2}}}
#' and the probability of failure is approximated by \eqn{p_f = \Phi(-\beta)}.
#'
#' Failure is defined as \eqn{g(x) < 0}. If \eqn{g(\mu) < 0}, the
#' reliability index is negative.
#'
#' Only the mean values and standard deviations of the basic variables are
#' used. The distribution types do not influence the result. The method is
#' exact for linear limit-state functions of normal variables. For nonlinear
#' limit-state functions the result depends on how the limit-state function is
#' formulated (e.g. \eqn{R - S} or \eqn{R/S - 1} give different values of
#' \eqn{\beta} for the same failure domain).
#'
#' @details
#' When MVFOSM is run through a \code{PROB_MACHINE} (\code{fCall = "MVFOSM"}),
#' \code{lsf} and \code{lDistr} are supplied automatically by \code{SYS_PROB}
#' or \code{SYS_PARAM} and must not be part of the machine's \code{options}.
#' Only \code{h}, \code{isExpression} and \code{debug.level} are passed as
#' options. With \code{isExpression = TRUE} the system passes the symbolic
#' expression of the limit-state function instead of the function.
#'
#' The mean values and standard deviations are taken from the distribution
#' objects (fields \code{mean} and \code{Sd}), not from the native distribution
#' parameters. All basic variables must therefore have a finite mean and
#' standard deviation. This excludes for example the log-Student-t
#' distribution (\code{"lt"}), for which the moments are not defined. In this
#' case the function stops with an error.
#'
#' With \code{isExpression = FALSE} the partial derivatives are computed by
#' central finite differences with the step \eqn{h \sigma_i} for variable
#' \eqn{i}. The step is relative to the standard deviation, so that variables
#' of very different magnitude (e.g. a section modulus in m^3 and a yield
#' strength in kN/m^2) are handled equally. For a variable with
#' \eqn{\sigma_i = 0} the step \eqn{h \max(|\mu_i|, 1)} is used. With
#' \code{isExpression = TRUE} the derivatives are computed symbolically with
#' \code{\link[stats]{D}}.
#'
#' @param lsf Limit-state function, either a function or an expression,
#'   as indicated by \code{isExpression}. Supplied automatically by a SYS_
#'   object (SYS_PROB or SYS_PARAM) when the method is run through a
#'   PROB_MACHINE. Do not provide it in the machine's options.
#' @param lDistr List of distributions as returned by \code{getlDistr()} of
#'   the basic variables. Each variable must have a finite mean and standard
#'   deviation. Supplied automatically as well.
#' @param h Step of the central finite differences, relative to the standard
#'   deviation of each variable (default 1e-4). Only used if
#'   \code{isExpression = FALSE}.
#' @param isExpression Logical. If \code{TRUE}, \code{lsf} must be an
#'   expression and the derivatives are computed symbolically. Otherwise
#'   \code{lsf} must be a function.
#' @param debug.level Verbosity level: 0 = silent, 1 = basic info,
#'   2 = detailed output.
#'
#' @return MVFOSM returns a list with the following elements:
#' \itemize{
#'   \item \code{method}: Character string "MVFOSM".
#'   \item \code{beta}: Reliability index.
#'   \item \code{pf}: Probability of failure \eqn{\Phi(-\beta)}.
#'   \item \code{design.point.x}: Approximate design point in the original
#'         \eqn{x}-space, \eqn{\mu_i + \alpha_i \beta \sigma_i}.
#'   \item \code{alpha}: Direction cosines (sensitivity factors).
#'   \item \code{runtime}: Runtime as returned by \code{proc.time()}.
#' }
#' @references
#' Freudenthal, A. M. (1956).
#' Safety and the probability of structural failure.
#' \emph{Transactions of the American Society of Civil Engineers},
#' 121, 1337-1397.
#'
#' Nowak, A. S., & Collins, K. R. (2000).
#' \emph{Reliability of Structures}. McGraw-Hill.
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#' @export

MVFOSM <- function(lsf,
                   lDistr,
                   h = 0.0001,
                   isExpression = FALSE,
                   debug.level = 0) {
  debug.TAG <- "MVFOSM"
  debug.print(debug.level, debug.TAG, c(TRUE),
              msg = "Mean Value First Order Second Moment Method started...")
  tic <- proc.time()

  n_vars <- length(lDistr)
  m <- numeric(n_vars)
  sd <- numeric(n_vars)
  dx <- numeric(n_vars)

  # mean and standard deviation from the distribution object.
  # lDistr[[i]][[2]] holds the native distribution parameters
  # (e.g. meanlog/sdlog for lnorm), not Mean/Sd.
  for (i in seq_len(n_vars)) {
    m[i] <- as.numeric(lDistr[[i]][[1]]$mean %||% NA_real_)
    sd[i] <- as.numeric(lDistr[[i]][[1]]$Sd %||% NA_real_)
  }
  if (any(!is.finite(m)) || any(!is.finite(sd))) {
    stop("MVFOSM requires finite Mean and Sd for all basic variables ",
         "(not available e.g. for 'lt').", call. = FALSE)
  }

  if (isExpression) {
    # first assign all variables, then differentiate
    for (i in 1:n_vars) assign(lDistr[[i]][[1]]$name, m[i])
    s <- 0
    for (i in 1:n_vars) {
      d <- stats::D(lsf, lDistr[[i]][[1]]$name)
      dx[i] <- eval(d)
      s <- s + (dx[i]^2) * (sd[i]^2)
    }
    sum_sd <- sqrt(s)
    lsf_mean <- eval(lsf)
  } else {
    # central differences with a step relative to the standard deviation:
    # the error of dx[i] * sd[i] is then about eps * |g| / h, independent
    # of the scale of the variable (e.g. Z = 1e-4 m^3 or Fy = 6e5 kN/m^2)
    for (i in 1:n_vars) {
      step <- h * if (sd[i] > 0) sd[i] else max(abs(m[i]), 1)
      x_lo <- m
      x_hi <- m
      x_lo[i] <- m[i] - step
      x_hi[i] <- m[i] + step
      dx[i] <- (lsf(x_hi) - lsf(x_lo)) / (2 * step)
    }
    sum_sd <- sqrt(sum((dx^2) * (sd^2)))
    lsf_mean <- lsf(m)
  }

  beta <- lsf_mean / sum_sd
  pf <- stats::pnorm(-beta, 0, 1)

  alpha <- -(dx * sd) / sum_sd
  designPoint <- m + alpha * beta * sd

  duration <- proc.time() - tic

  output <- list(
    "method" = "MVFOSM",
    "beta" = beta,
    "pf" = pf,
    "design.point.x" = designPoint,
    "alpha" = alpha,
    "runtime" = duration[1:5]
  )

  debug.print(debug.level, debug.TAG, c(duration),
              msg = "Mean Value First Order Second Moment Method finished in: ")
  return(output)
}
