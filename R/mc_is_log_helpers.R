#' Compute failure-probability and variance estimates from log sums
#'
#' Derives the failure-probability estimate, reliability index, variance
#' estimate and coefficient of variation from accumulated importance-sampling
#' sums represented on the logarithmic scale.
#'
#' @param log_sum_Iw Logarithm of the sum of the importance weights for
#'   samples in the failure domain, \eqn{\log\left(\sum_i I_i w_i\right)}.
#'   Use \code{-Inf} when no failure-weight mass was observed.
#' @param log_sum_Iw2 Logarithm of the sum of squared importance weights for
#'   samples in the failure domain, \eqn{\log\left(\sum_i (I_i w_i)^2\right)}.
#'   Use \code{-Inf} when no failure-weight mass was observed.
#' @param n Positive finite number of samples included in the estimator.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{\code{log_pf}}{Logarithm of the estimated failure probability.
#'       It is \code{-Inf} when no failure-weight mass was observed.}
#'     \item{\code{pf}}{Estimated failure probability on the ordinary scale.
#'       It may underflow to zero even when \code{log_pf} is finite.}
#'     \item{\code{beta}}{Reliability index derived from \code{log_pf}; it
#'       may remain finite when \code{pf} underflows.}
#'     \item{\code{log_var}}{Logarithm of the estimated variance of the
#'       failure-probability estimator. It is \code{NA} when no failure-weight
#'       mass was observed and \code{-Inf} when the estimated variance is zero.}
#'     \item{\code{var}}{Estimated variance on the ordinary scale. It may
#'       underflow to zero when \code{log_var} is finite.}
#'     \item{\code{log_cov}}{Logarithm of the coefficient of variation. This
#'       value is intended for stable internal comparisons with a CoV target.}
#'     \item{\code{cov}}{Coefficient of variation on the ordinary scale. It
#'       may underflow or overflow for extreme estimates.}
#'   }
#'
#' @details
#' Let \eqn{Y_i = I_i w_i}, where \eqn{I_i} is the failure indicator and
#' \eqn{w_i} is the importance weight. The estimated probability and
#' variance are
#'
#' \deqn{
#' \hat{p}_f = \frac{1}{n}\sum_{i=1}^n Y_i,
#' \qquad
#' \widehat{\operatorname{Var}}(\hat{p}_f)
#' = \frac{1}{n}\left(
#'   \frac{1}{n}\sum_{i=1}^n Y_i^2 - \hat{p}_f^2
#' \right).
#' }
#'
#' The function evaluates these quantities using the accumulated log sums
#' and stable log-domain arithmetic. It does not change the underlying
#' variance estimator.
#'
#' If no failure-weight mass was observed, the function returns
#' \code{log_pf = -Inf}, \code{pf = 0}, \code{beta = NA_real_},
#' \code{log_var = NA_real_}, \code{var = NA_real_},
#' \code{log_cov = Inf} and \code{cov = Inf}. In this case, the variance
#' cannot be estimated from the observed failure contributions.
#'
#' @keywords internal
#' @noRd
.mc_is_log_estimates <- function(log_sum_Iw, log_sum_Iw2, n) {
    if (length(n) != 1L || !is.finite(n) || n <= 0) {
        stop("'n' must be a finite positive scalar.", call. = FALSE)
    }

    if (length(log_sum_Iw) != 1L || length(log_sum_Iw2) != 1L) {
        stop("'log_sum_Iw' and 'log_sum_Iw2' must be scalars.", call. = FALSE)
    }

    # No observed failure mass: pf is zero, but its variance is not
    # estimable from the observed failure contributions.
    if (is.infinite(log_sum_Iw) && log_sum_Iw < 0) {
        return(list(
            log_pf  = -Inf,
            pf      = 0,
            beta    = NA_real_,
            log_var = NA_real_,
            var     = NA_real_,
            log_cov = Inf,
            cov     = Inf
        ))
    }

    if (!is.finite(log_sum_Iw)) {
        return(list(
            log_pf  = NA_real_,
            pf      = NA_real_,
            beta    = NA_real_,
            log_var = NA_real_,
            var     = NA_real_,
            log_cov = NA_real_,
            cov     = NA_real_
        ))
    }

    log_pf <- log_sum_Iw - log(n)
    pf <- exp(log_pf)
    beta <- -stats::qnorm(log_pf, log.p = TRUE)

    # A finite first failure-weight sum with no finite second-moment sum
    # is not sufficient to estimate the variance.
    if (is.na(log_sum_Iw2) ||
        (is.infinite(log_sum_Iw2) && log_sum_Iw2 < 0)) {
        log_var <- NA_real_
        var <- NA_real_
        log_cov <- NA_real_
        cov <- NA_real_
    } else {
        log_second_moment <- log_sum_Iw2 - log(n)

        # Mathematically, 2 * log(pf) <= log(E[Y^2]).
        # Clamp a possible positive round-off error to zero.
        delta <- min(2 * log_pf - log_second_moment, 0)

        # log(exp(a) - exp(b)) evaluated stably as
        # a + log(1 - exp(b - a)).
        log_var <- log_second_moment +
            log(-expm1(delta)) -
            log(n)

        var <- exp(log_var)

        log_cov <- 0.5 * log_var - log_pf
        cov <- exp(log_cov)
    }

    list(
        log_pf  = log_pf,
        pf      = pf,
        beta    = beta,
        log_var = log_var,
        var     = var,
        log_cov = log_cov,
        cov     = cov
    )
}
