# Utility helpers
#' Ensure that a required package is available
#'
#' @param pkg Character string with the package name
#' @return Invisible TRUE (used only for its side-effect)
#' @note The function stops with an informative error if the required
#'   package is not installed.
#' @keywords internal
ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf(
        "Package \"%s\" is required but not installed.\n",
        pkg
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}


#' Smooth approximation of the minimum function
#'
#' Computes a differentiable approximation of \code{min(x)}
#' using a log-sum-exp formulation:
#'
#' \deqn{
#'   \min(x) \approx -\frac{1}{\kappa} \log\left( \sum_i e^{-\kappa x_i} \right)
#' }
#'
#' The approximation becomes exact as \code{kappa -> Inf}.
#'
#' @param x Numeric vector.
#' @param kappa Positive smoothing parameter controlling
#'   the sharpness of the approximation (default = 20).
#'
#' @return Numeric scalar.
#'
#' @details
#' This function is used internally for smooth system
#' limit-state aggregation (serial systems) to improve
#' numerical stability and differentiability.
#'
#' @keywords internal
#' @noRd
smooth_min <- function(x, kappa = 20) {
  # shift by the minimum to avoid overflow of exp() for large |x|
  m <- min(x)
  m - log(sum(exp(-kappa * (x - m)))) / kappa
}

#' Smooth approximation of the maximum function
#'
#' Computes a differentiable approximation of \code{max(x)}
#' using a log-sum-exp formulation:
#'
#' \deqn{
#'   \max(x) \approx \frac{1}{\kappa} \log\left( \sum_i e^{\kappa x_i} \right)
#' }
#'
#' The approximation becomes exact as \code{kappa -> Inf}.
#'
#' @param x Numeric vector.
#' @param kappa Positive smoothing parameter controlling
#'   the sharpness of the approximation (default = 20).
#'
#' @return Numeric scalar.
#'
#' @details
#' This function is used internally for smooth system
#' limit-state aggregation (parallel systems) to improve
#' numerical stability and differentiability.
#'
#' @keywords internal
#' @noRd
smooth_max <- function(x, kappa = 20) {
  # shift by the maximum to avoid overflow of exp() for large |x|
  m <- max(x)
  m + log(sum(exp(kappa * (x - m)))) / kappa
}
