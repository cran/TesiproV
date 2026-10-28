#' Evaluate one chunk in single-limit-state importance sampling
#'
#' Transforms a chunk of samples from standard-normal space to physical
#' space, evaluates the limit-state function and calculates the log-domain
#' sums required by the single-limit-state Monte Carlo importance-sampling
#' estimator. This helper is intended to be called by \code{MC_IS_single()}
#' through the parallel dispatcher.
#'
#' The function receives each sample chunk and its dependencies as explicit
#' arguments. It is defined at package level so that it does not need to
#' capture the simulation's enclosing environment in a local closure.
#'
#' @param u_chunk Numeric matrix containing the samples for this worker,
#'   with one sample per row and one random variable per column.
#' @param lDistr_fun List of marginal distribution objects used to transform
#'   the samples to physical space.
#' @param u_to_x_worker Function that transforms standard-normal values to
#'   physical-space values for one marginal distribution.
#' @param dp Numeric vector giving the importance-sampling proposal centre
#'   in standard-normal space.
#' @param Sigma_chol_inv Inverse Cholesky factor used to evaluate the
#'   importance-sampling proposal density.
#' @param log_det_Sigma_chol Log determinant term for the proposal-density
#'   calculation.
#' @param n_vars Number of random variables.
#' @param lsf Limit-state function evaluated on the transformed samples.
#' @param lsf_vectorized Logical flag indicating whether \code{lsf} accepts
#'   matrix input and returns one value per sample.
#' @param debug.level Integer verbosity level: 0 = silent, 1 = basic
#'   information, 2 = detailed diagnostics.
#'
#' @return A list containing the log-domain sums for the chunk:
#'   \describe{
#'     \item{\code{log_sum_Iw}}{Log of the sum of importance weights for
#'       samples in the failure domain.}
#'     \item{\code{log_sum_Iw2}}{Log of the sum of squared importance
#'       weights for samples in the failure domain.}
#'     \item{\code{log_sum_w}}{Log of the sum of all importance weights.}
#'     \item{\code{log_sum_w2}}{Log of the sum of squared importance
#'       weights.}
#'     \item{\code{n}}{Number of samples processed in the chunk.}
#'   }
#'
#' @keywords internal
#' @noRd
.mc_is_single_worker <- function(
  u_chunk,
  lDistr_fun,
  u_to_x_worker,
  dp,
  Sigma_chol_inv,
  log_det_Sigma_chol,
  n_vars,
  lsf,
  lsf_vectorized,
  debug.level
) {
    log_sum_exp_local <- function(x) {
        if (length(x) == 0L) {
            return(-Inf)
        }

        m <- max(x)

        if (!is.finite(m)) {
            return(m)
        }

        m + log(sum(exp(x - m)))
    }

    # n_local <- nrow(u_chunk)

    # --------------------------------------------------
    # Transform u -> x (vectorized)
    # --------------------------------------------------

    # code is slow
    # p_mat <- pnorm(u_chunk)
    # p_mat <- pmin(pmax(p_mat, 1e-14), 1 - 1e-14)

    # q_list <- lapply(lDistr, function(d) d$q)

    # x_mat <- matrix(NA_real_, n_local, n_vars)

    # for (j in seq_len(n_vars)) {
    #   x_mat[, j] <- q_list[[j]](p_mat[, j])
    # }

    # Fix

    x_cols <- vector("list", n_vars)

    for (j in seq_len(n_vars)) {
        x_cols[[j]] <- u_to_x_worker(
            u = u_chunk[, j],
            marginal = lDistr_fun[[j]]
        )
    }

    x_mat <- do.call(cbind, x_cols)

    # --------------------------------------------------
    # log phi(u)  (standard normal)
    # --------------------------------------------------

    log_phi <- -0.5 * rowSums(u_chunk^2) -
        0.5 * n_vars * log(2 * pi)

    diff <- sweep(u_chunk, 2, dp, "-")

    v <- diff %*% Sigma_chol_inv

    quad <- rowSums(v * v)

    log_h <- -0.5 * quad -
        0.5 * n_vars * log(2 * pi) -
        log_det_Sigma_chol

    log_w <- log_phi - log_h
    # cat("mean(exp(log_fx - log_phi)) =", mean(exp(log_fx - log_phi)), "\n")

    # debug print
    # cat("mean(log_w) =", mean(log_w), "\n")
    # cat("var(log_w)  =", var(log_w), "\n")
    # cat("mean(w) =", mean(exp(log_w)), "\n")

    # --------------------------------------------------
    # Vectorized LSF evaluation
    # --------------------------------------------------

    # IMPORTANT: lsf must accept matrix input
    # g_vals <- lsf(x_mat)
    # g_vals <- apply(x_mat, 1, lsf)
    if (lsf_vectorized) {
        g_vals <- lsf(x_mat)
    } else {
        g_vals <- apply(x_mat, 1, lsf)
    }

    I_vals <- as.numeric(g_vals < 0)


    if (debug.level >= 2) {
        message(sprintf(
            "[MC_IS_single] u range: [%.3f, %.3f]",
            min(u_chunk),
            max(u_chunk)
        ))

        message(sprintf(
            "[MC_IS_single] x range: [%.6f, %.6f]",
            min(x_mat),
            max(x_mat)
        ))

        message(sprintf(
            "[MC_IS_single] g range: [%.6f, %.6f]",
            min(g_vals),
            max(g_vals)
        ))

        message(sprintf(
            "[MC_IS_single] failures: %d / %d",
            sum(I_vals),
            length(I_vals)
        ))

        message(sprintf(
            "[MC_IS_single] log_w range: [%.3f, %.3f]",
            min(log_w),
            max(log_w)
        ))
    }

    # Debug print
    # cat("mean(I*w) =", mean(exp(log_w) * I_vals), "\n")
    # cat("mean(u_chunk) =", colMeans(u_chunk), "\n")

    # list(
    #   log_w  = log_w,
    #   I_vals = I_vals
    # )


    # --------------------------------------------------
    # Stable aggregation on the log scale
    # --------------------------------------------------

    log_sum_w <- log_sum_exp_local(log_w)

    if (any(I_vals == 1L)) {
        log_sum_Iw <- log_sum_exp_local(log_w[I_vals == 1L])
    } else {
        log_sum_Iw <- -Inf
    }

    log_sum_w2 <- log_sum_exp_local(2 * log_w)

    # second moment of I * w for the variance of the estimator
    if (any(I_vals == 1L)) {
        log_sum_Iw2 <- log_sum_exp_local(2 * log_w[I_vals == 1L])
    } else {
        log_sum_Iw2 <- -Inf
    }

    list(
        log_sum_Iw = log_sum_Iw,
        log_sum_Iw2 = log_sum_Iw2,
        log_sum_w = log_sum_w,
        log_sum_w2 = log_sum_w2,
        n = length(log_w)
    )
}
