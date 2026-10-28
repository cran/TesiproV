# =============================================================================
#  OC_helpers.R
#  Internal helper functions for OC-curve based Bayesian updating.
#  NOT exported - called exclusively by called by PROB_QUALCON and PROB_BASEVAR$updateBasevar().
#
#  Based on original implementations by TesiproV authors.
#  Refactored for OOP integration and CRAN compliance:
#    - library() calls removed (CRAN policy violation)
#    - Rsolnp::solnp() replaced by a lightweight numerical root search
#    - while loops in Monte-Carlo sections replaced by vectorised operations
# =============================================================================


#
# # =============================================================================
# #  Internal helper - extract full 4-element hyperparameter vector
# # =============================================================================

#' @keywords internal
# .extract_hyper <- function(var) {
#   hp <- var$DistributionParameters
#
#   if (length(hp) < 4) {
#     cache <- var$.cache
#     if (!is.null(cache)) {
#       for (key in c("hyper_posterior", "hyper_prior")) {
#         if (exists(key, envir = cache, inherits = FALSE)) {
#           hp <- get(key, envir = cache)
#           break
#         }
#       }
#     }
#     if (length(hp) < 4) {
#       hp <- c(hp, rep(NA_real_, 4 - length(hp)))
#       hp[3] <- 3; hp[4] <- 10
#       message("Note: n' and v' not found - using JCSS defaults (n'=3, v'=10).")
#     }
#   }
#
#   names(hp) <- c("m'", "s'", "n'", "v'")
#   return(hp)
# }


# -----------------------------------------------------------------------------
#  .oc_check_limits
#  Checks the defect limits (DefectLimitLower / DefectLimitUpper). They define
#  when a single unit counts as defective. The acceptance criterion of the
#  sampling plan is independent of them.
#    fck_EN-206: DefectLimitLower = f_ck (finite), DefectLimitUpper = Inf
#    dim_k||n  : arbitrary limits, at least one finite, Lower < Upper
# -----------------------------------------------------------------------------
.oc_check_limits <- function(conf.type, lower, upper) {
  if (length(lower) != 1L || length(upper) != 1L || is.na(lower) || is.na(upper)) {
    stop("DefectLimitLower and DefectLimitUpper must each be a single number ",
         "(-Inf / Inf for an open side).", call. = FALSE)
  }
  if (conf.type == "fck_EN-206") {
    if (!is.finite(lower)) {
      stop("For conf.type = 'fck_EN-206', DefectLimitLower (= f_ck) must be set to a finite value. ",
           "This also applies to a manually specified OC curve.", call. = FALSE)
    }
    if (is.finite(upper)) {
      stop("For conf.type = 'fck_EN-206' only DefectLimitLower (= f_ck) is allowed, ",
           "DefectLimitUpper must remain Inf.", call. = FALSE)
    }
  } else {
    if (lower >= upper) {
      stop("For conf.type = 'dim_k||n', DefectLimitLower < DefectLimitUpper must hold",
           call. = FALSE)
    }
    if (!is.finite(lower) && !is.finite(upper)) {
      stop("For conf.type = 'dim_k||n', at least one limit ",
           "(DefectLimitLower or DefectLimitUpper) must be finite.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
#  .oc_compute_pa
#  Computes acceptance probabilities P_a for a vector of defect fractions theta.
#  Faithfully implements the original get_OC() logic, vectorised where possible.
#
#  @param theta        Numeric vector of defect fractions (0 < theta <= 1)
#  @param conf.type   "fck_EN-206" or "dim_k||n"
#  @param conf.param  c(n) for "fck_EN-206", c(n, c_accept) for "dim_k||n"
#  @param dist.type   "norm", "st" or "lt"
#  @param dist.param  Hyperparameter vector c(m', s', n', v')
#  @param lower       Lower defect limit. For "fck_EN-206" this is f_ck.
#  @param upper       Upper defect limit (must be Inf for "fck_EN-206").
#  @param N           Monte-Carlo sample size per defect fraction (only used
#                     for "fck_EN-206", default 100 000).
#
#  @return Numeric vector of P_a values (same length as theta)
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#  .oc_compute_pa
# -----------------------------------------------------------------------------
.oc_compute_pa <- function(theta,
                           conf.type,
                           conf.param,
                           dist.type,
                           dist.param,
                           lower = -Inf,
                           upper = Inf,
                           N = 100000
) {

  P_a <- rep(0, length(theta))

  if (conf.type == "dim_k||n") {
    if (length(conf.param) < 2 || anyNA(conf.param[1:2])) {
      stop("conf.param for conf.type = 'dim_k||n' must contain c(n, c_accept).",
           call. = FALSE)
    }
    for (j in seq_along(theta))
      P_a[j] <- stats::pbinom(conf.param[2], size = conf.param[1], prob = theta[j])
    return(P_a)
  }
  else if (conf.type == "fck_EN-206") {
    if (length(conf.param) < 1 || is.na(conf.param[1])) {
      stop("conf.param for conf.type = 'fck_EN-206' must contain n: c(n).",
           call. = FALSE)
    }
    if (length(conf.param) > 1) {
      stop("conf.param for conf.type = 'fck_EN-206' must be c(n). ",
           "The Monte-Carlo sample size is given by the argument N.",
           call. = FALSE)
    }
    if (length(N) != 1L || !is.finite(N) || N < 1) {
      stop("N must be a single positive number.", call. = FALSE)
    }
    .oc_check_limits(conf.type, lower, upper)
    n    <- conf.param[1]
    f_ck <- lower

    if (dist.type == "norm") {
      for (j in seq_along(theta)) {
        m      <- f_ck + stats::qnorm(1 - theta[j]) * dist.param[2]
        x      <- matrix(stats::rnorm(N * n, m, dist.param[2]), nrow = N)
        P_a[j] <- mean(rowMeans(x) >= f_ck + 1.48 * dist.param[2] &
                         apply(x, 1, min) >= f_ck - 4)
      }

    } else if (dist.type == "st") {
      s    <- dist.param[2] / sqrt(dist.param[3] / (dist.param[3] + 1))
      sd   <- s * sqrt(dist.param[4] / (dist.param[4] - 2))
      for (j in seq_along(theta)) {
        m_mod  <- f_ck + stats::qt(1 - theta[j], dist.param[4]) * s
        x      <- matrix(stats::rt(N * n, dist.param[4]) * s + m_mod, nrow = N)
        P_a[j] <- mean(rowMeans(x) >= f_ck + 1.48 * sd &
                         apply(x, 1, min) >= f_ck - 4)
      }

    } else if (dist.type == "lnorm") {
      sig  <- sqrt(log((dist.param[2] / dist.param[1])^2 + 1))
      sd   <- dist.param[2]
      for (j in seq_along(theta)) {
        f_root  <- function(sig_v) sd^2 - (exp(sig_v^2) - 1) * f_ck^2 *
          exp(2 * sig_v * stats::qnorm(1 - theta[j]) + sig_v^2)
        sig_mod <- tryCatch(
          as.numeric(stats::uniroot(f_root, c(0.001, sd), tol = 1e-7)$root),
          error = function(e) sig)
        mue_mod <- log(f_ck) + stats::qnorm(1 - theta[j]) * sig_mod
        x       <- matrix(stats::rlnorm(N * n, mue_mod, sig_mod), nrow = N)
        P_a[j]  <- mean(rowMeans(x) >= f_ck + 1.48 * sd &
                          apply(x, 1, min) >= f_ck - 4)
      }

    } else if (dist.type == "lt") {
      .get_sd_lt <- function(hp, n_sim = 50000L) {
        x <- rlt(n_sim, hp)
        x <- x[is.finite(x)]
        if (length(x) < 100) return(NA_real_)
        stats::sd(x)
      }

      .get_s_mod <- function(m_mod, hp, tol = 0.01, max_iter = 50L) {
        sd_goal <- .get_sd_lt(hp)
        if (!is.finite(sd_goal) || sd_goal <= 0)
          stop("'.get_s_mod': sd_goal is not finite/positive - check the hyperparameters (v' > 2 required).",
               call. = FALSE)

        hp2    <- hp[2]
        sd_new <- .get_sd_lt(c(m_mod, hp2, hp[3], hp[4]))

        iter <- 0L
        while (is.finite(sd_new) && abs(sd_new - sd_goal) > tol && iter < max_iter) {
          hp2    <- hp2 * abs(sd_goal / sd_new)
          if (!is.finite(hp2) || hp2 <= 0) break
          sd_new <- .get_sd_lt(c(m_mod, hp2, hp[3], hp[4]))
          iter   <- iter + 1L
        }

        if (!is.finite(sd_new) || !is.finite(hp2) || hp2 <= 0) {
          warning(sprintf("'.get_s_mod': no convergence for m_mod = %.4f - using hp[2] as fallback.", m_mod),
                  call. = FALSE)
          hp2 <- hp[2]
        } else if (iter >= max_iter) {
          warning(sprintf("'.get_s_mod': max_iter (%d) reached for m_mod = %.4f, |sd_new - sd_goal| = %.4f.",
                          max_iter, m_mod, abs(sd_new - sd_goal)),
                  call. = FALSE)
        }

        hp2 / sqrt(hp[3] / (hp[3] + 1))
      }
      sd   <- .get_sd_lt(dist.param)
      s    <- dist.param[2] / sqrt(dist.param[3] / (dist.param[3] + 1))
      for (j in seq_along(theta)) {
        m_mod     <- log(f_ck) + stats::qt(1 - theta[j], dist.param[4]) * s
        s_mod     <- .get_s_mod(m_mod, dist.param)
        hyper_mod <- c(m_mod, s_mod * sqrt(dist.param[3]/(dist.param[3]+1)),
                       dist.param[3], dist.param[4])
        x      <- matrix(rlt(N * n, hyper_mod), nrow = N)
        P_a[j] <- mean(rowMeans(x) >= f_ck + 1.48 * sd &
                         apply(x, 1, min) >= f_ck - 4)
      }
    }
    return(P_a)
  }
  else {
    stop("conf.type '", conf.type, "' not known, P_a not computable")
  }
}


# -----------------------------------------------------------------------------
#  .oc_get_outgoing_pdf
# -----------------------------------------------------------------------------
.oc_get_outgoing_pdf <- function(theta, P_a, conf.type, conf.param, dist.type, dist.param,
                                 lower = -Inf, upper = Inf) {

  lM <- 100L; lS <- 100L; lX <- 100L

  u_t <- stats::qt(theta, dist.param[4])
  fun_P <- function(x) {
    if (x > u_t[length(u_t)]) return(P_a[length(P_a)])
    if (x < u_t[1]) return(P_a[1])

    lo <- max(which(u_t <= x))
    hi <- lo + 1L

    if (is.infinite(u_t[lo]) || is.infinite(u_t[hi]))
      return(P_a[lo])

    P_a[lo] + (x - u_t[lo]) / (u_t[hi] - u_t[lo]) * (P_a[hi] - P_a[lo])
  }

  if (conf.type %in% c("fck_EN-206", "dim_k||n")) {
    # Limits that define a defect. For fck_EN-206 the OC curve was computed
    # with f_ck = DefectLimitLower. The limits must also be set if the OC
    # curve is specified manually.
    .oc_check_limits(conf.type, lower, upper)
    lim_lo <- lower
    lim_hi <- upper

    # Acceptance probability per grid point (M, S). 'scale' is the scaling
    # factor of the distribution family (st: sqrt(v/(v-2)), lt: 1).
    P_grid <- function(M, S, scale) {
      U_lo <- (lim_lo - M) %*% (1 / t(S)) / scale
      if (conf.type == "fck_EN-206") {
        return(apply(U_lo, c(1, 2), fun_P))
      }
      U_hi <- (lim_hi - M) %*% (1 / t(S)) / scale
      th   <- stats::pt(U_lo, dist.param[4]) + 1 - stats::pt(U_hi, dist.param[4])
      matrix(stats::approx(theta, P_a, xout = as.vector(th), rule = 2,
                           ties = "ordered")$y, nrow = nrow(th))
    }

    if (dist.type == "st") {
      m    <- dist.param[1]
      s    <- dist.param[2] / sqrt(dist.param[3] / (dist.param[3] + 1))
      M_incr <- (6 * s) / lM; S_incr <- (2 * s) / lS; X_incr <- (12 * s) / lX
      M  <- matrix(seq(m - 3 * s, m + 3 * s, by = M_incr), ncol = 1)
      S  <- matrix(seq(0.001, 2 * s, by = S_incr), ncol = 1)
      X  <- matrix(seq(m - 6 * s, m + 6 * s, by = X_incr), ncol = 1)
      P  <- P_grid(M, S, sqrt(dist.param[4] / (dist.param[4] - 2)))
      Z2 <- (M - m) %*% (1 / t(S)) * sqrt(dist.param[3])
      f_strich_MS <- sweep(
        sweep(exp(-0.5 * Z2^2), 2, sqrt(dist.param[3]) / sqrt(2 * pi) * (1 / t(S)), "*"), 2,
        (0.5 * dist.param[4] * (dist.param[2] * (1 / t(S)))^2)^(dist.param[4] / 2 - 1) *
          exp(-0.5 * dist.param[4] * (dist.param[2] * (1 / t(S)))^2) /
          gamma(dist.param[4] / 2) * dist.param[4] * dist.param[2]^2 * (1 / t(S))^3, "*")
      f_2strich_MS <- f_strich_MS * P
      f_2strich_MS <- f_2strich_MS / (sum(f_2strich_MS) * M_incr * S_incr)
      fun_f_o_star <- function(x) {
        Z1 <- (x - M) %*% (1 / t(S))
        sum(sweep(stats::dnorm(Z1), 2, 1 / t(S), "*") * f_2strich_MS)
      }

    } else if (dist.type == "lt") {
      m_n  <- exp(dist.param[1] + dist.param[2]^2 / 2)
      s_n  <- sqrt((exp(dist.param[2]^2) - 1) * exp(2 * dist.param[1] + dist.param[2]^2))

      # Grid in the positive range only: m_n - 8 s_n becomes negative for
      # CoV > 0.125, where log(M) and log(x) are not defined.
      x_lo <- max(m_n - 8 * s_n, 0.05 * m_n)
      x_hi <- m_n + 8 * s_n
      s_hi <- 8 * s_n
      M  <- matrix(seq(x_lo, x_hi, length.out = lM + 1L), ncol = 1)
      S  <- matrix(seq(0.001, s_hi, length.out = lS + 1L), ncol = 1)
      X  <- matrix(seq(x_lo, x_hi, length.out = lX + 1L), ncol = 1)
      M_incr <- (x_hi - x_lo) / lM
      S_incr <- (s_hi - 0.001) / lS
      X_incr <- (x_hi - x_lo) / lX

      P  <- P_grid(M, S, 1)
      delta   <- t(S %*% (1 / t(M)))
      sig_lnX <- sqrt(log(delta^2 + 1))
      Z2 <- sweep(1 / sig_lnX * sqrt(dist.param[3]), 1, log(M) - dist.param[1], "*")
      gamma_term <- (0.5 * dist.param[4] * (dist.param[2] / sig_lnX)^2)^(dist.param[4] / 2 - 1) *
        exp(-0.5 * dist.param[4] * (dist.param[2] / sig_lnX)^2) /
        gamma(dist.param[4] / 2) * dist.param[4] * dist.param[2]^2 / sig_lnX^3

      f_strich_MS <- exp(-0.5 * Z2^2) * sqrt(dist.param[3]) / sqrt(2 * pi) *
        (1 / sweep(sig_lnX, 1, M, "*")) * gamma_term

      f_2strich_MS <- f_strich_MS * P

      fun_f_o_star <- function(x) {
        mue_lnX <- log(sweep(1 / sqrt(outer(M[, 1]^2, S[, 1]^2, "+")), 1, M^2, "*"))
        Z1 <- (log(x) - mue_lnX) / sig_lnX
        sum((stats::dnorm(Z1) / sig_lnX / x) * f_2strich_MS)
      }

    } else {
      warning("dist.type '", dist.type, "' not known"); return(NULL)
    }

  } else {
    warning("conf.type '", conf.type, "' not supported"); return(NULL)
  }

  f_disc <- apply(X, 1, fun_f_o_star)
  f_o    <- f_disc / sum(f_disc) / X_incr
  m_o    <- sum(X * f_o * X_incr)
  s_o    <- sqrt(sum((X - m_o)^2 * f_o * X_incr))
  list(X = as.numeric(X), f_o = f_o, m_o = m_o, s_o = s_o,
       f_2strich_MS = f_2strich_MS, M = M, S = S)
}


# # -----------------------------------------------------------------------------
# #  .oc_fit_posterior_hyper
# # -----------------------------------------------------------------------------
# .oc_fit_posterior_hyper <- function(f_o_result, dist.type = "t") {
#   if (dist.type != "t") {
#     warning("dist.type '", dist.type, "' not yet supported"); return(NULL)
#   }
#   M <- f_o_result$M; S <- f_o_result$S
#   f_2strich_MS <- f_o_result$f_2strich_MS
#   m_o <- f_o_result$m_o; s_o <- f_o_result$s_o
#
#   fun_f <- function(mue, sig, hp) {
#     m <- hp[1]; s <- hp[2]; n <- hp[3]; v <- hp[4]
#     sqrt(n)/sqrt(2*pi)/sig * exp(-0.5*((mue-m)/(sig/sqrt(n)))^2) *
#       (0.5*v*(s/sig)^2)^(v/2-1) * exp(-0.5*v*(s/sig)^2) /
#       gamma(v/2) * v*s^2/sig^3
#   }
#   fit_fn <- function(hp, mu_v, sg_v, y)
#     sum((y - outer(mu_v, sg_v, Vectorize(function(mu,sg) fun_f(mu,sg,hp))))^2)
#
#   fit        <- stats::optim(c(m_o,s_o,10,10), fit_fn,
#                              mu_v=as.numeric(M), sg_v=as.numeric(S), y=f_2strich_MS)
#   hp_post        <- fit$par
#   names(hp_post) <- c("m''","s''","n''","v''")
#   return(hp_post)
# }

#' @title Fit parametric posterior hyperparameters to the outgoing distribution
#' @description Tries to fit a Student-t ("st") or log-Student-t ("lt")
#'   distribution to the discretised outgoing PDF f_o(X) by least squares.
#'   Returns NULL if the fit does not converge or if its quality (R^2) is
#'   below min_r2. The caller then falls back to the empirical distribution.
#' @keywords internal
.oc_fit_posterior_hyper <- function(f_o_result, dist.type, min_r2) {
  if (!dist.type %in% c("st", "lt")) {
    warning("dist.type '", dist.type, "' not supported for parametric fit.")
    return(NULL)
  }

  X   <- f_o_result$X
  f_o <- f_o_result$f_o
  m_o <- f_o_result$m_o
  s_o <- f_o_result$s_o

  d_fun <- if (dist.type == "st") dst else dlt

  if (dist.type == "st") {
    hp0 <- c(m_o, s_o, 5, 10)
  } else {
    sig0 <- sqrt(log(1 + (s_o / m_o)^2))
    mu0  <- log(m_o) - sig0^2 / 2
    hp0  <- c(mu0, sig0, 5, 10)
  }

  keep <- f_o > 1e-8 * max(f_o)   # ignore negligible mass

  fit_fn <- function(hp) {
    if (hp[3] <= 0 || hp[4] <= 2) return(1e10)   # n' > 0, v' > 2
    pred <- tryCatch(d_fun(X[keep], hp), error = function(e) rep(NA_real_, sum(keep)))
    if (anyNA(pred)) return(1e10)
    sum((f_o[keep] - pred)^2)
  }

  fit <- tryCatch(
    stats::optim(hp0, fit_fn, method = "Nelder-Mead",
                 control = list(maxit = 2000, reltol = 1e-10)),
    error = function(e) NULL
  )

  if (is.null(fit) || fit$convergence != 0) {
    message("Parametric fit did not converge - using the empirical distribution.")
    return(NULL)
  }

  hp_post <- fit$par
  pred    <- d_fun(X[keep], hp_post)
  ss_res  <- sum((f_o[keep] - pred)^2)
  ss_tot  <- sum((f_o[keep] - mean(f_o[keep]))^2)
  r2      <- 1 - ss_res / ss_tot

  if (!is.finite(r2) || r2 < min_r2) {
    message(sprintf("Parametric fit insufficient (R^2 = %.4f < %.4f) - using the empirical distribution.",
                    r2, min_r2))
    return(NULL)
  }

  names(hp_post) <- c("m'", "s'", "n'", "v'")
  attr(hp_post, "r2") <- r2
  hp_post
}
