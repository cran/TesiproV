#' @name FORM
#' @title First-Order Reliability Method (FORM)
#' @description
#' First-Order Reliability Method (FORM) for the approximation
#' of failure probabilities in structural reliability analysis.
#'
#' The FORM estimates the probability of failure by transforming
#' the basic random variables into standard normal space
#' and approximating the limit-state function by a first-order
#' (linear) Taylor expansion at the design point.
#'
#' The reliability index is defined as the minimum distance
#' from the origin to the limit-state surface in standard normal space.
#' The probability of failure is then approximated using
#' the standard normal cumulative distribution function.
#'
#' @param lsf Objective function representing the limit-state, e.g. \code{function(R,E){R-E}}.
#'   Supplied automatically by a SYS_ object - do not provide manually.
#' @param lDistr List of distribution objects created by TesiproV. Supplied automatically by
#'   a SYS_ object - do not provide manually.
#' @param n_optim Number of optimization cycles (not required for Lagrangian algorithms).
#' @param loctol Local tolerance for convergence of the solver algorithm.
#' @param optim_type Optimization type: `"auglag"` (Augmented Lagrangian) or `"rackfies"`
#'   (Rackwitz-Fiessler iterative scheme).
#' @param debug.level Verbosity level: 0 = silent, 1 = basic info, 2 = detailed output.
#'
#' @return A list containing:
#'   * `beta` - Hasofer-Lind reliability index
#'   * `pf` - probability of failure
#'   * `x_points` - design point in physical space
#'   * `dy` - gradient vector at design point
#'
#' @references
#' Hasofer, A. M., & Lind, N. C. (1974).
#' An exact and invariant first-order reliability format.
#' \emph{Journal of the Engineering Mechanics Division, ASCE},
#' 100(1), 111-121.
#'
#' Rackwitz, R., & Fieszler, B. (1978).
#' Structural reliability under combined random load sequences.
#' \emph{Computers & Structures}, 9(5), 489-494.
#'
#' @import nloptr
#'
#' @author (C) 2021-2026 K. Nille-Hauf, J.P. Schulze-Ardey, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#' @export
#'
#'
FORM <- function(lsf, lDistr, n_optim = 10, loctol = 1e-2, optim_type = "rackfies", debug.level = 0) {
  debug.TAG <- "FORM_OP"
  debug.print(debug.level, debug.TAG, c(optim_type), msg = "FORM Algorithm started with Solvingtype:")
  tic <- proc.time()

  # options(error = recover)

  ## message("[DEBUG] typeof(lsf)=", typeof(lsf))
  ## message("[DEBUG] is.function(lsf)=", is.function(lsf))
  ## message("[DEBUG] class(lsf)=", paste(class(lsf), collapse = ","))

  # --- Default settings -------------------------------------------------
  if (is.null(loctol) || is.na(loctol)) loctol <- 1e-4
  if (is.null(n_optim) || is.na(n_optim)) n_optim <- 20L

  # Local control variables ----------------------------------------------
  n_sim <- 0L # iteration counter
  delta_x <- Inf

  n_vars <- length(lDistr)

  # --- Transformation between physical and standard Gaussian domains ------------#
  gausTrans <- function(a, toGaussian = TRUE) {
    x <- vector(length = n_vars)
    y <- vector(length = n_vars)
    b <- vector(length = n_vars)

    if (toGaussian) {
      # Entspricht Y zu X
      y <- a
      for (i in 1:n_vars) {
        p <- pmin(pmax(stats::pnorm(y[i]), 1e-12), 1 - 1e-12)
        b[i] <- lDistr[[i]][[1]]$q(p)
      }
    } else {
      # Entspricht X zu Y
      x <- a
      for (i in 1:n_vars) {
        x_ <- x[i]
        b[i] <- stats::qnorm(lDistr[[i]][[1]]$p(x_))
      }
    }
    return(b)
  }


  # lsf in gaussian domain
  h <- function(y) {
    return(lsf(gausTrans(y, TRUE)))
  }


  beta <- vector("numeric", n_optim)
  par_u <- matrix(ncol = n_vars, nrow = n_optim)
  par_x <- matrix(ncol = n_vars, nrow = n_optim)
  dy <- matrix(ncol = n_vars, nrow = n_optim)


  ###############################################################################
  ## Augmented Lagrangian variant ----------------------------------------------
  ###############################################################################
  if (optim_type == "auglag") {
    for (i in 1:n_optim) {
      # Zufallsstartvektor erzeugen
      q_rand <- stats::qnorm(stats::runif(n_vars))

      # Standart-Optimierungsfunktion
      unorm <- function(u) {
        sqrt(sum(u^2))
      }
      info.print(debug.TAG, debug.level, c("optim_run_i", "/"), c(i, n_optim))
      res <- nloptr::auglag(
        x0 = q_rand,
        fn = unorm,
        heq = h,
        localtol = loctol
      )

      beta[i] <- res$value
      par_u[i, ] <- res$par
      par_x[i, ] <- gausTrans(res$par, FALSE)
      dy[i, ] <- grad(h, res$par)
    }
    if (res$convergence < 0) {
      warning("Error, Optimization failed!")
      return()
    }

    duration <- proc.time() - tic
    output <- list(
      "beta" = min(beta),
      "pf" = stats::pnorm(-abs(min(beta))),
      "u_points" = par_u[which.min(beta), ],
      "x_points" = par_x[which.min(beta), ],
      "dy" = dy[which.min(beta), ],
      "optim_type" = optim_type,
      "runtime" = duration
    )

    ###############################################################################
    ## Rackwitz-Fiessler iterative scheme ----------------------------------------
    ###############################################################################
  } else if (identical(tolower(optim_type), "rackfies")) {
    n_sim <- 0L
    x <- numeric(n_vars)
    y <- numeric(n_vars)
    dy <- numeric(n_vars)
    m_s <- numeric(n_vars)
    sd_s <- numeric(n_vars)

    ## Safe wrappers for p/d calls ---------------------------------------
    safe_p <- function(fun, x) {
      val <- fun(x)
      if (is.na(val) || val <= 0 || val >= 1) {
        val <- min(max(val, .Machine$double.eps), 1 - .Machine$double.eps)
      }
      val
    }
    safe_d <- function(fun, x) {
      val <- fun(x)
      if (is.na(val) || val == 0) val <- .Machine$double.eps
      val
    }

    ## Numerical differentiation ----------------------------------------
    finite_diff <- function(fun_obj, x, h = 1e-5) {
      dx <- numeric(length(x))
      for (i in seq_along(x)) {
        x_lowerbound <- x
        x_upperbound <- x
        x_lowerbound[i] <- x[i] - h
        x_upperbound[i] <- x[i] + h

        # message("[DEBUG] calling lsf() at point A")
        dx_low <- fun_obj(x_lowerbound)
        # message("[DEBUG] calling lsf() at point B")
        dx_high <- fun_obj(x_upperbound)
        dx[i] <- (dx_high - dx_low) / (2 * h)
      }
      dx
    }

    ## Initial values (mean if available, otherwise median) -------------
    for (i in seq_len(n_vars)) {
      ld <- lDistr[[i]][[1]]

      if (!is.null(ld$mean) && is.finite(ld$mean)) {
        x[i] <- ld$mean
      } else {
        # fallback for heavy-tail distributions (e.g. Log-Student-t)
        x[i] <- ld$q(0.5)
        if (debug.level >= 1) {
          message(sprintf(
            "[FORM] Using median as start value for variable %d (type=%s)",
            i, ld$DistributionType
          ))
        }
      }
    }

    if (debug.level >= 2) {
      message("\n[DEBUG] Initial mean values of basic variables:")
      print(sapply(lDistr, function(ld) {
        c(
          mean = ld[[1]]$mean,
          sd = ld[[1]]$Sd,
          type = ld[[1]]$DistributionType
        )
      }))
    }
    # debug
    # print(str(lDistr[[1]][[1]]))
    # print(str(lDistr[[2]][[1]]))

    # lokale Kopie verhindern Ueberschreibung durch finite_diff()
    lsf_local <- lsf

    repeat {
      sd_s <- m_s <- y <- numeric(n_vars)

      ## Transformation: physical -> standardized -------------------------
      for (i in seq_len(n_vars)) {
        ld <- lDistr[[i]][[1]]
        x_ <- x[i]

        p_val <- safe_p(ld$p, x_)
        d_val <- safe_d(ld$d, x_)

        sd_s[i] <- stats::dnorm(stats::qnorm(p_val)) / d_val
        m_s[i] <- x[i] - sd_s[i] * stats::qnorm(p_val)

        y[i] <- (x[i] - m_s[i]) / sd_s[i]

        # Debug information
        if (debug.level >= 2) {
          message(sprintf(
            "[DEBUG][%d] p=%g d=%g sd_s=%g m_s=%g y=%g",
            i, p_val, d_val, sd_s[i], m_s[i], y[i]
          ))
        }
      }

      h_y <- lsf_local(x)

      # message("[DEBUG] calling lsf() at point D")
      dy <- finite_diff(lsf_local, x) * sd_s

      # --- Richtungscosinus alpha --------------------------------------
      alpha <- numeric(n_vars)

      for (i in seq_len(n_vars)) {
        denom <- sqrt(sum(dy^2))
        alpha[i] <- ifelse(denom == 0, -dy[i], -dy[i] / denom)
      }

      beta <- (h_y - sum(y * dy)) / sqrt(sum(dy^2))

      if (!is.finite(beta)) {
        warning(sprintf(
          "[FORM/rackfies]: Beta not finite! h_y=%s sum(y*dy)=%s",
          h_y, sum(y * dy)
        ))
        beta <- NA_real_
      }

      x_neu <- numeric(n_vars)

      for (i in seq_len(n_vars)) {
        x_neu[i] <- if (!is.na(alpha[i]) && !is.na(sd_s[i]) && !is.na(m_s[i])) {
          m_s[i] + alpha[i] * sd_s[i] * beta
        } else {
          x[i]
        }
      }

      delta_x <- sqrt(sum((x_neu - x)^2))

      x <- x_neu
      n_sim <- n_sim + 1L

      debug.print(
        debug.level, debug.TAG, c("n_sim", "beta", "delta_x", "loctol"),
        c(n_sim, beta, delta_x, loctol)
      )

      info.print(
        debug.TAG, debug.level, c("n_sim", "beta", "delta_x", "loctol"),
        c(n_sim, beta, delta_x, loctol)
      )

      if ((isTRUE(delta_x <= loctol)) || isTRUE(n_sim > n_optim)) break
    } # repeat loop end

    duration <- proc.time() - tic

    output <- list(
      method = "FORM",
      beta = beta,
      pf = stats::pnorm(-beta),
      x_points = x,
      dy = dy,
      alpha = alpha,
      diff = delta_x,
      n_sim = n_sim,
      optim_type = optim_type,
      runtime = duration[1:5]
    )
  }
  debug.print(debug.level, debug.TAG, c(duration), msg = "\nFORM Algorithm finished in  [s]: ")

  return(output)
}
