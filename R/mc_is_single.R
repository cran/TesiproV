#' Stable log-sum-exp calculation
#'
#' Computes
#'
#' \deqn{
#' \log\left(\sum_{i=1}^{n} \exp(x_i)\right)
#' }
#'
#' in a numerically stable way by subtracting the maximum input value
#' before exponentiation. This prevents numerical overflow for large
#' values and underflow for very negative values.
#'
#' For an empty input vector, the function returns \code{-Inf}, which
#' corresponds to the logarithm of a sum equal to zero.
#'
#' @param x Numeric vector.
#'
#' @return A numeric scalar containing the stable log-sum-exp value.
#'
#' @details
#' The calculation uses
#'
#' \deqn{
#' m + \log\left(\sum_{i=1}^{n} \exp(x_i - m)\right),
#' \quad m = \max_i(x_i).
#' }
#'
#' This helper is used internally by Monte-Carlo importance sampling
#' routines to accumulate likelihood weights on the logarithmic scale.
#'
#' @keywords internal
#' @noRd
log_sum_exp <- function(x) {
  if (length(x) == 0L) {
    return(-Inf)
  }

  m <- max(x)

  if (!is.finite(m)) {
    return(m)
  }

  m + log(sum(exp(x - m)))
}


#' Stable logarithmic addition
#'
#' Computes
#'
#' \deqn{
#' \log\left(\exp(a) + \exp(b)\right)
#' }
#'
#' in a numerically stable manner.
#'
#' @param a Numeric scalar representing a value on the logarithmic scale.
#' @param b Numeric scalar representing a value on the logarithmic scale.
#'
#' @return A numeric scalar equal to
#'   \code{log(exp(a) + exp(b))}, evaluated stably.
#'
#' @details
#' This function handles \code{-Inf} as the logarithm of zero:
#'
#' \deqn{
#' \log\left(\exp(-\infty) + \exp(b)\right) = b.
#' }
#'
#' It is used internally to combine log-sums across worker chunks and
#' Monte-Carlo batches without converting importance weights back to
#' the ordinary probability scale.
#'
#' @keywords internal
#' @noRd
log_add_exp <- function(a, b) {
  if (identical(a, -Inf)) {
    return(b)
  }

  if (identical(b, -Inf)) {
    return(a)
  }

  m <- max(a, b)

  m + log(exp(a - m) + exp(b - m))
}


#' Transform standard-normal values to physical space
#'
#' Transforms samples from the standard normal space to the physical
#' space of one marginal distribution.
#'
#' For normal marginals, the transformation is evaluated directly as
#'
#' \deqn{
#' x = \mu + \sigma u.
#' }
#'
#' This is algebraically equivalent to
#'
#' \deqn{
#' x = F_X^{-1}(\Phi(u)),
#' }
#'
#' but remains numerically stable for extreme values of \code{u}.
#'
#' For all other distributions, the generic isoprobabilistic
#' transformation \code{q(pnorm(u))} is used. The probabilities are
#' clipped to avoid invalid inputs to the quantile function.
#'
#' @param u Numeric vector of values in standard normal space.
#' @param marginal Distribution funlist containing at least \code{$q}.
#'   For the direct normal transformation, \code{$DistributionType},
#'   \code{$mean}, and \code{$Sd} are additionally required.
#'
#' @return Numeric vector of values in physical space.
#'
#' @details
#' The direct normal-distribution branch preserves tail values, such as
#' \code{u = 40}, that would otherwise be lost when
#' \code{pnorm(u)} rounds to \code{1}.
#'
#' @keywords internal
#' @noRd
u_to_x <- function(u, marginal) {
  if (!is.numeric(u)) {
    stop("'u' must be numeric.", call. = FALSE)
  }

  if (!is.list(marginal) || !is.function(marginal$q)) {
    stop(
      "'marginal' must be a distribution funlist containing a quantile function '$q'.",
      call. = FALSE
    )
  }

  dist_type <- ""

  if (!is.null(marginal$DistributionType) &&
    length(marginal$DistributionType) == 1L &&
    !is.na(marginal$DistributionType)) {
    dist_type <- tolower(marginal$DistributionType)
  }

  # Direct and tail-stable transformation for normal marginals:
  #
  # qnorm(pnorm(u), mean = mu, sd = sigma) = mu + sigma * u
  #
  if (
    identical(dist_type, "norm") &&
      !is.null(marginal$mean) &&
      !is.null(marginal$Sd) &&
      length(marginal$mean) == 1L &&
      length(marginal$Sd) == 1L &&
      is.finite(marginal$mean) &&
      is.finite(marginal$Sd)
  ) {
    return(marginal$mean + marginal$Sd * u)
  }

  # Generic fallback for non-normal distributions.
  p <- stats::pnorm(u)
  p <- pmin(pmax(p, 1e-14), 1 - 1e-14)

  marginal$q(p)
}


#' @title Internal Monte Carlo Importance Sampling Routine (Single Limit State Function)
#'
#' @description
#' Performs the core Monte Carlo simulation with importance sampling for a single limit state function.
#' This internal routine is called by the public interface [`MC_IS()`] when the input argument `lsf`
#' is a single function rather than a list of system functions.
#'
#' The algorithm uses the FORM design point as center of the sampling density in standard normal space.
#' Each random variable is transformed between physical and normal spaces using its distribution object
#' (`PROB_BASEVAR$getlDistr()`), which provides parameterized CDF and inverse-CDF functions:
#'
#' \deqn{u = \Phi^{-1}(F_X(x))}
#' \deqn{x = F_X^{-1}(\Phi(u))}
#'
#' Parallel execution is supported via either base R's **parallel** package or the **future** backend.
#'
#' @param lsf Objective (limit state) function in form \code{function(x) \{ ... \}}.
#' @param lDistr List of input variable distributions created from TesiproV base variable objects.
#' @param cov_user Target coefficient of variation for convergence criterion.
#' @param n_batch Number of samples per batch (recommended > 1000).
#' @param n_max Maximum total number of samples before forced stop.
#' @param use_threads Number of threads/workers to use for parallel computation.
#' @param dataRecord Logical; if `TRUE`, all intermediate results are recorded and returned as data frame.
#' @param densityType Character; type of sampling density ("origin" keeps original input distributions).
#' @param dps Optional vector of design points supplied by user instead of FORM result.
#' @param debug.level Integer verbosity level: 0 = silent, 1 = info, 2 = detailed output during calculation.
#' @param streams RNG stream object created by [`init_rng_master()`].
#' @param libPaths_local Local library paths passed to cluster workers for package availability.
#' @param seed Optional integer seed for reproducibility across parallel runs. If `NULL`, a random seed is used.
#' @param backend Character; either `"parallel"` or `"future"` backend for parallelization on multiple cores.
#' @param adaptive_alpha Logical; if \code{TRUE}, adaptive mixture weights are used
#'   (system case only).
#'
#' @param alpha_update_rate Damping factor \eqn{\lambda \in [0,1]} controlling
#'   adaptation speed of mixture weights.
#'
#' @param adaptive_batch Logical; if TRUE, the batch size is adapted dynamically
#'   based on the current effective sample size (ESS) and target RSE.
#'
#' @param batch_control List controlling adaptive batch behaviour:
#'   \describe{
#'     \item{n_min}{Minimum batch size}
#'     \item{n_max}{Maximum batch size}
#'     \item{K_future}{Desired number of remaining iterations}
#'   }
#'
#' @param min_adapt_samples Integer. Minimum total number of MC samples required before adaptive mixtures are updated.
#'    The adaption mechanism is activated only after this threshold has been exceeded to ensure statistically stable estimates.
#'
#' @param alpha_min Numeric scalar in \eqn{(0,1)}.
#'   Lower bound for mixture weights \eqn{\alpha_i} during adaptive
#'   updating. This prevents individual components of the importance
#'   sampling density from collapsing to zero and ensures numerical
#'   stability in multimodal system analyses.
#'
#' @param stability_mode Character string specifying the numerical
#'   stabilization strategy used for importance sampling weight
#'   accumulation.
#'
#'   \describe{
#'     \item{"robust"}{
#'       Default. Uses full log-domain accumulation of weights
#'       across Monte Carlo iterations. This improves numerical
#'       stability for high-dimensional systems and very small
#'       failure probabilities.
#'     }
#'     \item{"fast"}{
#'       Uses shifted exponentiation with batch-level scaling.
#'       Computationally equivalent under normal engineering
#'       conditions but may be slightly less stable in extreme
#'       rare-event scenarios.
#'     }
#'   }
#'
#'   Both modes are mathematically equivalent for typical
#'   reliability problems. The \code{"robust"} mode is recommended
#'   for research applications and extreme reliability levels.
#'
#' @keywords internal
#'
MC_IS_single <- function(lsf, lDistr,
                         cov_user,
                         n_batch,
                         n_max,
                         use_threads,
                         backend = "future",
                         dataRecord,
                         densityType = "norm",
                         dps = NULL,
                         debug.level,
                         streams,
                         libPaths_local,
                         seed = NULL,
                         adaptive_alpha = FALSE,
                         alpha_update_rate = 0.1,
                         adaptive_batch = FALSE,
                         batch_control = list(),
                         min_adapt_samples = 10000,
                         alpha_min = 0.02,
                         stability_mode) {
  if (debug.level >= 1) {
    cat("DEBUG LEVEL INSIDE MC_IS_single:", debug.level, "\n")
  }

  ## --- Input validation ------------------------------------------------------
  if (missing(debug.level) || is.null(debug.level)) {
    debug.level <- 0
  }
  if (!is.function(lsf)) {
    stop("MC_IS_single: 'lsf' must be a function.", call. = FALSE)
  }
  if (!is.list(lDistr) || length(lDistr) < 1) {
    stop("MC_IS_single: 'lDistr' must be a list of marginals.", call. = FALSE)
  }

  if (!is.numeric(n_batch) || length(n_batch) != 1L || n_batch <= 0L) {
    stop("`n_batch` must be a positive integer (>= 1).")
  }

  if (!is.numeric(n_max) || length(n_max) != 1L || n_max <= 0L) {
    stop("`n_max` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(use_threads) || length(use_threads) != 1L || use_threads < 1L) {
    stop("`use_threads` must be an integer >= 1.", call. = FALSE)
  }

  debug.TAG <- "MC_IS_single"
  if (debug.level >= 1) message("[MC_IS_single] Starting single-LSF MC_IS...")

  ## -------------------------------------------------
  ## Deterministic RNG manager (thread-independent)
  ## -------------------------------------------------
  rng_mgr <- create_rng_manager(seed)

  ## -------------------------------------------------
  ## Cluster setup for the Windows parallel backend
  ##
  ## On Unix-like systems, run_parallel() uses mclapply().
  ## Therefore, no persistent cluster is required there.
  ## -------------------------------------------------

  cl <- NULL

  if (
    backend == "parallel" &&
      use_threads > 1L &&
      .Platform$OS.type == "windows"
  ) {
    cl <- make_parallel_cluster(
      use_threads = use_threads,
      seed = seed
    )

    on.exit(parallel::stopCluster(cl), add = TRUE)
  }

  # ## Master RNG / seed handling
  # RNGkind("L'Ecuyer-CMRG")

  # # Normalise seed for future.apply and cluster
  # seed_for_future <- NULL # integer or 7-length stream or NULL
  # streams_master <- NULL # 7-length stream (for parallel::clusterSetRNGStream)

  # if (!is.null(seed)) {
  #   if (is.numeric(seed) && length(seed) == 1L && !is.na(seed)) {
  #     # integer seed -> reproducible for future.apply
  #     set.seed(as.integer(seed))
  #     # derive a full CMRG stream for cluster backend if needed
  #     streams_master <- parallel::nextRNGStream(.Random.seed)
  #     seed_for_future <- as.integer(seed)
  #     if (debug.level >= 1) message(sprintf("[MC_IS_single] Fixed master seed set to %d", as.integer(seed)))
  #   } else if (is.numeric(seed) && length(seed) == 7L) {
  #     # user supplied full L'Ecuyer stream
  #     # assign(".Random.seed", seed, envir = .GlobalEnv)
  #     set.seed(seed)
  #     streams_master <- seed
  #     seed_for_future <- seed
  #     if (debug.level >= 1) message("[MC_IS_single] Using supplied L'Ecuyer-CMRG stream.")
  #   } else {
  #     stop("Seed must be a single integer or a L'Ecuyer-CMRG stream (length 7).", call. = FALSE)
  #   }
  # } else {
  #   # no fixed seed -> random behaviour; for future.apply use future.seed = TRUE (resolve_future_seed(NULL) -> TRUE)
  #   streams_master <- parallel::nextRNGStream(.Random.seed)
  #   seed_for_future <- NULL
  #   if (debug.level >= 1) message("[MC_IS_single] No fixed seed supplied; using random streams.")
  # }

  # Valid density types - origin means: the original input distribution
  # (from the package defined in PROB_BASEVAR) is used.


  allowed_density <- c(
    "norm", "lnorm", "weibull", "gamma", "beta",
    "exponential", "cauchy", "logis", "origin"
  )

  if (!densityType %in% allowed_density) {
    stop(
      sprintf(
        "`densityType = \"%s\"` is not supported. Supported types are: %s. ",
        densityType,
        paste(allowed_density[allowed_density != "origin"], collapse = ", ")
      ),
      "Use \"origin\" to keep the original input distribution."
    )
  }

  # check if lDistr is properly defined
  .check_lDistr(lDistr, "MC_IS")

  tic <- proc.time()

  n_vars <- length(lDistr)
  # -------------------------------------------------
  # Normalize marginal structure for FORM and MC
  # -------------------------------------------------

  if (!is.null(lDistr[[1]]$p) && is.function(lDistr[[1]]$p)) {
    # Funlist input (e.g. build_gsys)
    lDistr_fun <- lDistr

    lDistr_form <- lapply(lDistr_fun, function(fun) {
      params <- c(NA_real_, NA_real_)

      if (!is.null(fun$mean)) params[1] <- fun$mean
      if (!is.null(fun$Sd)) params[2] <- fun$Sd

      list(fun, params)
    })
  } else {
    # SYS machine format
    lDistr_form <- lDistr
    lDistr_fun <- lapply(lDistr_form, function(d) d[[1]])
  }


  ### cat("total_batch =", total_batch, "\n")

  if (n_vars < 1) stop("[MC_IS_single] lDistr must contain at least one marginal.", call. = FALSE)
  if (n_batch > n_max) n_batch <- n_max

  # ---------------------------------------------------------------------------
  # Design Point calculation using FORM
  # ---------------------------------------------------------------------------

  if (is.numeric(dps)) {
    dp <- as.numeric(dps)

    if (length(dp) == 1L && n_vars > 1L) {
      dp <- rep(dp, n_vars)
    }

    if (length(dp) != n_vars) {
      stop("[MC_IS_single] dps must have length 1 or n_vars", call. = FALSE)
    }

    res_form <- list(
      beta = NA_real_,
      pf = NA_real_,
      x_points = dp
    )

    if (debug.level >= 1) {
      message("[MC_IS_single] Using user supplied design point")
    }
  } else {
    # --- use EXACT same structure as SYS / FORM
    # FORM only provides the start point for importance sampling. Missing
    # convergence (e.g. noisy gradients of empirical distributions) is
    # harmless as long as the design point is finite (checked below).
    res_form <- suppressWarnings(TesiproV::FORM(
      lsf,
      lDistr_form,
      n_optim = 20,
      loctol = 1e-4
    ))

    if (is.null(res_form$x_points) || any(!is.finite(res_form$x_points))) {
      stop("[MC_IS_single] FORM failed to compute valid design point.", call. = FALSE)
    }

    if (debug.level >= 1) {
      message(sprintf(
        "[MC_IS_single] FORM beta = %.4f",
        res_form$beta
      ))
    }

    # ------------------------------------------------
    # Transform x* -> u*
    # ------------------------------------------------

    x_vec <- as.numeric(res_form$x_points)

    dp <- numeric(n_vars)

    for (i in seq_len(n_vars)) {
      Fxi <- lDistr_fun[[i]]$p(x_vec[i])

      Fxi <- min(max(Fxi, 1e-14), 1 - 1e-14)

      dp[i] <- qnorm(Fxi)
    }

    if (debug.level >= 2) {
      message(sprintf(
        "[MC_IS_single] Design point in u-space: %s",
        paste(round(dp, 4), collapse = ", ")
      ))
    }

    # ------------------------------------------------------------
    # Stabilized design-point shift
    # ------------------------------------------------------------

    ### To do: implement switch use_shift_damping, but default = FALSE

    ### shift_damping <- min(1, 5 / sqrt(length(dp)))

    ### dp <- shift_damping * dp

    dp <- dp
  }

  alpha_vec <- dp / sqrt(sum(dp^2)) # unit vector in direction of dp
  beta_dp <- if (is.finite(res_form$beta)) res_form$beta else sqrt(sum(dp^2))
  c_scale <- beta_dp^2
  Sigma_chol <- chol(diag(n_vars) + (c_scale - 1) * outer(alpha_vec, alpha_vec))
  Sigma_chol_inv <- solve(Sigma_chol)
  log_det_Sigma_chol <- sum(log(diag(Sigma_chol)))

  # Sigma_chol <- diag(n_vars)
  # Sigma_chol_inv <- diag(n_vars)
  # log_det_Sigma_chol <- 0

  # ---------------------------------------------------------------------------
  # Detect whether LSF supports vectorized matrix input
  # ---------------------------------------------------------------------------

  lsf_vectorized <- FALSE

  if (debug.level >= 1) {
    message("[MC_IS_single] Testing if LSF supports vectorized evaluation...")
  }

  # Create small test matrix (2 samples)
  # x-space design point
  x_star <- as.numeric(unlist(res_form$x_points))

  test_mat <- matrix(
    rep(x_star, 2),
    nrow = 2,
    byrow = TRUE
  )
  # test_mat <- matrix(dp, nrow = 2, ncol = n_vars, byrow = TRUE)

  test_result <- try(lsf(test_mat), silent = TRUE)

  if (!inherits(test_result, "try-error") &&
    is.numeric(test_result) &&
    length(test_result) == 2 &&
    all(is.finite(test_result))) {
    lsf_vectorized <- TRUE

    if (debug.level >= 1) {
      message("[MC_IS_single] LSF detected as vectorized.")
    }
  } else {
    if (debug.level >= 1) {
      message("[MC_IS_single] LSF is NOT vectorized - falling back to row-wise evaluation.")
    }
  }


  # ---------------------------------------------------------------------------
  # Local Monte Carlo routine executed on each worker/thread
  # ---------------------------------------------------------------------------
  # Worker-safe copy: do not capture MC_IS_single's local environment.
  u_to_x_worker <- u_to_x
  environment(u_to_x_worker) <- baseenv()

  mc_local <- function(u_mat) {
    n_local <- nrow(u_mat)

    # --------------------------------------------------
    # Transform u -> x (vectorized)
    # --------------------------------------------------

    # code is slow
    # p_mat <- pnorm(u_mat)
    # p_mat <- pmin(pmax(p_mat, 1e-14), 1 - 1e-14)

    # q_list <- lapply(lDistr, function(d) d$q)

    # x_mat <- matrix(NA_real_, n_local, n_vars)

    # for (j in seq_len(n_vars)) {
    #   x_mat[, j] <- q_list[[j]](p_mat[, j])
    # }

    # Fix

    x_cols <- vector("list", n_vars)

    for (j in seq_len(n_vars)) {
      x_cols[[j]] <- u_to_x(
        u = u_mat[, j],
        marginal = lDistr_fun[[j]]
      )
    }

    x_mat <- do.call(cbind, x_cols)

    # --------------------------------------------------
    # log phi(u)  (standard normal)
    # --------------------------------------------------

    log_phi <- -0.5 * rowSums(u_mat^2) -
      0.5 * n_vars * log(2 * pi)

    diff <- sweep(u_mat, 2, dp, "-")

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
        min(u_mat),
        max(u_mat)
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
    # cat("mean(u_mat) =", colMeans(u_mat), "\n")

    # list(
    #   log_w  = log_w,
    #   I_vals = I_vals
    # )


    # --------------------------------------------------
    # Stable aggregation on the log scale
    # --------------------------------------------------

    log_sum_w <- log_sum_exp(log_w)

    if (any(I_vals == 1L)) {
      log_sum_Iw <- log_sum_exp(log_w[I_vals == 1L])
    } else {
      log_sum_Iw <- -Inf
    }

    log_sum_w2 <- log_sum_exp(2 * log_w)

    # second moment of I * w for the variance of the estimator
    if (any(I_vals == 1L)) {
      log_sum_Iw2 <- log_sum_exp(2 * log_w[I_vals == 1L])
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


  # ## ---- Cluster setup --------------------------------------------------------
  # cl <- NULL # placeholder in case no cluster is needed

  # ## ---- 1a. Determine seed for the cluster -------------------------------
  # ## `streams` is the vector returned by `init_rng_master()`.
  # ## The first element is a valid integer seed.
  # ## seed_used <- as.integer(streams[1])

  # # streams <- init_rng_master(NULL)    # NULL -> new random seed per call
  # # seed_used <- as.integer(streams[1])

  # ## -------------------------------------------------------------------------
  # ## 2) Choice of parallel backend
  # ## -------------------------------------------------------------------------
  # if (backend == "parallel" && Sys.info()[["sysname"]] != "Windows") {
  #   ## --- Prepare RNG stream --------------------------------------------
  #   ## The parallel backend needs the complete CMRG stream
  #   if (is.null(streams_master)) streams_master <- parallel::nextRNGStream(.Random.seed)

  #   cl <- make_ready_cluster(
  #     use_threads = use_threads,
  #     libPaths_local = libPaths_local,
  #     seed = streams_master,
  #     export_objs = c(
  #       "libPaths_local", "mc_local", "lsf", "lDistr", "dp",
  #       "n_batch", "densityType", "n_vars"
  #     ),
  #     debug.level = debug.level
  #   )

  #   on.exit(parallel::stopCluster(cl), add = TRUE)
  # } else if (backend == "future") {
  #   if (!requireNamespace("future", quietly = TRUE)) {
  #     stop("Package 'future' needed for backend = 'future'", call. = FALSE)
  #   }

  #   if (!requireNamespace("future.apply", quietly = TRUE)) {
  #     stop("Package 'future.apply' needed for backend = 'future'", call. = FALSE)
  #   }

  #   # ------------------------------------------------------------
  #   # CRAN-conform future handling
  #   # ------------------------------------------------------------
  #   # Save current plan
  #   old_plan <- future::plan()

  #   # Restore original plan when function exits
  #   on.exit(
  #     {
  #       future::plan(old_plan)
  #     },
  #     add = TRUE
  #   )

  #   # Only set a parallel plan if the user did not already set one
  #   if (inherits(old_plan, "sequential") && use_threads > 1) {
  #     if (.Platform$OS.type == "windows") {
  #       future::plan(future::multisession,
  #         workers = use_threads
  #       )
  #     } else {
  #       future::plan(future::multicore,
  #         workers = use_threads
  #       )
  #     }

  #     if (debug.level >= 1) {
  #       message(sprintf(
  #         "[MC_IS_single] Temporary future plan set with %d workers.",
  #         use_threads
  #       ))
  #     }
  #   } else {
  #     if (debug.level >= 1) {
  #       message("[MC_IS_single] Using existing future plan (no override).")
  #     }
  #   }
  # } else {
  #   # stop("Unsupported backend")
  #   if (backend != "parallel" && backend != "future") stop("Unsupported backend", call. = FALSE)
  # }

  ## -------------------------------------------------------------------------
  ## 3. Monte Carlo loop (identical for both backends)
  ## -------------------------------------------------------------------------
  n_sim <- 0L
  pf_mc <- 0
  var_mc <- 0
  w_mc <- 0
  I_mc <- 0
  var <- NA_real_
  cov <- 1
  pf <- NA_real_

  # -------------------------------------------------
  # Stable global statistics on the log scale
  # -------------------------------------------------
  log_total_sum_Iw <- -Inf
  log_total_sum_w <- -Inf
  log_total_sum_w2 <- -Inf
  log_total_sum_Iw2 <- -Inf

  total_n <- 0L


  if (dataRecord) {
    # rec <- make_recorders(ceiling(n_max / (use_threads * n_batch)))
    rec <- make_recorders(ceiling(n_max / (n_batch)))
  }

  ### cat(
  ###   "Initial RNG head:",
  ###   paste(.Random.seed[1:5], collapse = ", "),
  ###   "\n"
  ### )

  iter <- 0L

  while (TRUE) {
    iter <- iter + 1L

    # Limit each batch to the samples still allowed by n_max.
    remaining <- floor(n_max - total_n)

    if (remaining <= 0) {
      break
    }

    total_batch <- as.integer(min(floor(n_batch), remaining))

    # ============================================================
    # 3a. Deterministic sampling block (RNG only in master)
    # ============================================================

    # 1) Generate ALL samples in one deterministic block
    # u_mat <- matrix(rnorm(total_batch * n_vars),
    #   ncol = n_vars
    # )
    if (debug.level >= 1) {
      t_sample_start <- proc.time()
    }
    ## -------------------------------------------------
    ## Deterministic sampling block (thread independent)
    ## -------------------------------------------------

    # u_mat <- matrix(u_vec, ncol = n_vars)

    z_mat <- matrix(rng_mgr$generate_norm_block(total_batch * n_vars), ncol = n_vars)
    u_mat <- sweep(z_mat %*% Sigma_chol, 2, dp, "+")

    if (debug.level >= 1) {
      t_sample <- (proc.time() - t_sample_start)[3]
    }

    ### cat(
    ###   "After sampling RNG head:",
    ###   paste(.Random.seed[1:5], collapse = ", "),
    ###   "\n"
    ### )

    # u_mat <- sweep(u_mat, 2, dp, "+")

    # 2) Deterministic split into worker chunks
    # if (use_threads <= 1) {
    #   chunk_indices <- list(seq_len(total_batch))
    # } else {
    #   chunk_indices <- split(
    #     seq_len(total_batch),
    #     cut(seq_len(total_batch),
    #       breaks = use_threads,
    #       labels = FALSE
    #     )
    #   )
    # }

    if (use_threads == 1) {
      chunk_indices <- list(seq_len(total_batch))
    } else {
      chunk_indices <- split(
        seq_len(total_batch),
        cut(seq_len(total_batch), breaks = use_threads)
      )
    }

    ### cat(
    #   "Chunk sizes:",
    #   sapply(chunk_indices, length),
    #   "\n"
    # )

    # 3) Worker function (NO RNG inside!)
    worker_fun <- function(idx) {
      mc_local(u_mat[idx, , drop = FALSE])
    }

    # # 4) Parallel evaluation (RNG disabled!)
    # if (use_threads == 1) {
    #   res_list <- lapply(chunk_indices, worker_fun)
    # } else if (backend == "parallel") {
    #   if (.Platform$OS.type == "windows") {
    #     res_list <- lapply(chunk_indices, worker_fun)
    #   } else {
    #     res_list <- parallel::parLapply(
    #       cl,
    #       chunk_indices,
    #       worker_fun
    #     )
    #   }
    # } else if (backend == "future") {
    #   res_list <- future.apply::future_lapply(
    #     chunk_indices,
    #     worker_fun,
    #     future.seed = FALSE,
    #     future.globals = FALSE # CRITICAL
    #   )
    # }
    if (debug.level >= 1) {
      t_parallel_start <- proc.time()
    }
    ## -------------------------------------------------
    ## Unified parallel evaluation (no RNG in workers)
    ## -------------------------------------------------
    # res_list <- run_parallel(
    #   X = seq_len(use_threads),
    #   FUN = function(worker_id) {
    #     idx <- chunk_indices[[worker_id]]
    #     worker_fun(idx)
    #   },
    #   backend = backend,
    #   use_threads = use_threads,
    #   seed = NULL # IMPORTANT: workers must not use RNG
    # )
    # res_list <- run_parallel(
    #   X = chunk_indices,
    #   FUN = function(idx) {
    #     mc_local(u_mat[idx, , drop = FALSE])
    #   },
    #   backend = backend,
    #   use_threads = use_threads,
    #   cl = cl
    # )

    u_chunks <- lapply(
      chunk_indices,
      function(idx) u_mat[idx, , drop = FALSE]
    )

    res_list <- run_parallel(
      X = u_chunks,
      FUN = .mc_is_single_worker,
      backend = backend,
      use_threads = use_threads,
      cl = cl,
      lDistr_fun = lDistr_fun,
      u_to_x_worker = u_to_x_worker,
      dp = dp,
      Sigma_chol_inv = Sigma_chol_inv,
      log_det_Sigma_chol = log_det_Sigma_chol,
      n_vars = n_vars,
      lsf = lsf,
      lsf_vectorized = lsf_vectorized,
      debug.level = debug.level
    )

    if (debug.level >= 1) {
      t_parallel <- (proc.time() - t_parallel_start)[3]
    }

    if (debug.level >= 1) {
      cat(sprintf(
        "Sampling: %.3f s | Parallel: %.3f s\n",
        t_sample, t_parallel
      ))
    }


    # -------------------------------------------------
    # 3b. Stable aggregation across worker chunks
    # -------------------------------------------------

    log_sum_Iw_iter <- -Inf
    log_sum_w_iter <- -Inf
    log_sum_w2_iter <- -Inf
    log_sum_Iw2_iter <- -Inf
    n_iter <- 0L

    for (worker_res in res_list) {
      log_sum_Iw_iter <- log_add_exp(
        log_sum_Iw_iter,
        worker_res$log_sum_Iw
      )

      log_sum_w_iter <- log_add_exp(
        log_sum_w_iter,
        worker_res$log_sum_w
      )

      log_sum_w2_iter <- log_add_exp(
        log_sum_w2_iter,
        worker_res$log_sum_w2
      )

      log_sum_Iw2_iter <- log_add_exp(
        log_sum_Iw2_iter,
        worker_res$log_sum_Iw2
      )

      n_iter <- n_iter + worker_res$n
    }

    # -------------------------------------------------
    # Stable accumulation across MC iterations
    # -------------------------------------------------

    log_total_sum_Iw <- log_add_exp(
      log_total_sum_Iw,
      log_sum_Iw_iter
    )

    log_total_sum_w <- log_add_exp(
      log_total_sum_w,
      log_sum_w_iter
    )

    log_total_sum_w2 <- log_add_exp(
      log_total_sum_w2,
      log_sum_w2_iter
    )

    log_total_sum_Iw2 <- log_add_exp(
      log_total_sum_Iw2,
      log_sum_Iw2_iter
    )

    total_n <- total_n + n_iter


    # --------------------------------------------------------
    # FAST MODE
    # --------------------------------------------------------
    # if (stability_mode == "fast") {
    #   max_log_w <- max(all_log_w)
    #   w_all <- exp(all_log_w - max_log_w)

    #   # accumulate using base sum (fast, stable enough here)
    #   sum_w_iter <- sum(w_all)
    #   sum_Iw_iter <- sum(w_all[all_I_vals == 1])

    #   # second moment for ESS
    #   sum_w2_iter <- sum(w_all * w_all)
    # }

    # # --------------------------------------------------------
    # # ROBUST MODE (numerically stable weight accumulation)
    # # --------------------------------------------------------
    # if (stability_mode == "robust") {
    #   # Compute global shift (log-sum-exp reference)
    #   # This prevents overflow when exponentiating large log-weights
    #   max_log_w <- max(all_log_w)

    #   # Shift log-weights before exponentiation
    #   # This is mathematically equivalent to fast-mode scaling
    #   shifted_weights <- exp(all_log_w - max_log_w)

    #   # First moment: sum of weights
    #   sum_w_iter <- sum(shifted_weights)

    #   # Weighted failure indicator
    #   if (any(all_I_vals == 1, na.rm = TRUE)) {
    #     sum_Iw_iter <- sum(shifted_weights[all_I_vals == 1])
    #   } else {
    #     sum_Iw_iter <- 0
    #   }

    #   # Second moment: required for ESS computation
    #   sum_w2_iter <- sum(shifted_weights * shifted_weights)
    # }

    # --- 3c. Recording ------------------------------------------------------------
    # record totals for this outer iteration

    n_sim <- total_n


    # -------------------------------------------------
    # 3d. Estimate failure probability and variance
    #     from log-domain sums
    # -------------------------------------------------

    log_pf <- -Inf
    pf <- 0
    log_var <- NA_real_
    var_est <- NA_real_
    log_cov <- Inf
    cov <- Inf
    ESS <- 0

    if (total_n > 0L) {
      est <- .mc_is_log_estimates(
        log_sum_Iw   = log_total_sum_Iw,
        log_sum_Iw2  = log_total_sum_Iw2,
        n            = total_n
      )

      log_pf <- est$log_pf
      pf <- est$pf
      log_var <- est$log_var
      var_est <- est$var
      log_cov <- est$log_cov
      cov <- est$cov

      # Preserve MC_IS_single's existing ESS definition:
      # ESS = (sum of all weights)^2 / sum of squared all weights.
      # It is calculated only when failure-weight mass was observed,
      # matching the previous behavior.
      if (is.finite(log_total_sum_Iw) &&
        is.finite(log_total_sum_w) &&
        is.finite(log_total_sum_w2)) {
        log_ESS <- 2 * log_total_sum_w - log_total_sum_w2
        ESS <- min(exp(log_ESS), total_n)
      }
    }


    # -------------------------------------------------
    # 3e. Record the values calculated in this iteration
    # -------------------------------------------------

    if (dataRecord) {
      rec <- record_step(
        rec = rec,
        idx = iter,
        pf_mc = pf,
        var_mc = var_est,
        I_mc = NA_real_,
        n_sim = n_sim,
        tic = tic
      )

      rec$cov[iter] <- cov
    }

    # if (total_sum_w > 0 && total_n > 0) {
    #
    #   pf_hat <- total_sum_Iw / total_sum_w
    #
    #   # Effective Sample Size
    #   ESS <- (total_sum_w^2) / total_sum_w2

    # if (pf_hat > 0 && ESS > 1) {
    #   # Correct self-normalized IS variance
    #   # Var[pf_hat] ~= (1/N) * (E[I^2*w^2]/E[w]^2 - pf_hat^2)
    #   # simplified with ESS:
    #   #var_est <- pf_hat * (1 - pf_hat) / ESS
    #   var_est <- (1 / total_n) *
    #     (total_sum_Iw2 / total_sum_w^2 - pf_hat^2)
    #   cov     <- sqrt(var_est) / pf_hat
    # } else if (pf_hat == 0) {
    #   var_est <- NA_real_
    #   cov     <- Inf
    # } else {
    #   var_est <- NA_real_
    #   cov     <- Inf
    # }

    #   if (pf_hat > 0 && ESS > 1) {
    #
    #     inner <- total_sum_Iw2 / total_sum_w^2 - pf_hat^2
    #
    #     if (is.finite(inner) && inner > 0) {
    #       var_est <- inner / total_n
    #       cov     <- sqrt(var_est) / pf_hat
    #     } else {
    #       var_est <- NA_real_
    #       cov     <- Inf
    #     }
    #
    #   } else if (pf_hat == 0 || total_sum_Iw == 0) {
    #     var_est <- NA_real_
    #     cov     <- Inf
    #   } else {
    #     var_est <- NA_real_
    #     cov     <- Inf
    #   }
    #
    #   pf <- pf_hat
    #
    # } else {
    #   pf      <- 0
    #   var_est <- NA_real_
    #   cov     <- Inf
    # }

    # if (total_sum_w > 0 && total_n > 0) {
    #
    #   pf_hat <- total_sum_Iw / total_sum_w
    #
    #   ESS <- (total_sum_w^2) / total_sum_w2
    #
    #   if (pf_hat > 0 && ESS > 1) {
    #     var_est <- pf_hat * (1 - pf_hat) / ESS
    #     cov     <- sqrt(var_est) / pf_hat
    #   } else if (pf_hat == 0) {
    #     var_est <- NA_real_
    #     cov     <- Inf
    #   } else {
    #     var_est <- NA_real_
    #     cov     <- Inf
    #   }
    #
    #   pf <- pf_hat
    #
    # } else {
    #   pf      <- 0
    #   var_est <- NA_real_
    #   cov     <- Inf
    #   ESS     <- 0
    # }
    # if (total_sum_w > 0 && total_n > 0) {
    #   pf_hat <- total_sum_Iw / total_sum_w

    #   S_w2 <- total_sum_w2
    #   S_w <- total_sum_w
    #   N <- total_n

    #   if (total_sum_w > 0) {
    #     pf_hat <- total_sum_Iw / total_sum_w

    #     # Effective Sample Size
    #     ESS <- (total_sum_w^2) / total_sum_w2

    #     # Stable ESS-based variance approximation
    #     if (ESS > 0 && pf_hat > 0) {
    #       var_est <- pf_hat * (1 - pf_hat) / ESS
    #       cov <- sqrt(var_est) / pf_hat
    #     } else {
    #       var_est <- NA_real_
    #       cov <- Inf
    #     }
    #   } else {
    #     pf_hat <- 0
    #     var_est <- NA_real_
    #     cov <- Inf
    #     ESS <- 0
    #   }

    #   pf <- pf_hat
    # } else {
    #   # no weight mass -> pf = 0
    #   pf <- 0
    #   var_est <- 0
    #   cov <- Inf
    # }

    # debug printing (optional)
    if (debug.level >= 1) {
      message(sprintf("[MC_IS_single] iter=%d, n_sim=%d, pf=%.12g, cov=%.4f, E=%.4f", iter, n_sim, pf, cov, exp(log_total_sum_w - log(total_n))))
    }
    # convergence checks
    if (length(cov_user) == 1L &&
      is.finite(cov_user) &&
      cov_user > 0 &&
      !is.na(log_cov) &&
      log_cov < log(cov_user)) {
      break
    }

    if (n_sim >= n_max) {
      break # maximum iterations reached
    }
  }


  ## -------------------------------------------------------------------------
  ## 4. Assemble results
  ## -------------------------------------------------------------------------
  duration <- proc.time() - tic

  df <- if (dataRecord) as_record_df(rec, iter) else data.frame()


  if (!is.numeric(pf) || is.na(pf)) pf <- NA_real_

  beta_val <- if (is.finite(log_pf)) {
    -stats::qnorm(log_pf, log.p = TRUE)
  } else {
    NA_real_
  }

  output <- list(
    method = "MCIS_Single",
    beta = beta_val,
    pf = pf,
    log_pf = log_pf,
    log_var = log_var,
    FORM_beta = res_form$beta,
    FORM_pf = res_form$pf,
    design_points = dp,
    var = var_est,
    cov_mc = cov,
    cov_user = cov_user,
    ESS = ESS,
    ESS_ratio = if (n_sim > 0) ESS / n_sim else NA_real_,
    n_mc = n_sim,
    n_max = n_max,
    data = df,
    runtime = duration[1:5]
  )

  if (debug.level >= 1) message(sprintf("[MC_IS_single] Monte-Carlo finished in [s]: %g", duration[1]))
  return(output)
}
