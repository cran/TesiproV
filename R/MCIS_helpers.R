#' @title Internal helper functions for Monte Carlo Importance Sampling (MC_IS)
#' @description
#' This file contains internal helper functions used by the main function
#' [`MC_IS()`] to perform Monte Carlo simulations with importance sampling.
#'
#' The helpers defined here support:
#' * random number generation and cluster setup,
#' * creation of recorder objects for iterative logging,
#' * single limit state simulation routine (`MC_IS_single()`),
#' * conversion of recorded data to data frames after completion.
#'
#' Additional helper sets will be provided for system calculations
#' (serial and parallel systems) in separate files.
#'
#' @details
#' These functions are not intended to be called directly by users.
#' They are exported internally to allow modular access within the package.
#'
#' @seealso [TesiproV::MC_IS()] for the public interface that calls these helpers.
#'
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#' @name Internal helper functions for MC_IS()
NULL


## ---------------------------------------------------------------
## Random number generator utilities
## ---------------------------------------------------------------

#' @title Initialize master RNG stream
#' @description Sets RNG type `"L'Ecuyer-CMRG"` and optionally a fixed seed.
#'
#' @param seed Optional integer seed; if `NULL`, a random base seed is generated once per call.
#' @param debug.level Integer controlling verbosity level.
#'   Values >= 1 produce diagnostic console output.
#' @return A master RNG stream object (numeric vector of length 7)
#'         as required by [parallel::clusterSetRNGStream()] or usable
#'         directly in [future.apply::future_lapply()].
#'
## @keywords internal
init_rng_master <- function(seed = NULL, debug.level = 0) {
  RNGkind("L'Ecuyer-CMRG")

  if (!is.null(seed)) {
    # fixed user seed -> reproducible
    set.seed(as.integer(seed))
    if (debug.level >= 0L) message(sprintf("[init_rng_master] User seed set to %d", seed))
  } else {
    base_seed <- sample.int(.Machine$integer.max, 1L)
    set.seed(base_seed)
    if (debug.level >= 0L) message(sprintf("[init_rng_master] Random RNG seed drawn: %d", base_seed))
  }

  # produce full L'Ecuyer-CMRG stream and return it
  stream <- parallel::nextRNGStream(.Random.seed) # Vector length 7
  if (debug.level >= 1L) {
    message(sprintf("[init_rng_master] L'Ecuyer stream head: %d (use this for cluster RNG)", as.integer(stream[1])))
  }
  return(stream)
}

#' @title Set cluster RNG streams
#' @description Assigns identical RNG streams to each worker in an existing cluster.
#'
#' @param cl Cluster object created by `parallel::makeCluster()`.
#' @param seed Integer base seed used on the master process.
#'
#' @keywords internal
set_cluster_rng <- function(cl, seed) {
  ## Es muss kein neuer Seed erzeugt werden - wir benutzen denselben
  ## Basis-Seed, den `init_rng_master()` bereits gesetzt hat.
  parallel::clusterSetRNGStream(cl, iseed = seed)
}


## ---------------------------------------------------------------
## Recording / logging utilities for iterative MC loops
## ---------------------------------------------------------------

#' @title Create empty recorder vectors for MC iterations
#' @description Preallocates numeric vectors and time tracking data frames used to record intermediate results during simulation loops.
#'
#' @param n_iter Expected maximum number of iterations/batches in the simulation loop.
#'
#' @return A list containing numeric vectors (`pf`, `var`, `cov`, etc.) and a nested time data frame with user/system/elapsed times per iteration.
#'
#' @keywords internal
make_recorders <- function(n_iter) {
  list(
    pf = numeric(n_iter),
    var = numeric(n_iter),
    cov = numeric(n_iter),
    I_sum = numeric(n_iter),
    n_sim = numeric(n_iter),
    time = data.frame(
      user = numeric(n_iter),
      sys = numeric(n_iter),
      elapsed = numeric(n_iter),
      child_user = numeric(n_iter),
      child_sys = numeric(n_iter)
    )
  )
}


#' @title Internal Parallel Dispatcher for Importance Sampling
#'
#' @description
#' Internal helper function used by \code{MC_IS_single()} to distribute
#' Monte Carlo sample blocks across worker threads in a reproducible and
#' backend-independent manner.
#'
#' The function ensures worker-invariant execution by:
#' \itemize{
#'   \item Performing deterministic chunking of sample indices,
#'   \item Avoiding per-worker random number generation,
#'   \item Supporting both \code{"future"} and \code{"parallel"} backends.
#' }
#'
#' This abstraction isolates parallel execution from the statistical
#' core of the Importance Sampling algorithm.
#'
#' @param chunk_indices List of integer vectors defining sample indices
#'   assigned to each worker.
#'
#' @param worker_fun Function that evaluates a given subset of samples.
#'
#' @param u_mat Matrix of samples in standard normal space.
#'
#' @param use_threads Number of worker threads.
#'
#' @param backend Character string specifying backend
#'   (\code{"future"} or \code{"parallel"}).
#'
#' @param debug.level Verbosity level for diagnostic output.
#'
#' @return A list containing worker evaluation results.
#'
#' @keywords internal
parallel_dispatch <- function(chunk_indices,
                              worker_fun,
                              u_mat,
                              use_threads,
                              backend,
                              debug.level = 0) {
  # --------------------------------------------------
  # Single-thread fallback
  # --------------------------------------------------
  if (use_threads == 1) {
    if (debug.level >= 2) {
      message("[MC_IS_system] Using lapply (single thread).")
    }

    return(lapply(chunk_indices, worker_fun, u_mat))
  }

  # --------------------------------------------------
  # Unix-like systems (Linux/macOS)
  # --------------------------------------------------
  if (.Platform$OS.type != "windows") {
    if (debug.level >= 2) {
      message("[MC_IS_system] Using mclapply (fork).")
    }

    return(
      parallel::mclapply(
        chunk_indices,
        worker_fun,
        u_mat,
        mc.cores = use_threads
      )
    )
  }

  # --------------------------------------------------
  # Windows fallback -> future
  # --------------------------------------------------
  if (debug.level >= 2) {
    message("[MC_IS_system] Using future_lapply (Windows fallback).")
  }

  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Package 'future.apply' required for Windows parallel execution.")
  }

  return(
    future.apply::future_lapply(
      chunk_indices,
      worker_fun,
      u_mat,
      future.seed = FALSE,
      future.globals = FALSE
    )
  )
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
  ## --- Input validation ------------------------------------------------------
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

  ## Master RNG / seed handling
  RNGkind("L'Ecuyer-CMRG")

  # Normalise seed for future.apply and cluster
  seed_for_future <- NULL # integer or 7-length stream or NULL
  streams_master <- NULL # 7-length stream (for parallel::clusterSetRNGStream)

  if (!is.null(seed)) {
    if (is.numeric(seed) && length(seed) == 1L && !is.na(seed)) {
      # integer seed -> reproducible for future.apply
      set.seed(as.integer(seed))
      # derive a full CMRG stream for cluster backend if needed
      streams_master <- parallel::nextRNGStream(.Random.seed)
      seed_for_future <- as.integer(seed)
      if (debug.level >= 1) message(sprintf("[MC_IS_single] Fixed master seed set to %d", as.integer(seed)))
    } else if (is.numeric(seed) && length(seed) == 7L) {
      # user supplied full L'Ecuyer stream
      # assign(".Random.seed", seed, envir = .GlobalEnv)
      set.seed(seed)
      streams_master <- seed
      seed_for_future <- seed
      if (debug.level >= 1) message("[MC_IS_single] Using supplied L'Ecuyer-CMRG stream.")
    } else {
      stop("Seed must be a single integer or a L'Ecuyer-CMRG stream (length 7).", call. = FALSE)
    }
  } else {
    # no fixed seed -> random behaviour; for future.apply use future.seed = TRUE (resolve_future_seed(NULL) -> TRUE)
    streams_master <- parallel::nextRNGStream(.Random.seed)
    seed_for_future <- NULL
    if (debug.level >= 1) message("[MC_IS_single] No fixed seed supplied; using random streams.")
  }

  ## Gueltige density-Typen - origin bedeutet: die originale Eingabeverteilung
  ## (aus dem Paket, das im PROB_BASEVAR definiert ist) wird verwendet.
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


  tic <- proc.time()

  n_vars <- length(lDistr)
  total_batch <- n_batch # use_threads * n_batch

  ### cat("total_batch =", total_batch, "\n")

  if (n_vars < 1) stop("[MC_IS_single] lDistr must contain at least one marginal.", call. = FALSE)
  if (n_batch > n_max) n_batch <- n_max

  # ---------------------------------------------------------------------------
  # Design Point calculation using FORM
  # ---------------------------------------------------------------------------

  if (is.numeric(dps)) {
    dp <- as.numeric(dps)
    if (length(dp) == 1L && n_vars > 1L) dp <- rep(dp, n_vars)
    if (length(dp) != n_vars) stop("[MC_IS_single] dps must have length 1 or n_vars", call. = FALSE)
    res_form <- list(beta = NA_real_, pf = NA_real_, x_points = dp)
    if (debug.level >= 1) message("[MC_IS_single] Using user-supplied design points")
  } else {
    # If FORM expects the wrapper (funlist + params), build it
    # assume lDistr_funlist is list of funlists (with $p,$q,...)
    # --- Ensure we have funlists (one per marginal) -------------------------
    # lDistr_funlist <- result of normalization
    lDistr_funlist <- lapply(lDistr, function(dd) {
      # dd kann funlist oder wrapper(list(funlist, params)) sein
      if (is.list(dd) && !is.null(dd$p) && is.function(dd$p)) {
        return(dd)
      } # already funlist
      if (is.list(dd) && length(dd) >= 1 && is.list(dd[[1]]) && !is.null(dd[[1]]$p)) {
        return(dd[[1]])
      } # wrapper -> funlist
      # fallback: search for element that is a funlist
      if (is.list(dd)) {
        for (el in dd) {
          if (is.list(el) && !is.null(el$p) && is.function(el$p)) {
            return(el)
          }
        }
      }
      stop("[MC_IS_single] Cannot extract funlist from lDistr element", call. = FALSE)
    })

    # --- Build wrapper list expected by FORM: list(funlist, params) ----------
    # lDistr_for_FORM <- wrapper for FORM
    lDistr_for_FORM <- lapply(lDistr_funlist, function(fun) {
      params <- c(NA_real_, NA_real_)
      if (!is.null(fun$mean)) params[1] <- as.numeric(fun$mean)
      if (!is.null(fun$Sd)) params[2] <- as.numeric(fun$Sd)
      list(fun, params)
    })

    # Now call FORM with wrapper:
    res_form <- TesiproV::FORM(lsf, lDistr_for_FORM, n_optim = 20, loctol = 1e-4)
    if (is.null(res_form$x_points) || length(res_form$x_points) == 0L) {
      stop("[MC_IS_single] FORM did not return valid design point x_points.", call. = FALSE)
    }
    # res_form <- TesiproV::FORM(lsf, lDistr, n_optim = 20, loctol = 1e-4)
    # after building lDistr_funlist and lDistr_for_FORM and calling FORM
    # use funlist for MC:
    lDistr <- lDistr_funlist

    # Debug-Ausgabe zur Kontrolle der Struktur
    if (debug.level >= 1) {
      message("[MC_IS_single] Structure of FORM result:")
      str(res_form)
      message(sprintf(
        "Designpoint (x-space): %s",
        paste(round(as.numeric(unlist(res_form$x_points)), 3), collapse = ", ")
      ))
      message(sprintf("[MC_IS_single] Beta from FORM: %.3f", res_form$beta))
    }

    # transform FORM x* to u*
    x_vec <- as.numeric(unlist(res_form$x_points))
    if (length(x_vec) == 0L || any(is.na(x_vec))) {
      stop("[MC_IS_single] FORM returned no valid design point (x_points).", call. = FALSE)
    }
    # ensure length
    if (length(x_vec) < n_vars) {
      warning("FORM returned fewer x* components than variables; padding with first value.")
      x_vec <- rep(x_vec[1], n_vars)
    }
    dp <- numeric(n_vars)
    for (i in seq_len(n_vars)) {
      x_star_i <- as.numeric(x_vec[i])
      Fxi <- tryCatch(lDistr[[i]]$p(x_star_i), error = function(e) NA_real_)
      Fxi <- min(max(as.numeric(Fxi), 1e-14), 1 - 1e-14)
      dp[i] <- qnorm(Fxi)
      if (!is.function(lDistr[[i]]$p) || !is.function(lDistr[[i]]$q)) stop("[MC_IS_single] marginal missing p/q", call. = FALSE)
    }
    # check marginals or missing p or q
  }

  if (!is.numeric(dp) || length(dp) != n_vars) {
    stop("[MC_IS_single] 'dp' must be numeric vector of length n_vars", call. = FALSE)
  }
  if (debug.level >= 2) {
    message(sprintf("[MC_IS_single] Using design point (u-space): %s", paste(round(dp, 4), collapse = ", ")))
  }

  ## ---- Debug-Ausgabe des Designpunkts --------------------------------------
  if (debug.level >= 2) {
    message(sprintf(
      "[MC_IS_single] Using design point: %s",
      paste(round(dp, 3), collapse = ", ")
    ))
  }

  # ------------------------------------------------------------
  # Stabilized design-point shift
  # ------------------------------------------------------------

  shift_damping <- min(1, 5 / sqrt(length(dp)))

  dp <- shift_damping * dp

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
  mc_local <- function(u_mat) {
    n_local <- nrow(u_mat)

    # --------------------------------------------------
    # Transform u -> x (vectorized)
    # --------------------------------------------------

    p_mat <- pnorm(u_mat)
    p_mat <- pmin(pmax(p_mat, 1e-12), 1 - 1e-12)

    q_list <- lapply(lDistr, function(d) d$q)

    x_mat <- matrix(NA_real_, n_local, n_vars)

    for (j in seq_len(n_vars)) {
      x_mat[, j] <- q_list[[j]](p_mat[, j])
    }

    # --------------------------------------------------
    # log phi(u)  (standard normal)
    # --------------------------------------------------

    log_phi <- -0.5 * rowSums(u_mat^2) -
      0.5 * n_vars * log(2 * pi)

    # --------------------------------------------------
    # log f(x)
    # --------------------------------------------------

    d_list <- lapply(lDistr, function(d) d$d)

    log_fx_mat <- matrix(0, n_local, n_vars)

    for (j in seq_len(n_vars)) {
      dens_vals <- d_list[[j]](x_mat[, j])
      log_fx_mat[, j] <- log(pmax(dens_vals, .Machine$double.xmin))
    }

    log_fx <- rowSums(log_fx_mat)

    # --------------------------------------------------
    # log h(u)  (shifted normal N(dp, I))
    # --------------------------------------------------

    # broadcast dp to matrix form
    diff <- u_mat - matrix(dp,
      nrow = n_local,
      ncol = n_vars,
      byrow = TRUE
    )

    log_h <- -0.5 * rowSums(diff^2) -
      0.5 * n_vars * log(2 * pi)

    # --------------------------------------------------
    # log weights
    # --------------------------------------------------
    log_w <- log_fx - log_h
    log_w[!is.finite(log_w)] <- -Inf

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

    list(
      log_w  = log_w,
      I_vals = I_vals
    )
  }

  ## ---- Cluster setup --------------------------------------------------------
  cl <- NULL # Platzhalter, falls kein Cluster gebraucht wird

  ## ---- 1a. Seed fuer den Cluster bestimmen -------------------------------
  ## `streams` ist der Vektor, den `init_rng_master()` zurueckgegeben hat.
  ## Das erste Element ist ein gueltiger Integer-Seed.
  ## seed_used <- as.integer(streams[1])

  # streams <- init_rng_master(NULL)    # NULL -> neuer zufaelliger Seed pro Aufruf
  # seed_used <- as.integer(streams[1])

  ## -------------------------------------------------------------------------
  ## 2) Wahl des Parallel-Backends
  ## -------------------------------------------------------------------------
  if (backend == "parallel" && Sys.info()[["sysname"]] != "Windows") {
    ## --- RNG-Stream vorbereiten --------------------------------------------
    ## Fuer das Parallel-Backend brauchen wir den vollstaendigen CMRG-Stream
    if (is.null(streams_master)) streams_master <- parallel::nextRNGStream(.Random.seed)

    cl <- make_ready_cluster(
      use_threads = use_threads,
      libPaths_local = libPaths_local,
      seed = streams_master,
      export_objs = c(
        "libPaths_local", "mc_local", "lsf", "lDistr", "dp",
        "n_batch", "densityType", "n_vars"
      ),
      debug.level = debug.level
    )

    on.exit(parallel::stopCluster(cl), add = TRUE)
  } else if (backend == "future") {
    if (!requireNamespace("future", quietly = TRUE)) {
      stop("Package 'future' needed for backend = 'future'", call. = FALSE)
    }

    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Package 'future.apply' needed for backend = 'future'", call. = FALSE)
    }

    # ------------------------------------------------------------
    # CRAN-conform future handling
    # ------------------------------------------------------------
    # Save current plan
    old_plan <- future::plan()

    # Restore original plan when function exits
    on.exit(
      {
        future::plan(old_plan)
      },
      add = TRUE
    )

    # Only set a parallel plan if the user did not already set one
    if (inherits(old_plan, "sequential") && use_threads > 1) {
      if (.Platform$OS.type == "windows") {
        future::plan(future::multisession,
          workers = use_threads
        )
      } else {
        future::plan(future::multicore,
          workers = use_threads
        )
      }

      if (debug.level >= 1) {
        message(sprintf(
          "[MC_IS_single] Temporary future plan set with %d workers.",
          use_threads
        ))
      }
    } else {
      if (debug.level >= 1) {
        message("[MC_IS_single] Using existing future plan (no override).")
      }
    }
  } else {
    # stop("Unsupported backend")
    if (backend != "parallel" && backend != "future") stop("Unsupported backend", call. = FALSE)
  }

  ## -------------------------------------------------------------------------
  ## 3. Monte-Carlo-Loop (identisch fuer beide Back-ends)
  ## -------------------------------------------------------------------------
  n_sim <- 0L
  pf_mc <- 0
  var_mc <- 0
  w_mc <- 0
  I_mc <- 0
  var <- NA_real_
  cov <- 1
  pf <- NA_real_

  # Statistic container
  total_sum_Iw <- 0
  total_sum_w <- 0
  total_sum_w2 <- 0
  total_I_unw <- 0
  total_n <- 0

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

    # ============================================================
    # 3a. Deterministic sampling block (RNG only in master)
    # ============================================================

    # 1) Generate ALL samples in one deterministic block
    u_mat <- matrix(rnorm(total_batch * n_vars),
      ncol = n_vars
    )

    ### cat(
    ###   "After sampling RNG head:",
    ###   paste(.Random.seed[1:5], collapse = ", "),
    ###   "\n"
    ### )

    u_mat <- sweep(u_mat, 2, dp, "+")

    # 2) Deterministic split into worker chunks
    if (use_threads <= 1) {
      chunk_indices <- list(seq_len(total_batch))
    } else {
      chunk_indices <- split(
        seq_len(total_batch),
        cut(seq_len(total_batch),
          breaks = use_threads,
          labels = FALSE
        )
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

    # 4) Parallel evaluation (RNG disabled!)
    if (use_threads == 1) {
      res_list <- lapply(chunk_indices, worker_fun)
    } else if (backend == "parallel") {
      if (.Platform$OS.type == "windows") {
        res_list <- lapply(chunk_indices, worker_fun)
      } else {
        res_list <- parallel::parLapply(
          cl,
          chunk_indices,
          worker_fun
        )
      }
    } else if (backend == "future") {
      res_list <- future.apply::future_lapply(
        chunk_indices,
        worker_fun,
        future.seed = FALSE,
        future.globals = FALSE # CRITICAL
      )
    }

    # --- 3b. Aufsummieren ---------------------------------------------------------
    all_log_w <- unlist(lapply(res_list, `[[`, "log_w"))
    all_I_vals <- unlist(lapply(res_list, `[[`, "I_vals"))

    # --------------------------------------------------------
    # FAST MODE
    # --------------------------------------------------------
    if (stability_mode == "fast") {
      max_log_w <- max(all_log_w)
      w_all <- exp(all_log_w - max_log_w)

      # accumulate using base sum (fast, stable enough here)
      sum_w_iter <- sum(w_all)
      sum_Iw_iter <- sum(w_all[all_I_vals == 1])

      # second moment for ESS
      sum_w2_iter <- sum(w_all * w_all)
    }

    # --------------------------------------------------------
    # ROBUST MODE (numerically stable weight accumulation)
    # --------------------------------------------------------
    if (stability_mode == "robust") {
      # Compute global shift (log-sum-exp reference)
      # This prevents overflow when exponentiating large log-weights
      max_log_w <- max(all_log_w)

      # Shift log-weights before exponentiation
      # This is mathematically equivalent to fast-mode scaling
      shifted_weights <- exp(all_log_w - max_log_w)

      # First moment: sum of weights
      sum_w_iter <- sum(shifted_weights)

      # Weighted failure indicator
      if (any(all_I_vals == 1)) {
        sum_Iw_iter <- sum(shifted_weights[all_I_vals == 1])
      } else {
        sum_Iw_iter <- 0
      }

      # Second moment: required for ESS computation
      sum_w2_iter <- sum(shifted_weights * shifted_weights)
    }

    # global accumulation
    total_sum_w <- total_sum_w + sum_w_iter
    total_sum_Iw <- total_sum_Iw + sum_Iw_iter
    total_sum_w2 <- total_sum_w2 + sum_w2_iter

    total_n <- total_n + length(all_log_w)

    # --- 3c. Recording ------------------------------------------------------------
    # record totals for this outer iteration
    n_sim <- total_n
    if (dataRecord) {
      rec <- record_step(rec, iter, total_sum_Iw, total_sum_w2, total_I_unw, n_sim, tic)
    }

    # ------ 3d. Schaetzer berechnen ----------------------
    pf <- NA_real_
    var_est <- NA_real_
    cov <- Inf
    ESS <- 0

    if (total_sum_w > 0 && total_n > 0) {
      pf_hat <- total_sum_Iw / total_sum_w

      S_w2 <- total_sum_w2
      S_w <- total_sum_w
      N <- total_n

      if (total_sum_w > 0) {
        pf_hat <- total_sum_Iw / total_sum_w

        # Effective Sample Size
        ESS <- (total_sum_w^2) / total_sum_w2

        # Stable ESS-based variance approximation
        if (ESS > 0 && pf_hat > 0) {
          var_est <- pf_hat * (1 - pf_hat) / ESS
          cov <- sqrt(var_est) / pf_hat
        } else {
          var_est <- NA_real_
          cov <- Inf
        }
      } else {
        pf_hat <- 0
        var_est <- NA_real_
        cov <- Inf
        ESS <- 0
      }

      pf <- pf_hat
    } else {
      # no weight mass -> pf = 0
      pf <- 0
      var_est <- 0
      cov <- Inf
    }

    # debug printing (optional)
    if (debug.level >= 1) {
      message(sprintf("[MC_IS_single] iter=%d, n_sim=%d, pf=%.12g, cov=%.4f", iter, n_sim, pf, cov))
    }
    # convergence checks
    if (!is.na(cov) && is.finite(cov) && cov < cov_user) {
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
  beta_val <- if (!is.na(pf) && pf > 0) -stats::qnorm(pf) else NA_real_

  output <- list(
    method = "MCIS_Single",
    beta = beta_val, #-stats::qnorm(pf),
    pf = pf,
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


#' @title Update one recorder entry during main loop iteration
#' @description
#' Records current cumulative sums and timing information at iteration index \code{idx}.
#'
#' @inheritParams make_recorders
#'
#' @param idx Current iteration index within main loop.
#' @param pf_mc,var_mc,I_mc,n_sim Current cumulative sums from Monte Carlo batches.
#' @param tic Start time reference from \code{proc.time()} at simulation start.
#'
#' @return Updated recorder list object with new values inserted at position \code{idx}.
#'
#' @keywords internal
record_step <- function(rec, idx, pf_mc, var_mc, I_mc, n_sim, tic) {
  rec$pf[idx] <- pf_mc
  rec$var[idx] <- var_mc
  rec$I_sum[idx] <- I_mc
  rec$n_sim[idx] <- n_sim
  p_stamp <- proc.time() - tic
  rec$time$user[idx] <- p_stamp[1]
  rec$time$sys[idx] <- p_stamp[2]
  rec$time$elapsed[idx] <- p_stamp[3]
  rec$time$child_user[idx] <- p_stamp[4]
  rec$time$child_sys[idx] <- p_stamp[5]
  rec
}


#' @title Convert recorder lists into tidy data frame after completion
#' @description#'
#' Converts all recorded variables up to iteration \code{k} into a compact table format suitable for plotting or result export.
#'
#' Used internally after convergence or stop criterion is reached in Monte Carlo simulations.
#'
#' @param rec Recorder list created by [make_recorders()]
#' @param k Number of completed iterations/batches to include in output table.
#' @return Data frame containing columns: n_sim,pf,var,cov,I_sum,time_user,time_sys,time_elapsed,...
#' @keywords internal
#'
as_record_df <- function(rec, k) {
  if (k <= 0) {
    return(data.frame())
  }
  data.frame(
    n_sim = rec$n_sim[1:k],
    pf = rec$pf[1:k],
    var = rec$var[1:k],
    cov = rec$cov[1:k],
    I_sum = rec$I_sum[1:k],
    time_user = rec$time$user[1:k],
    time_sys = rec$time$sys[1:k],
    time_elapsed = rec$time$elapsed[1:k],
    time_user_child = rec$time$child_user[1:k],
    time_sys_child = rec$time$child_sys[1:k]
  )
}


## ---------------------------------------------------------------
## Parallel cluster setup utility used by MC_IS_single()
## ---------------------------------------------------------------

#' @title Prepare parallel cluster environment for MC Importance Sampling
#' @description Creates an R parallel cluster with specified number of workers,
#' exports required objects/functions and initializes reproducible RNG streams on each worker.
#'
#' @param use_threads Number of worker cores/threads to create via \link[parallel]{makeCluster}.
#' @param libPaths_local Library paths propagated to each worker via `.libPaths()`.
#' @param seed Integer base seed used for reproducible random numbers across workers.
#' @param export_objs Character vector naming objects/functions that should be exported from parent environment.
#'
#' @return Ready-to-use parallel cluster object configured with correct libraries and RNG streams.
#'
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport clusterEvalQ
#' @keywords internal
make_ready_cluster <- function(use_threads,
                               libPaths_local,
                               seed = NULL,
                               export_objs = character(0),
                               debug.level = 0) {
  # -------------------------------------------------
  # 1) Validate input
  # -------------------------------------------------
  stopifnot(is.numeric(use_threads), use_threads >= 1)
  stopifnot(is.character(export_objs))

  # -------------------------------------------------
  # 2) Choose cluster type depending on OS
  # -------------------------------------------------
  cl_type <- if (.Platform$OS.type == "windows") "PSOCK" else "FORK"

  cl <- parallel::makeCluster(use_threads, type = cl_type)
  # cl <- parallel::makeCluster(use_threads)

  # Ensure cluster is stopped on error
  on.exit(
    {
      try(parallel::stopCluster(cl), silent = TRUE)
    },
    add = TRUE
  )

  # -------------------------------------------------
  # 3) RNG handling
  # -------------------------------------------------
  if (!is.null(seed)) {
    ## seed ist ein kompletter CMRG-Stream (Laenge 7)
    parallel::clusterSetRNGStream(cl, iseed = seed)
    if (debug.level >= 1) {
      message("[make_ready_cluster] Fixed CMRG stream set on workers (backend parallel)")
    }
  } else {
    ## zufaellige Sub-Streams erzeugen
    parallel::clusterSetRNGStream(cl) # erzeugt automatisch neue Streams
    if (debug.level >= 1) {
      message("[make_ready_cluster] Random sub-streams generated (backend parallel)")
    }
  }

  #   ## 2) Export aller benoetigten Objekte
  #   all_export <- unique(c(export_objs, "libPaths_local"))
  #   parallel::clusterExport(cl,
  #     varlist = all_export,
  #     envir = parent.frame()
  #   )

  #   ## 3) Packages auf den Workern laden
  #   parallel::clusterEvalQ(cl, {
  #     .libPaths(libPaths_local)
  #     suppressMessages({
  #       if (!requireNamespace("TesiproV", quietly = TRUE)) library(TesiproV)
  #       if (!requireNamespace("edfun", quietly = TRUE)) library(edfun)
  #       if (!requireNamespace("evd", quietly = TRUE)) library(evd)
  #     })
  #   })
  #   cl
  # }

  # -------------------------------------------------
  # 4) Build clean export environment
  # -------------------------------------------------
  export_env <- new.env(parent = emptyenv())

  # Always export libPaths_local explicitly
  assign("libPaths_local", libPaths_local, envir = export_env)

  # Export additional requested objects if they exist
  if (length(export_objs) > 0) {
    for (obj in export_objs) {
      if (exists(obj, envir = parent.frame(), inherits = TRUE)) {
        assign(obj,
          get(obj, envir = parent.frame(), inherits = TRUE),
          envir = export_env
        )
      } else {
        warning(sprintf(
          "[make_ready_cluster] Object '%s' not found in calling environment - not exported.",
          obj
        ))
      }
    }
  }

  # Perform export
  parallel::clusterExport(
    cl,
    varlist = ls(export_env),
    envir = export_env
  )

  # -------------------------------------------------
  # 5) Initialize worker environments
  # -------------------------------------------------
  parallel::clusterEvalQ(cl, {
    .libPaths(libPaths_local)

    suppressMessages({
      if (!requireNamespace("TesiproV", quietly = TRUE)) {
        library(TesiproV)
      }

      if (!requireNamespace("edfun", quietly = TRUE)) {
        library(edfun)
      }

      if (!requireNamespace("evd", quietly = TRUE)) {
        library(evd)
      }
    })

    NULL
  })

  # Remove automatic stop handler (cluster valid)
  on.exit(NULL, add = FALSE)

  return(cl)
}

#' @title Resolve and validate RNG seed for Future backend
#' @description
#' Internal helper to normalize the \code{seed} argument used in
#' \link[future.apply]{future_lapply} calls within Monte Carlo simulations.
#'
#' The function ensures that the supplied seed is of a valid type:
#' either a single integer (for deterministic runs),
#' a complete L'Ecuyer-CMRG RNG stream (as returned by \link[parallel]{nextRNGStream}),
#' or \code{NULL}/\code{TRUE} for fully random behavior.
#'
#' This avoids common errors such as:
#' \itemize{
#'   \item "Argument 'seed' must be L'Ecuyer-CMRG RNG seed ..." when passing incomplete streams,
#'   \item identical random numbers on all workers due to improper seeding,
#'   \item non-reproducible results across different backends.
#' }
#'
#' @param seed Optional seed value provided by user or upstream function.
#' @description
#' May be one of:
#' \itemize{
#'   \item A single numeric/integer value (e.g. 1234) - fixed reproducible run.
#'   \item A numeric vector of length 7 (L'Ecuyer-CMRG stream) - substreams per worker.
#'   \item NULL - automatically use fully random streams (\code{future.seed = TRUE}).
#'  }
#'
#' @return A normalized object suitable for direct use in the argument
#'         \code{future.seed = resolve_future_seed(seed)} inside
#'         [future.apply::future_lapply()].
#'
#' @seealso [TesiproV::MC_IS()], [parallel::nextRNGStream()], [future.apply::future_lapply()]
#'
#' @author (C) 2026 M. Ricker - TU Dortmund University - Chair of Structural Concrete
#'
#' @keywords internal
resolve_future_seed <- function(seed) {
  ## --- Case 0: zero-length or NULL -> random automatic seeding ---
  if (is.null(seed) || length(seed) == 0L) {
    return(TRUE)
  }

  ## --- Case 1: Single integer/numeric ----------------------------------------
  if (is.numeric(seed) && length(seed) == 1L && !is.na(seed)) {
    return(as.integer(seed))
  }

  ## --- Case 2: Full L'Ecuyer-CMRG RNG stream ---------------------------------
  if (is.numeric(seed) && length(seed) == 7L && all(is.finite(seed))) {
    return(seed)
  }

  ## --- Otherwise invalid -----------------------------------------------------
  stop(
    sprintf(
      "Invalid 'seed': must be NULL/TRUE, single integer or L'Ecuyer-CMRG stream.\nGot type=%s length=%d.",
      typeof(seed), length(seed)
    ),
    call. = FALSE
  )
}


# =====================================================================
# MC_IS_system.R - Importance Sampling for a system of LSFs
# =====================================================================

#' Monte-Carlo Importance Sampling for a system of limit-state functions
#'
#' @inheritParams MC_IS_single
#' @param sys_type  Either "parallel" (default) or "serial"
#' @param beta_l  Optional threshold for reliability index screening.
#'   Limit-state functions with \eqn{\beta > beta_l} may be excluded
#'   from mixture construction in system analysis.
#' @export
MC_IS_system <- function(lsf,
                         lDistr,
                         cov_user,
                         n_batch,
                         n_max,
                         use_threads,
                         sys_type,
                         dataRecord,
                         beta_l,
                         densityType,
                         dps,
                         debug.level,
                         streams = NULL,
                         libPaths_local = .libPaths(),
                         seed = NULL,
                         backend = "future",
                         adaptive_alpha = FALSE,
                         alpha_update_rate = 0.1,
                         adaptive_batch = FALSE,
                         batch_control = list(),
                         min_adapt_samples = 10000,
                         alpha_min = 0.02,
                         stability_mode = c("robust", "fast")) {
  # Basic input checks
  stability_mode <- match.arg(stability_mode)

  # ---- CRAN safe future handling ----
  if (backend == "future") {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)

    if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
      future::plan(sequential)
    }
  }

  if (!is.list(lsf) || length(lsf) == 0L) {
    stop("MC_IS_system: 'lsf' must be a non-empty list of functions.")
  }

  if (!is.list(lDistr) || length(lDistr) != length(lsf)) {
    stop("MC_IS_system: 'lDistr' must be a list with one entry per LSF.")
  }

  if (!sys_type %in% c("serial", "parallel")) {
    stop("MC_IS_system: 'sys_type' must be 'serial' or 'parallel'.")
  }

  debug.TAG <- "MC_IS_System"
  if (debug.level >= 1) message("[MC_IS_system] starting system importance sampling...")

  tic <- proc.time()


  if (length(batch_control) == 0) {
    batch_control <- list(
      batch_min = 2000,
      batch_max = 100000,
      K_future = 5
    )
  }

  n_lsfs <- length(lsf)

  n_sim <- 0L
  total_sum_w <- 0
  total_sum_Iw <- 0
  total_sum_w2 <- 0

  ESS_global <- 0
  pf_hat <- NA_real_
  var_est <- NA_real_
  cov <- Inf

  # --------------------------------------------------
  # Normalize LSF input: accept functions or SYS_LSF objects
  # --------------------------------------------------
  lsf_fun <- lapply(lsf, function(obj) {
    if (is.function(obj)) {
      return(obj)
    } else if (!is.null(obj$func) && is.function(obj$func)) {
      return(obj$func)
    } else {
      stop("Invalid LSF format: must be function or SYS_LSF object.")
    }
  })


  # seeds setzen
  if (!is.null(seed)) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(as.integer(seed))
  }


  ## ------------------------------------------------------------------
  ## Normalize lDistr -> ensure each marginal is a funlist (with $d/$p/$q/$r)
  ## ------------------------------------------------------------------
  normalize_marginal <- function(m) {
    if (is.list(m) && !is.null(m$p) && is.function(m$p)) {
      return(m)
    } # already funlist
    if (is.list(m) && length(m) >= 1 && is.list(m[[1]]) && !is.null(m[[1]]$p) && is.function(m[[1]]$p)) {
      return(m[[1]])
    }
    stop("MC_IS_system: lDistr contains unexpected marginal format", call. = FALSE)
  }

  lDistr <- lapply(lDistr, function(per_lsf) {
    lapply(per_lsf, normalize_marginal)
  })

  ## ------------------------------------------------------------------
  ## 1) Build unique variable list and mappings
  ## ------------------------------------------------------------------

  var.names <- character(0)
  vars <- list() # store funlist per unique variable
  var.map <- integer(0)
  var.per.lsf <- integer(n_lsfs)

  for (i in seq_len(n_lsfs)) {
    nd <- length(lDistr[[i]])
    var.per.lsf[i] <- nd

    for (j in seq_len(nd)) {
      funlist <- lDistr[[i]][[j]]
      if (!is.list(funlist) || is.null(funlist$p) || !is.function(funlist$p)) {
        stop(sprintf("MC_IS_system: lDistr[[%d]][[%d]] is not a valid marginal (funlist)", i, j), call. = FALSE)
      }
      name_j <- funlist$name
      if (is.null(name_j) || name_j == "") name_j <- paste0("var_", length(var.names) + 1L)
      idx <- match(name_j, var.names, nomatch = 0L)
      if (idx == 0L) {
        var.names <- c(var.names, name_j)
        vars[[length(var.names)]] <- funlist
        var.map <- c(var.map, length(var.names))
      } else {
        var.map <- c(var.map, idx)
      }
    }
  }

  n_vars_unique <- length(vars)

  # var.in.lsf
  var.in.lsf <- vector("list", n_lsfs)
  k <- 1L
  for (i in seq_len(n_lsfs)) {
    var.in.lsf[[i]] <- var.map[k:(k + var.per.lsf[i] - 1L)]
    k <- k + var.per.lsf[i]
  }

  # lsf.in.var
  lsf.in.var <- vector("list", n_vars_unique)
  for (p in seq_len(n_vars_unique)) lsf.in.var[[p]] <- which(vapply(var.in.lsf, function(v) p %in% v, logical(1)))

  ## ------------------------------------------------------------------
  ## 2) FORM analysis -> design points per LSF
  ## ------------------------------------------------------------------
  res_form <- vector("list", n_lsfs)
  res_form.beta <- rep(NA_real_, n_lsfs)
  res_form.dp <- vector("list", n_lsfs)
  dp_u <- vector("list", n_lsfs)

  for (i in seq_len(n_lsfs)) {
    if (is.numeric(dps) && length(dps) == n_lsfs) {
      # --- user supplied design point ---
      x_star <- as.numeric(dps[i])

      Fxi <- tryCatch(lDistr[[i]][[1]]$p(x_star),
        error = function(e) NA_real_
      )

      Fxi <- min(max(as.numeric(Fxi), 1e-12), 1 - 1e-12)

      dp_u[[i]] <- qnorm(Fxi)

      res_form[[i]] <- list(
        beta = NA_real_,
        pf = NA_real_,
        x_points = x_star
      )

      res_form.beta[i] <- NA_real_
      res_form.dp[[i]] <- as.numeric(x_star)
    } else {
      # --- Funlists extrahieren ---
      lDistr_funlist_i <- lapply(lDistr[[i]], function(dd) {
        if (is.list(dd) && !is.null(dd$p) && is.function(dd$p)) {
          return(dd)
        }
        if (is.list(dd) && length(dd) >= 1 &&
          is.list(dd[[1]]) && !is.null(dd[[1]]$p)) {
          return(dd[[1]])
        }
        stop("Invalid marginal format for FORM")
      })

      # --- Wrapper bauen ---
      lDistr_for_FORM_i <- lapply(lDistr_funlist_i, function(fun) {
        params <- c(NA_real_, NA_real_)
        if (!is.null(fun$mean)) params[1] <- as.numeric(fun$mean)
        if (!is.null(fun$Sd)) params[2] <- as.numeric(fun$Sd)
        list(fun, params)
      })

      # --- FORM aufrufen ---
      res_form[[i]] <- TesiproV::FORM(
        lsf_fun[[i]],
        lDistr_for_FORM_i,
        n_optim = 20,
        loctol  = 1e-4
      )

      res_form.beta[i] <- res_form[[i]]$beta
      res_form.dp[[i]] <- as.numeric(unlist(res_form[[i]]$x_points))

      # --- U-Transformation ---
      x_star_vec <- res_form.dp[[i]]

      dp_u[[i]] <- numeric(length(lDistr_funlist_i))

      for (kidx in seq_along(lDistr_funlist_i)) {
        xsk <- as.numeric(x_star_vec[kidx])

        Fxi <- tryCatch(
          lDistr_funlist_i[[kidx]]$p(xsk),
          error = function(e) NA_real_
        )

        Fxi <- min(max(as.numeric(Fxi), 1e-12), 1 - 1e-12)

        dp_u[[i]][kidx] <- qnorm(Fxi)
      }
    }
  }

  # ------------------------------------------------------------
  # Automatic vectorization detection for each LSF
  # ------------------------------------------------------------

  lsf_vectorized_list <- logical(n_lsfs)

  for (k in seq_len(n_lsfs)) {
    # small test matrix (2 samples, correct variable subset)
    idxs <- var.in.lsf[[k]]

    x_star <- res_form.dp[[k]]

    test_mat <- matrix(
      rep(x_star, 2),
      nrow = 2,
      byrow = TRUE
    )

    test_result <- try(lsf_fun[[k]](test_mat), silent = TRUE)

    if (!inherits(test_result, "try-error") &&
      is.numeric(test_result) &&
      length(test_result) == 2 &&
      all(is.finite(test_result))) {
      lsf_vectorized_list[k] <- TRUE

      if (debug.level >= 1) {
        message(sprintf("[MC_IS_system] LSF %d detected as vectorized.", k))
      }
    }
  }

  # --------------------------------------------------
  # Globale Einbettung der Designpunkte im U-Raum
  # --------------------------------------------------

  dp_u_global <- vector("list", n_lsfs)

  for (k in seq_len(n_lsfs)) {
    u_star <- numeric(n_vars_unique)
    idxs <- var.in.lsf[[k]]
    u_star[idxs] <- dp_u[[k]]
    dp_u_global[[k]] <- u_star
  }

  if (debug.level >= 2) {
    dp_glob_str <- paste(
      sapply(seq_along(dp_u_global), function(i) {
        paste0(
          "LSF", i, ": (",
          paste(round(dp_u_global[[i]], 4), collapse = ", "),
          ")"
        )
      }),
      collapse = " | "
    )

    message(sprintf("[MC_IS_system] Global U-shifts: %s", dp_glob_str))
  }

  # ------------------------------------------------------------------
  # Parallel system: compute system design point
  # ------------------------------------------------------------------
  if (sys_type == "parallel") {
    g_sys <- function(x) {
      vals <- numeric(n_lsfs)
      for (i in seq_len(n_lsfs)) {
        vals[i] <- lsf_fun[[i]](x[var.in.lsf[[i]]])
      }
      min(vals)
    }

    lDistr_global <- lapply(vars, function(fun) {
      list(fun, c(fun$mean, fun$Sd))
    })

    res_form_sys <- TesiproV::FORM(
      g_sys,
      lDistr_global,
      n_optim = 20,
      loctol = 1e-4
    )

    x_star_sys <- as.numeric(res_form_sys$x_points)

    dp_u_sys <- numeric(n_vars_unique)
    for (j in seq_len(n_vars_unique)) {
      Fxi <- vars[[j]]$p(x_star_sys[j])
      dp_u_sys[j] <- qnorm(Fxi)
    }
  }

  # ------------------------------------------------------------------
  # Initial mixture weights (serial only)
  # ------------------------------------------------------------------
  if (sys_type == "serial") {
    pf_single <- pnorm(-res_form.beta)
    alpha_i <- pf_single / sum(pf_single)
  }

  ## ------------------------------------------------------------------
  ## 3) Build dp matrix (global)
  ## ------------------------------------------------------------------
  dp <- matrix(NA_real_, nrow = n_vars_unique, ncol = n_lsfs)
  for (i in seq_len(n_lsfs)) {
    local_vars <- var.in.lsf[[i]]
    if (!is.null(res_form.dp[[i]])) {
      for (kl in seq_along(local_vars)) {
        global_id <- local_vars[kl]
        if (kl <= length(res_form.dp[[i]])) dp[global_id, i] <- res_form.dp[[i]][kl]
      }
    }
  }

  dp.means <- rowMeans(dp, na.rm = TRUE)

  ## ------------------------------------------------------------------
  ## 4) ai weights
  ## ------------------------------------------------------------------
  ai <- matrix(0, nrow = n_vars_unique, ncol = n_lsfs)
  for (p in seq_len(n_vars_unique)) {
    comps <- lsf.in.var[[p]]
    if (length(comps) == 0L) next
    betas <- res_form.beta[comps]
    betas[!is.finite(betas) | betas <= 0] <- Inf
    inv_b <- 1 / betas
    s <- sum(inv_b, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) ai[p, comps] <- 1 / length(comps) else ai[p, comps] <- inv_b / s
  }

  ## ------------------------------------------------------------------
  ## 5) Build hv mixtures per global variable
  ## ------------------------------------------------------------------
  hv_list <- vector("list", n_vars_unique)
  for (p in seq_len(n_vars_unique)) {
    comps <- lsf.in.var[[p]]
    if (length(comps) == 0L) {
      hv_list[[p]] <- list(rfun = function(n) rep(NA_real_, n), dfun = function(x) rep(0, length(x)))
      next
    }
    orig_funlist <- vars[[p]] # funlist (d,p,q,r,mean,Sd,...)
    sd_p <- orig_funlist$Sd
    comp_dfun <- list()
    comp_rfun <- list()
    comp_probs <- numeric(length(comps))
    for (ci in seq_along(comps)) {
      j <- comps[ci]
      mu_j <- dp[p, j]
      sd_loc <- ifelse(is.finite(sd_p) && sd_p > 0, sd_p, max(1e-8, abs(mu_j) * 0.1, 1e-3))
      dist_type <- tolower(if (is.null(orig_funlist$DistributionType) || orig_funlist$DistributionType == "") "norm" else orig_funlist$DistributionType)
      if (dist_type %in% c("lnorm", "lt", "weibull")) {
        if (!is.finite(mu_j) || mu_j <= 0) {
          mn <- mu_j
          sn <- sd_loc
          comp_rfun[[ci]] <- function(n, m = mn, s = sn) stats::rnorm(n, mean = m, sd = s)
          comp_dfun[[ci]] <- function(x, m = mn, s = sn) stats::dnorm(x, mean = m, sd = s)
        } else {
          sd_rel <- sd_loc / mu_j
          sn <- sqrt(log(1 + sd_rel^2))
          mn <- log(mu_j / sqrt(1 + sd_rel^2))
          comp_rfun[[ci]] <- function(n, mn = mn, sn = sn) stats::rlnorm(n, meanlog = mn, sdlog = sn)
          comp_dfun[[ci]] <- function(x, mn = mn, sn = sn) stats::dlnorm(x, meanlog = mn, sdlog = sn)
        }
      } else {
        comp_rfun[[ci]] <- function(n, m = mu_j, s = sd_loc) stats::rnorm(n, mean = m, sd = s)
        comp_dfun[[ci]] <- function(x, m = mu_j, s = sd_loc) stats::dnorm(x, mean = m, sd = s)
      }
      comp_probs[ci] <- ai[p, j]
    }
    if (!is.finite(sum(comp_probs)) || sum(comp_probs) <= 0) comp_probs <- rep(1 / length(comp_probs), length(comp_probs)) else comp_probs <- comp_probs / sum(comp_probs)
    hv_dfun <- function(x) {
      res <- numeric(length(x))
      for (ci in seq_along(comp_dfun)) res <- res + comp_probs[ci] * comp_dfun[[ci]](x)
      res
    }
    hv_rfun <- function(n) {
      if (n <= 0) {
        return(numeric(0))
      }
      comps_drawn <- sample(seq_along(comp_probs), size = n, replace = TRUE, prob = comp_probs)
      res <- numeric(n)
      for (ci in unique(comps_drawn)) {
        idx <- which(comps_drawn == ci)
        res[idx] <- comp_rfun[[ci]](length(idx))
      }
      res
    }
    hv_list[[p]] <- list(rfun = hv_rfun, dfun = hv_dfun, comps = comps, probs = comp_probs)
  }

  ## ------------------------------------------------------------------
  ## 6) Monte Carlo Loop (single-sample per iter)
  ## ------------------------------------------------------------------
  I_sum <- 0L
  pf_acc <- 0
  var_acc <- 0
  n_sim <- 0L
  cov <- Inf


  if (dataRecord) {
    max_iter_est <- ceiling(n_max)
    rec_pf <- numeric(max_iter_est)
    rec_var <- numeric(max_iter_est)
    rec_cov <- numeric(max_iter_est)
    rec_nsim <- integer(max_iter_est)
    rec_I_sum <- integer(max_iter_est)
    rec_time_user <- numeric(max_iter_est)
    rec_time_sys <- numeric(max_iter_est)
    rec_time_elapsed <- numeric(max_iter_est)
  }

  orig_d_list <- lapply(vars, function(v) v$d)


  # --------------------------------------------------------------
  # Monte Carlo Loop (worker-invariant, SNIS, global stabilization)
  # --------------------------------------------------------------

  # --------------------------------------------------------------
  # Monte Carlo Loop - Multimodal, worker-invariant, global log-stable
  # --------------------------------------------------------------

  if (!exists("dp_u")) {
    stop("dp_u not constructed before Monte Carlo loop")
  }

  n_sim <- 0L
  total_sum_Iw <- 0
  total_sum_w <- 0
  total_sum_w2 <- 0

  ESS_modes <- NULL


  if (stability_mode == "robust") {
    log_total_sum_w <- -Inf
    log_total_sum_Iw <- -Inf
    log_total_sum_w2 <- -Inf
  }

  total_batch <- n_batch
  iter <- 0L

  if (sys_type == "serial") {
    if (all(is.finite(res_form.beta))) {
      pf_single <- pnorm(-res_form.beta)
      alpha_i <- pf_single / sum(pf_single)
    } else {
      alpha_i <- rep(1 / n_lsfs, n_lsfs)
    }
  }

  state <- list(
    total_sum_w = total_sum_w,
    total_sum_Iw = total_sum_Iw,
    total_sum_w2 = total_sum_w2,
    log_total_sum_w = if (stability_mode == "robust") -Inf else NULL,
    log_total_sum_Iw = if (stability_mode == "robust") -Inf else NULL,
    log_total_sum_w2 = if (stability_mode == "robust") -Inf else NULL,
    n_sim = n_sim
  )


  while (TRUE) {
    iter <- iter + 1L

    # ==========================================================
    # Sampling Block
    # ==========================================================

    if (sys_type == "serial") {
      # Serial: multimodal sampling
      mode_ids <- sample(seq_len(n_lsfs),
        size = total_batch,
        replace = TRUE,
        prob = alpha_i
      )

      u_mat <- matrix(0, nrow = total_batch, ncol = n_vars_unique)

      for (i in seq_len(total_batch)) {
        shift <- dp_u_global[[mode_ids[i]]]
        u_mat[i, ] <- rnorm(n_vars_unique,
          mean = shift,
          sd = 1
        )
      }
    } else { # parallel

      # Parallel: single shift at system design point
      u_mat <- matrix(rnorm(total_batch * n_vars_unique),
        ncol = n_vars_unique
      )

      u_mat <- sweep(u_mat, 2, dp_u_sys, "+")

      # Dummy mode ids (needed for ESS calculation)
      mode_ids <- rep(1L, total_batch)
    }

    # ==========================================================
    # Parallel Evaluation Block
    # ==========================================================

    chunk_ids <- rep(seq_len(use_threads),
      length.out = total_batch
    )

    chunk_indices <- split(seq_len(total_batch), chunk_ids)


    worker_fun <- function(idx) {
      u_block <- u_mat[idx, , drop = FALSE]
      n_local <- nrow(u_block)

      # ----------------------------
      # 1) Transform to physical space
      # ----------------------------

      p_mat <- pnorm(u_block)

      x_mat <- matrix(NA_real_, n_local, n_vars_unique)

      for (j in seq_len(n_vars_unique)) {
        x_mat[, j] <- vars[[j]]$q(p_mat[, j])
      }

      # ----------------------------
      # 2) Evaluate LSFs
      # ----------------------------

      I_matrix <- matrix(0L, n_local, n_lsfs)

      for (k in seq_len(n_lsfs)) {
        idxs <- var.in.lsf[[k]]
        x_sub <- x_mat[, idxs, drop = FALSE]

        fun_k <- lsf_fun[[k]]

        if (lsf_vectorized_list[k]) {
          vals <- fun_k(x_sub)
        } else {
          vals <- apply(x_sub, 1, fun_k)
        }

        I_matrix[, k] <- as.integer(vals < 0)
      }

      if (sys_type == "serial") {
        I_vals <- as.integer(rowSums(I_matrix) > 0)
      } else {
        I_vals <- as.integer(rowSums(I_matrix) == n_lsfs)
      }

      # ----------------------------
      # 3) log phi(u)
      # ----------------------------

      log_phi <- rowSums(dnorm(u_block, log = TRUE))

      # ----------------------------
      # 4) log h(u)
      # ----------------------------

      if (sys_type == "serial") {
        log_h_mat <- matrix(NA_real_, n_local, n_lsfs)

        for (k in seq_len(n_lsfs)) {
          shift <- dp_u_global[[k]]
          shifted <- sweep(u_block, 2, shift, "-")

          log_h_mat[, k] <-
            log(alpha_i[k]) +
            rowSums(dnorm(shifted, log = TRUE))
        }

        max_log <- apply(log_h_mat, 1, max)

        log_h <- max_log +
          log(rowSums(exp(log_h_mat - max_log)))
      } else {
        shifted <- sweep(u_block, 2, dp_u_sys, "-")
        log_h <- rowSums(dnorm(shifted, log = TRUE))
      }

      log_w <- log_phi - log_h

      list(
        log_w = log_w,
        I_vals = I_vals,
        I_matrix = I_matrix,
        mode_ids = mode_ids[idx]
      )
    }


    # --------------------------------------------------
    # Parallele Auswertung
    # --------------------------------------------------

    if (backend == "future") {
      res_list <- future.apply::future_lapply(
        chunk_indices,
        worker_fun,
        future.seed = FALSE
      )
    } else {
      res_list <- lapply(chunk_indices, worker_fun)
    }

    # ==========================================================
    # Weight Accumulation Block (refactored)
    # ==========================================================

    all_log_w <- unlist(lapply(res_list, `[[`, "log_w"))
    all_I_vals <- unlist(lapply(res_list, `[[`, "I_vals"))
    all_I_matrix <- do.call(rbind, lapply(res_list, `[[`, "I_matrix"))
    all_mode_ids <- unlist(lapply(res_list, `[[`, "mode_ids"))

    state <- compute_weight_update(
      all_log_w,
      all_I_vals,
      stability_mode,
      sys_type,
      state
    )

    total_sum_w <- state$total_sum_w
    total_sum_Iw <- state$total_sum_Iw
    total_sum_w2 <- state$total_sum_w2
    log_total_sum_w <- state$log_total_sum_w
    log_total_sum_Iw <- state$log_total_sum_Iw
    log_total_sum_w2 <- state$log_total_sum_w2
    n_sim <- state$n_sim

    # ==========================================================
    # Estimator Update Block
    # ==========================================================

    if (sys_type == "serial" && stability_mode == "fast") {
      max_log_w <- max(all_log_w)
      shifted_weights <- exp(all_log_w - max_log_w)

      ESS_modes <- numeric(n_lsfs)

      for (k in seq_len(n_lsfs)) {
        idx_k <- which(all_mode_ids == k)

        if (length(idx_k) > 0) {
          w_k <- shifted_weights[idx_k]
          ESS_modes[k] <- (sum(w_k)^2) / sum(w_k^2)
        } else {
          ESS_modes[k] <- 0
        }
      }
    }
    # --------------------------------------------------
    # Adaptive batch size update
    # --------------------------------------------------
    if (adaptive_batch && ESS_global > 0) {
      RSE_target <- cov_user
      K_future <- batch_control$K_future
      n_min_batch <- batch_control$n_min
      n_max_batch <- batch_control$n_max

      # Estimate current efficiency factor
      kappa_hat <- ESS_global / n_sim

      if (kappa_hat > 0) {
        ESS_target <- 1 / (RSE_target^2)
        N_target <- ESS_target / kappa_hat
        N_remaining <- max(N_target - n_sim, 0)

        n_batch_new <- ceiling(N_remaining / K_future)

        # Ensure divisibility by threads
        n_batch_new <- ceiling(n_batch_new / use_threads) * use_threads

        # Apply bounds
        n_batch <- max(n_min_batch, min(n_batch_new, n_max_batch))

        if (debug.level >= 2) {
          message(sprintf(
            "[MC_IS_system] Adaptive batch updated to %d (ESS=%.2f)",
            n_batch, ESS_global
          ))
        }
      }
    }


    # ==========================================================
    # Estimator Update Block
    # ==========================================================

    # FAST MODE or SERIAL SYSTEM
    if (stability_mode == "fast" || sys_type == "serial") {
      if (total_sum_w > 0) {
        pf_hat <- total_sum_Iw / total_sum_w
        ESS_global <- (total_sum_w^2) / total_sum_w2
      } else {
        pf_hat <- 0
        ESS_global <- 0
      }
    }

    # ROBUST MODE (parallel only)
    if (stability_mode == "robust" && sys_type == "parallel") {
      if (log_total_sum_w > -Inf) {
        pf_hat <- exp(log_total_sum_Iw - log_total_sum_w)
        ESS_global <- exp(2 * log_total_sum_w - log_total_sum_w2)
      } else {
        pf_hat <- 0
        ESS_global <- 0
      }
    }

    # Variance (common)
    if (ESS_global > 0 && pf_hat > 0) {
      var_est <- pf_hat * (1 - pf_hat) / ESS_global
      cov <- sqrt(var_est) / pf_hat
    } else {
      var_est <- NA_real_
      cov <- Inf
    }

    if (ESS_global > 0) {
      RSE_est <- 1 / sqrt(ESS_global)
    } else {
      RSE_est <- Inf
    }


    if (adaptive_alpha && sys_type == "serial" && n_sim > min_adapt_samples) {
      pf_i_hat <- numeric(n_lsfs)

      if (stability_mode == "fast") {
        pf_sys_hat <- sum(shifted_weights[all_I_vals == 1])

        if (!is.na(pf_sys_hat) && is.finite(pf_sys_hat) && pf_sys_hat > 0) {
          for (k in seq_len(n_lsfs)) {
            pf_i_hat[k] <-
              sum((all_I_matrix[, k] * all_I_vals) *
                shifted_weights) / pf_sys_hat
          }
        }
      }

      if (stability_mode == "robust") {
        # pf_sys_hat in log-domain
        pf_sys_hat <- exp(log_total_sum_Iw - log_total_sum_w)

        if (!is.na(pf_sys_hat) && is.finite(pf_sys_hat) && pf_sys_hat > 0) {
          for (k in seq_len(n_lsfs)) {
            idx_k <- which(all_I_matrix[, k] == 1 &
              all_I_vals == 1)

            if (length(idx_k) > 0) {
              log_contrib <- all_log_w[idx_k]

              max_log_c <- max(log_contrib)

              pf_i_hat[k] <-
                exp(
                  max_log_c +
                    log(sum(exp(log_contrib - max_log_c))) -
                    log_total_sum_w
                )
            } else {
              pf_i_hat[k] <- 0
            }
          }
        }
      }

      # update alpha
      if (sum(pf_i_hat) > 0) {
        alpha_target <- pf_i_hat / sum(pf_i_hat)

        alpha_i <- (1 - alpha_update_rate) * alpha_i +
          alpha_update_rate * alpha_target

        alpha_i <- pmax(alpha_i, alpha_min)
        alpha_i <- alpha_i / sum(alpha_i)
      }
    }


    if (debug.level >= 2) {
      if (adaptive_alpha) {
        print(alpha_i)
      }
    }

    # ==========================================================
    # Convergence Check
    # ==========================================================
    if (!exists("ESS_global")) ESS_global <- 0

    ESS_min <- if (sys_type == "parallel") 50 else 30
    n_min <- 5 * n_batch

    stop_condition <- (
      !is.na(cov) &&
        is.finite(cov) &&
        cov < cov_user &&
        ESS_global > ESS_min &&
        n_sim > n_min
    )

    if (stop_condition) {
      break
    }

    if (n_sim >= n_max) {
      break
    }
  }

  beta_final <- if (pf_hat > 0) -qnorm(pf_hat) else NA_real_

  out <- list(
    method    = "MCIS_System",
    sys_type  = sys_type,
    beta      = beta_final,
    pf        = pf_hat,
    var       = var_est,
    cov_mc    = cov,
    RSE_est   = RSE_est,
    ESS       = ESS_global,
    ESS_modes = if (sys_type == "serial") ESS_modes else NULL,
    n_mc      = n_sim
  )

  if (sys_type == "parallel" && ESS_global < ESS_min) {
    warning(sprintf(
      "Parallel system: Effective sample size (ESS = %.1f) is low. Relative standard error ~= %.2f%%.",
      ESS_global,
      100 * RSE_est
    ))
  }

  if (debug.level >= 1) message("[MC_IS_system] finished.")

  return(out)
}


#' @title Internal weight accumulation for Monte Carlo Importance Sampling (System Case)
#'
#' @description
#' Internal helper used by `MC_IS_system()` to accumulate importance sampling
#' weights and update running Monte Carlo statistics.
#'
#' The function performs numerically stable accumulation of:
#' \itemize{
#'   \item \eqn{\sum w},
#'   \item \eqn{\sum I w},
#'   \item \eqn{\sum w^2},
#'   \item and (in robust mode) their log-domain equivalents.
#' }
#'
#' Two stabilization modes are supported:
#'
#' \describe{
#'   \item{\code{"fast"}}{
#'     Batch-level stabilization using a max-log shift:
#'     \deqn{w_i = \exp(\log w_i - \max(\log w))}
#'     Suitable for most engineering applications.
#'   }
#'   \item{\code{"robust"}}{
#'     Global log-sum-exp accumulation across iterations:
#'     \deqn{\log \sum w = m + \log \sum \exp(\log w - m)}
#'     Recommended for rare-event simulation or very small failure probabilities.
#'   }
#' }
#'
#' The function updates the running state object and returns it.
#'
#' @param all_log_w Numeric vector of log-weights for the current batch.
#' @param all_I_vals Numeric vector of failure indicators (0/1).
#' @param stability_mode Character string: `"fast"` or `"robust"`.
#' @param sys_type Character string: `"serial"` or `"parallel"`.
#' @param state List containing current accumulated sums:
#'   \describe{
#'     \item{total_sum_w}{Cumulative sum of weights.}
#'     \item{total_sum_Iw}{Cumulative weighted failure sum.}
#'     \item{total_sum_w2}{Cumulative sum of squared weights.}
#'     \item{log_total_sum_w}{Log-domain cumulative weight sum (robust mode).}
#'     \item{log_total_sum_Iw}{Log-domain cumulative failure-weight sum (robust mode).}
#'     \item{log_total_sum_w2}{Log-domain cumulative squared-weight sum (robust mode).}
#'     \item{n_sim}{Total number of samples processed.}
#'   }
#'
#' @return Updated `state` list with accumulated statistics.
#'
#' @details
#' This function is internal and not intended to be called directly by users.
#' It is part of the system-level importance sampling implementation
#' in `MC_IS_system()`.
#'
#' @keywords internal
compute_weight_update <- function(all_log_w,
                                  all_I_vals,
                                  stability_mode,
                                  sys_type,
                                  state) {
  if (stability_mode == "fast" || sys_type == "serial") {
    max_log_w <- max(all_log_w)
    shifted_weights <- exp(all_log_w - max_log_w)

    sum_w_iter <- sum(shifted_weights)
    sum_Iw_iter <- sum(shifted_weights[all_I_vals == 1])
    sum_w2_iter <- sum(shifted_weights^2)

    state$total_sum_w <- state$total_sum_w + sum_w_iter
    state$total_sum_Iw <- state$total_sum_Iw + sum_Iw_iter
    state$total_sum_w2 <- state$total_sum_w2 + sum_w2_iter
    state$n_sim <- state$n_sim + length(shifted_weights)

    return(state)
  }

  if (stability_mode == "robust" && sys_type == "parallel") {
    log_w_batch <- all_log_w

    max_log_w <- max(log_w_batch)
    log_sum_w_batch <- max_log_w +
      log(sum(exp(log_w_batch - max_log_w)))

    if (any(all_I_vals == 1)) {
      log_Iw_batch <- log_w_batch[all_I_vals == 1]
      max_log_Iw <- max(log_Iw_batch)
      log_sum_Iw_batch <- max_log_Iw +
        log(sum(exp(log_Iw_batch - max_log_Iw)))
    } else {
      log_sum_Iw_batch <- -Inf
    }

    log_w2_batch <- 2 * log_w_batch
    max_log_w2 <- max(log_w2_batch)
    log_sum_w2_batch <- max_log_w2 +
      log(sum(exp(log_w2_batch - max_log_w2)))

    max_global_w <- max(state$log_total_sum_w, log_sum_w_batch)
    state$log_total_sum_w <- max_global_w +
      log(exp(state$log_total_sum_w - max_global_w) +
        exp(log_sum_w_batch - max_global_w))

    max_global_Iw <- max(state$log_total_sum_Iw, log_sum_Iw_batch)
    state$log_total_sum_Iw <- max_global_Iw +
      log(exp(state$log_total_sum_Iw - max_global_Iw) +
        exp(log_sum_Iw_batch - max_global_Iw))

    max_global_w2 <- max(state$log_total_sum_w2, log_sum_w2_batch)
    state$log_total_sum_w2 <- max_global_w2 +
      log(exp(state$log_total_sum_w2 - max_global_w2) +
        exp(log_sum_w2_batch - max_global_w2))

    state$n_sim <- state$n_sim + length(log_w_batch)

    return(state)
  }

  return(state)
}
