#' @name MC_IS
#' @title Monte Carlo Simulation with Importance Sampling
#' @description
#' Performs Monte Carlo simulation with importance sampling (IS) to estimate
#' the failure probability of a single limit-state function or of a system of
#' limit-state functions.
#'
#' Failure is defined as \eqn{g(x) < 0}. Samples with \eqn{g(x) = 0} are
#' counted as safe.
#'
#' The samples are drawn from a sampling density \eqn{h} in standard normal
#' space that is centred at FORM design points. The sampling is corrected by
#' the likelihood ratio \eqn{w(u) = \varphi(u) / h(u)}, where \eqn{\varphi} is
#' the standard normal density.
#'
#' @details
#' \strong{Estimator}
#'
#' The failure probability and the variance of the estimator are
#' \deqn{\hat{P}_f = \frac{1}{N} \sum_{i=1}^{N} I_i w_i, \qquad
#'   \widehat{\mathrm{Var}}(\hat{P}_f) =
#'   \frac{1}{N}\left(\frac{1}{N}\sum_{i=1}^{N} I_i w_i^2 - \hat{P}_f^2\right),}
#' where \eqn{I_i = 1} if sample \eqn{i} lies in the failure domain and
#' \eqn{I_i = 0} otherwise. The coefficient of variation is
#' \eqn{\sqrt{\widehat{\mathrm{Var}}(\hat{P}_f)} / \hat{P}_f}.
#'
#' \strong{Single limit-state function}
#'
#' The failure domain is \eqn{\{g(x) < 0\}}. The sampling density is a
#' multivariate normal density centred at the FORM design point \eqn{u^*}
#' with the covariance matrix
#' \deqn{\Sigma = I + (\beta^2 - 1)\, \alpha \alpha^T, \qquad
#'   \alpha = u^* / \|u^*\|,}
#' where \eqn{\beta} is the FORM reliability index, or \eqn{\|u^*\|} if it is
#' not available. The variance is \eqn{\beta^2} in the direction of
#' \eqn{\alpha} and 1 perpendicular to it. FORM is only used to locate the
#' centre. A design point from a FORM run that did not converge is accepted as
#' long as it is finite.
#'
#' \strong{Systems}
#'
#' For limit-state functions \eqn{g_1, \dots, g_m} the failure domain is
#' \eqn{\bigcup_k \{g_k(x) < 0\}} for a serial system and
#' \eqn{\bigcap_k \{g_k(x) < 0\}} for a parallel system. Basic variables that
#' appear in several limit-state functions are identified by their name and
#' enter the simulation only once.
#'
#' The sampling density is a mixture of unit-variance normal densities
#' \deqn{h(u) = \sum_{j} \lambda_j \, \varphi(u - u_j^*), \qquad
#'   \sum_j \lambda_j = 1.}
#' \itemize{
#'   \item Serial systems: the centres \eqn{u_j^*} are the FORM design points
#'     of the components. The weights are proportional to
#'     \eqn{\Phi(-\beta_j)} if all FORM results are finite, otherwise they are
#'     equal. If the largest distance between the design points is smaller
#'     than 1.5 and the smallest cosine between their gradient directions is
#'     larger than 0.85, a single normal density at the mean of the design
#'     points is used instead of the mixture.
#'   \item Parallel systems: FORM is run for the joint failure domain
#'     \eqn{\max_k g_k(x) < 0}. If its design point is available, it is a
#'     centre with weight 0.7 and the remaining weight is divided equally
#'     among the component design points. Otherwise the component design
#'     points are used with equal weights.
#' }
#'
#' \strong{Stopping criterion}
#'
#' For a single limit-state function the simulation stops when the
#' coefficient of variation falls below \code{cov_user} or when \code{n_max}
#' samples are reached. The last batch is reduced so that the total does not
#' exceed \code{n_max}.
#'
#' For systems the simulation stops when the coefficient of variation falls
#' below \code{cov_user}, the effective number of failure samples exceeds 30
#' (serial) or 50 (parallel) and more than \code{5 * n_batch} samples have been
#' used. It also stops when \code{n_max} is reached. The batches are not
#' reduced, so \code{n_mc} can exceed \code{n_max} by less than one batch. For
#' parallel systems a warning is issued if the effective number of failure
#' samples is below 50 at the end.
#'
#' \strong{Numerical stability}
#'
#' The weighted sums are accumulated in the log domain. The result includes
#' \code{log_pf}. For extremely small probabilities \code{pf} may underflow
#' to zero, \code{beta} is then computed from \code{log_pf} and can remain
#' finite. This applies to single limit-state functions and to systems.
#'
#' \strong{Automatic detection of vectorized limit-state functions}
#'
#' Before the simulation the limit-state function is called once with a
#' matrix of two rows (the design point in physical space). If it returns a
#' finite numeric vector of length 2, it is evaluated for whole batches at
#' once. Otherwise it is evaluated row by row with \code{apply()}. Scalar
#' limit-state functions therefore need not be changed.
#'
#' \strong{Parallel execution}
#'
#' The samples of each batch are generated in the master process and split
#' into \code{use_threads} chunks, which are evaluated in parallel. The workers
#' do not generate random numbers.
#' \itemize{
#'   \item \code{backend = "future"}: the chunks are evaluated with
#'     \code{future.apply::future_lapply()} under the active \pkg{future}
#'     plan. The function never changes the plan. If the plan is sequential and
#'     \code{use_threads > 1}, a message is shown once per session and the
#'     calculation runs sequentially. The number of workers is determined by
#'     the plan. \code{use_threads} only determines the number of chunks.
#'   \item \code{backend = "parallel"}: on Unix-like systems
#'     \code{parallel::mclapply()} with \code{use_threads} cores is used. On
#'     Windows a PSOCK cluster with \code{use_threads} workers is created for
#'     the call and stopped afterwards. With \code{use_threads = 1}
#'     \code{lapply()} is used.
#' }
#'
#' \strong{Reproducibility}
#'
#' With a fixed \code{seed} the random numbers are generated from a
#' L'Ecuyer-CMRG stream in the master process. The sample therefore does not
#' depend on \code{use_threads}. Results for different values of
#' \code{use_threads} agree up to rounding errors, because the sums are
#' accumulated per chunk. The function sets \code{RNGkind("L'Ecuyer-CMRG")}
#' and changes \code{.Random.seed} in the global environment.
#'
#' @param lsf A limit-state function or a list of limit-state functions for a
#'   system. Must follow the LSF interface described below.
#'
#' @param lDistr List of marginal distribution objects (funlists with
#'   \code{$d}, \code{$p}, \code{$q}, \code{$r}). For a system, a list with one
#'   such list per limit-state function.
#'
#' @param cov_user Target coefficient of variation of the estimator.
#'
#' @param n_batch Number of samples per batch.
#'
#' @param n_max Maximum total number of samples (see stopping criterion).
#'
#' @param use_threads Number of chunks per batch. With
#'   \code{backend = "parallel"} also the number of cores or workers.
#'
#' @param backend Parallel backend, \code{"future"} or \code{"parallel"}. If
#'   \code{NULL} (default), \code{"parallel"} with one worker is used, i.e.
#'   the calculation runs sequentially on every system.
#'
#' @param sys_type \code{"serial"} or \code{"parallel"}. Only used for systems.
#'
#' @param dataRecord Logical. If \code{TRUE}, the results of each batch are
#'   returned in \code{data}. Single limit-state function only.
#'
#' @param beta_l Currently without effect.
#'
#' @param densityType Currently without effect. A normal sampling density in
#'   standard normal space is always used.
#'
#' @param dps Optional design points in physical space. If supplied, FORM is
#'   skipped for the components. Single limit-state function: a numeric vector
#'   with one value per variable (or one value for all). System: a list with
#'   one such vector per limit-state function. A numeric vector with one value
#'   per limit-state function is accepted if each function has exactly one
#'   variable.
#'
#' @param debug.level Verbosity level: 0 = silent, 1 = summary,
#'   2 = detailed output.
#'
#' @param seed Optional integer seed for reproducible results.
#'
#' @param adaptive_alpha Logical. Adaptive mixture weights are not
#'   implemented. Without effect for a single limit-state function. For serial
#'   multimodal systems the calculation stops with an error once
#'   \code{min_adapt_samples} is exceeded.
#'
#' @param alpha_update_rate Reserved for \code{adaptive_alpha}, currently
#'   without effect.
#'
#' @param adaptive_batch Logical. Systems only. If \code{TRUE}, the batch size
#'   is adapted from the current effective number of failure samples, so that
#'   \code{cov_user} is expected to be reached after \code{K_future} further
#'   batches.
#'
#' @param batch_control List for \code{adaptive_batch} with the elements
#'   \code{n_min} and \code{n_max} (bounds of the batch size) and
#'   \code{K_future} (number of remaining batches). Systems only.
#'
#' @param min_adapt_samples Reserved for \code{adaptive_alpha} (see there).
#'
#' @param alpha_min Reserved for \code{adaptive_alpha}, currently without
#'   effect.
#'
#' @param stability_mode \code{"robust"} or \code{"fast"}. Without influence
#'   on the result, all sums are accumulated in the log domain. For serial
#'   multimodal systems \code{"fast"} additionally returns the effective sample
#'   size per mixture component (\code{ESS_modes}).
#'
#' @return For a single limit-state function a list with
#' \itemize{
#'   \item \code{method} - \code{"MCIS_Single"}
#'   \item \code{beta} - reliability index \eqn{-\Phi^{-1}(\hat{P}_f)}
#'   \item \code{pf} - estimated failure probability
#'   \item \code{log_pf} - natural logarithm of \code{pf}
#'   \item \code{log_var} - natural logarithm of \code{var}
#'   \item \code{FORM_beta}, \code{FORM_pf} - FORM results (\code{NA} if
#'     \code{dps} was supplied)
#'   \item \code{design_points} - centre of the sampling density in standard
#'     normal space
#'   \item \code{var} - variance of the estimator
#'   \item \code{cov_mc} - coefficient of variation of the estimator
#'   \item \code{cov_user} - target coefficient of variation
#'   \item \code{ESS} - effective sample size
#'     \eqn{(\sum w_i)^2 / \sum w_i^2}
#'   \item \code{ESS_ratio} - \code{ESS / n_mc}
#'   \item \code{n_mc}, \code{n_max} - used and maximum number of samples
#'   \item \code{data} - results per batch (if \code{dataRecord = TRUE})
#'   \item \code{runtime} - runtime as returned by \code{proc.time()}
#' }
#' For a system a list with
#' \itemize{
#'   \item \code{method} - \code{"MCIS_System"}
#'   \item \code{sys_type} - \code{"serial"} or \code{"parallel"}
#'   \item \code{beta}, \code{pf}, \code{var}, \code{cov_mc} - as above
#'   \item \code{log_pf}, \code{log_var} - as above
#'   \item \code{RSE_est} - identical to \code{cov_mc}
#'   \item \code{ESS} - effective number of failure samples
#'     \eqn{(\sum I_i w_i)^2 / \sum (I_i w_i)^2}
#'   \item \code{ESS_modes} - effective sample size per mixture component
#'     (serial systems with \code{stability_mode = "fast"}, otherwise
#'     \code{NULL})
#'   \item \code{n_mc} - number of samples used
#' }
#'
#' @section Limit-State Function (LSF) Interface:
#'
#' The limit-state function follows the sign convention
#' \eqn{g(x) > 0} safe, \eqn{g(x) = 0} limit state and \eqn{g(x) < 0}
#' failure. It must use one of the following signatures:
#'
#' \describe{
#'   \item{\code{function(x)}}{
#'     A single numeric vector \code{x} containing all basic variables
#'     in the order defined in \code{lDistr}. This form also allows vectorized
#'     evaluation and is recommended.
#'   }
#'   \item{\code{function(Z, Fy, M, ...)}}{
#'     Explicitly named scalar arguments corresponding to the variable
#'     names defined in the probabilistic model.
#'   }
#' }
#'
#' The use of \code{function(...)} as the sole argument is not supported,
#' because it leads to inconsistent behaviour between FORM-based and Monte
#' Carlo-based reliability methods.
#'
#' @section Windows users:
#'
#' On Windows, \code{backend = "future"} with
#' \code{future::plan(multisession, workers = ...)} is recommended.
#' \code{backend = "parallel"} creates a new PSOCK cluster for every call and
#' issues a warning.
#'
#' @references
#' Ditlevsen, O., & Madsen, H. O. (1996).
#' \emph{Structural Reliability Methods}.
#' Wiley.
#'
#' Melchers, R. E., & Beck, A. T. (2018).
#' \emph{Structural Reliability Analysis and Prediction}.
#' Wiley.
#'
#' Spaethe, G. (1991).
#' \emph{Die Sicherheit tragender Baukonstruktionen}.
#' Springer.
#'
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#'
#' @import parallel
#' @import future
#' @import future.apply
#' @export
MC_IS <- function(lsf,
                  lDistr,
                  cov_user = 0.05,
                  n_batch = 5000,
                  n_max = 1e6,
                  use_threads = NULL,
                  backend = NULL,
                  sys_type = "serial",
                  dataRecord = TRUE,
                  beta_l = 100,
                  densityType = "norm",
                  dps = NULL,
                  debug.level = 0,
                  seed = NULL,
                  adaptive_alpha = FALSE,
                  alpha_update_rate = 0.1,
                  adaptive_batch = FALSE,
                  batch_control = list(),
                  min_adapt_samples = 10000,
                  alpha_min = 0.02,
                  stability_mode = c("robust", "fast")) {

  # Record whether the caller explicitly supplied use_threads.
  use_threads_was_supplied <- !missing(use_threads) &&
    !is.null(use_threads)

  stability_mode <- match.arg(stability_mode)


  ## -------------------------------------------------------------------------
  ## Backend selection with fallback and validation
  ## -------------------------------------------------------------------------
  valid_backends <- c("parallel", "future")

  # Keep NULL in the function signature so omitted arguments can be resolved here.
  if (is.null(backend)) {
    backend <- "parallel"
  }

  if (!is.character(backend) ||
    length(backend) != 1L ||
    is.na(backend) ||
    !backend %in% valid_backends) {
    stop(
      sprintf(
        "[TesiproV::MC_IS] Invalid backend '%s'. Valid options are: %s.",
        paste(backend, collapse = ", "),
        paste(valid_backends, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  ## -------------------------------------------------------------------------
  ## Argument checks 
  ## -------------------------------------------------------------------------
  if (!is.numeric(n_batch) || length(n_batch) != 1L ||
      !is.finite(n_batch) || n_batch < 1) {
    stop("`n_batch` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(n_max) || length(n_max) != 1L ||
      !is.finite(n_max) || n_max < 1) {
    stop("`n_max` must be a positive integer.", call. = FALSE)
  }

  

  # Validate an explicitly supplied worker count once.
  if (use_threads_was_supplied) {
    if (!is.numeric(use_threads) ||
      length(use_threads) != 1L ||
      !is.finite(use_threads) ||
      use_threads < 1 ||
      use_threads != floor(use_threads)) {
      stop(
        "`use_threads` must be a single positive integer when supplied.",
        call. = FALSE
      )
    }

    use_threads <- as.integer(use_threads)
  }

  if (backend == "parallel") {
    # Parallel backend defaults to one worker.
    use_threads_effective <- if (is.null(use_threads)) {
      1L
    } else {
      use_threads
    }
  } else {
    # The active Future plan determines the worker count.
    current_plan <- future::plan()
    plan_workers_raw <- future::nbrOfWorkers()

    if (!is.numeric(plan_workers_raw) ||
      length(plan_workers_raw) != 1L ||
      !is.finite(plan_workers_raw) ||
      plan_workers_raw < 1 ||
      plan_workers_raw != floor(plan_workers_raw)) {
      stop(
        paste0(
          "[TesiproV::MC_IS] Could not determine a finite positive integer ",
          "worker count from the active Future plan."
        ),
        call. = FALSE
      )
    }

    plan_workers <- as.integer(plan_workers_raw)
    use_threads_effective <- plan_workers

    if (use_threads_was_supplied) {
      warning(
        sprintf(
          paste0(
            "[TesiproV::MC_IS] backend = 'future' ignores 'use_threads = %d' ",
            "and uses the active Future plan with %d worker%s."
          ),
          use_threads,
          plan_workers,
          if (plan_workers == 1L) "" else "s"
        ),
        call. = FALSE
      )
    }

    if (inherits(current_plan, "sequential") &&
      plan_workers == 1L &&
      !isTRUE(.mc_is_state$plan_hint_shown)) {
      message(
        "[TesiproV::MC_IS] The active Future plan is sequential.\n",
        "Computation will run sequentially.\n",
        "To enable parallel execution, set a parallel Future plan, for example:\n",
        if (.Platform$OS.type == "windows") {
          "  future::plan(future::multisession, workers = 4)"
        } else {
          "  future::plan(future::multicore, workers = 4)"
        }
      )
      .mc_is_state$plan_hint_shown <- TRUE
    }
  }


  # -------------------------------------------------------------------------
  ## RNG initialization depending on backend
  ## -------------------------------------------------------------------------
  if (backend == "future") {
    if (is.null(seed)) {
      # No fixed seed provided -> create and set a random global seed
      if (debug.level >= 1) message("[TesiproV::MC_IS] No seed provided -- using random streams for 'future' backend.")

      rng_state <- init_rng_master(NULL, debug.level = debug.level) # <--- call the helper to draw and set a random seed!

      # The first integer element is sufficient as base seed for future.apply()
      seed_used <- NULL # must remain NULL for random seeding

      # --- Optional debug output -----------------------------------------------
      if (debug.level >= 1) {
        message("[TesiproV::MC_IS] Future backend running with automatic random substreams.")
        message(sprintf(
          "[TesiproV::MC_IS] Global RNG state head: %s",
          paste(.Random.seed[1:5], collapse = ", ")
        ))
      }
    } else {
      # Fixed deterministic seed -> reproducible run
      if (debug.level >= 1) message("[TesiproV::MC_IS] Fixed seed provided -- using reproducible streams for 'future' backend.")

      rng_state <- init_rng_master(seed, debug.level = debug.level)

      # The first integer element of the returned stream is sufficient as base seed
      # seed_used <- as.integer(rng_state[1])
      seed_used <- seed
    }
  } else if (backend == "parallel") {
    if (!is.null(seed)) {
      if (debug.level >= 1) message("[TesiproV::MC_IS] Fixed seed provided -- using reproducible streams for 'parallel' backend.")
      rng_state <- init_rng_master(seed, debug.level = debug.level)
      seed_used <- seed
    } else {
      if (debug.level >= 1) message("[TesiproV::MC_IS] No seed provided -- using random streams for 'parallel' backend.")
      rng_state <- init_rng_master(seed, debug.level = debug.level)
      seed_used <- NULL
    }
  }

  ## -------------------------------------------------------------------------
  ## -------------------------------------------------------------------------
  ## Windows fallback for threads
  ## -------------------------------------------------------------------------
  if (.Platform$OS.type == "windows") {
    if (debug.level >= 1) message("[MC_IS] Windows detected - PSOCK cluster will be used (if backend == 'parallel').")
  }
  # if (.Platform$OS.type == "windows") use_threads <- 1L

  ## -------------------------------------------------------------------------
  ## Cluster library paths setup
  ## -------------------------------------------------------------------------
  libPaths_local <- .libPaths()
  # assign("libPaths_local", libPaths_local, envir = .GlobalEnv)

  ## -------------------------------------------------------------------------
  ## Decide which mode to run: system or single limit state function
  ## -------------------------------------------------------------------------

  if (is.list(lsf)) {
    return(
      MC_IS_system(
        lsf = lsf, lDistr = lDistr,
        cov_user = cov_user, n_batch = n_batch, n_max = n_max,
        use_threads = use_threads_effective, sys_type = sys_type, dataRecord = dataRecord,
        beta_l = beta_l, densityType = densityType, dps = dps,
        debug.level = debug.level, streams = NULL, libPaths_local = libPaths_local,
        seed = seed_used, backend = backend, adaptive_alpha = adaptive_alpha,
        alpha_update_rate = alpha_update_rate,
        adaptive_batch = adaptive_batch,
        batch_control = batch_control,
        min_adapt_samples = min_adapt_samples,
        alpha_min = alpha_min,
        stability_mode = stability_mode
      )
    )
  } else if (is.function(lsf)) {
    return(
      MC_IS_single(
        lsf = lsf, lDistr = lDistr, cov_user = cov_user, n_batch = n_batch,
        n_max = n_max, use_threads = use_threads_effective, dataRecord = dataRecord,
        densityType = densityType, dps = dps,
        debug.level = debug.level, streams = NULL,
        libPaths_local = libPaths_local,
        seed = seed_used, backend = backend,
        adaptive_alpha = adaptive_alpha,
        alpha_update_rate = alpha_update_rate,
        min_adapt_samples = min_adapt_samples,
        alpha_min = alpha_min,
        stability_mode = stability_mode
      )
    )
  } else {
    stop("[TesiproV::MC_IS] Argument 'lsf' must be either a function (single LSF) or a list of functions (system).",
      call. = FALSE
    )
  }
}
