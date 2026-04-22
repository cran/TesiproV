#' @name MC_IS
#' @title Monte Carlo Simulation with Importance Sampling
#' @description
#' Performs Monte Carlo simulation with Importance Sampling (IS) to estimate
#' structural failure probabilities for single limit state functions or
#' systems of limit state functions.
#'
#' The method reduces the variance of crude Monte Carlo by sampling from a
#' shifted (or multimodal) density in standard normal space and correcting
#' via likelihood ratios. For systems, a multimodal sampling density based on
#' FORM design points is used (cf. Melchers & Beck, Chapter 5).
#'
#' The implementation supports:
#' \itemize{
#'   \item Single limit state functions (FORM-based IS shift)
#'   \item Serial and parallel systems (union / intersection events)
#'   \item Multimodal sampling densities for system reliability
#'   \item Optional adaptive mixture weights (adaptive \eqn{\alpha_i})
#'   \item Parallel execution (via \pkg{future} or \pkg{parallel})
#'   \item Deterministic reproducibility via fixed seeds
#'   \item Automatic detection of vectorized limit-state functions
#' }
#'
#' @details
#' \strong{Single limit state case:}
#' The importance sampling density is constructed in standard normal space
#' using the FORM design point \eqn{u^*}. Samples are generated from
#' \eqn{u \sim \mathcal{N}(u^*, I)} and transformed to physical space.
#'
#'
#' \strong{System case:}
#' For a system with limit states \eqn{g_i}, the failure event is
#' \eqn{F = \cup_i F_i} (serial system) or
#' \eqn{F = \cap_i F_i} (parallel system).
#'
#' The sampling density follows Melchers & Beck (Eq. 5.26):
#' \deqn{
#'   h(u) = \sum_{i=1}^m \alpha_i \, \phi(u - u_i^*)
#' }
#' where \eqn{u_i^*} are the FORM design points in standard normal space and
#' \eqn{\alpha_i} are mixture weights.
#'
#' By default, mixture weights are proportional to individual failure
#' probabilities:
#' \deqn{
#'   \alpha_i \propto P(F_i).
#' }
#'
#' If \code{adaptive_alpha = TRUE}, mixture weights are updated iteratively
#' based on weighted failure frequency estimates:
#' \deqn{
#'   \alpha_i^{new} =
#'   (1-\lambda)\alpha_i^{old} +
#'   \lambda \frac{\hat{P}(F_i)}{\sum_k \hat{P}(F_k)}
#' }
#' where \eqn{\lambda} is given by \code{alpha_update_rate}.
#'
#' The failure probability is estimated using a self-normalized IS estimator:
#' \deqn{
#'   \hat{P}_f =
#'   \frac{\sum I(u) w(u)}{\sum w(u)}
#' }
#' with likelihood ratio
#' \deqn{
#'   w(u) = \frac{\phi(u)}{h(u)}.
#' }
#'
#' Numerical stability is ensured via global log-weight stabilization.
#'
#'
#' \strong{Stopping Criterion}
#'
#' The simulation terminates when both:
#' \itemize{
#'   \item The estimated coefficient of variation (CoV) falls below \code{cov_user}, and
#'   \item The effective sample size (ESS) exceeds a minimum threshold.
#' }
#'
#' For parallel systems with very small failure probabilities, the ESS-based
#' safeguard prevents premature termination due to unstable variance estimates.
#'
#'
#' \strong{Numerical stability}
#'
#' The default \code{"robust"} mode performs global log-domain
#' accumulation of importance sampling weights using a
#' log-sum-exp formulation. This avoids overflow and underflow
#' effects in high-dimensional or rare-event settings.
#'
#' The alternative \code{"fast"} mode uses shifted exponentiation
#' at the batch level and is computationally equivalent for most
#' engineering reliability problems.
#'
#'
#' \strong{Automatic detection of vectorized limit-state functions}
#'
#' For performance reasons, the algorithm automatically checks whether the
#' supplied limit-state function (LSF) supports vectorized evaluation.
#'
#' A short internal test is performed before the simulation starts. If the LSF
#' accepts a matrix of input samples and returns a numeric vector of matching
#' length, it is evaluated in fully vectorized form:
#'
#' \deqn{g(X_{1:n}) \rightarrow \{g(x_1), \dots, g(x_n)\}}
#'
#' Otherwise, the algorithm falls back to row-wise evaluation using
#' \code{apply()}.
#'
#' This mechanism is fully automatic and backward compatible.
#' Users do not need to modify existing scalar LSF definitions.
#'
#' For computationally expensive LSFs (e.g., nonlinear models or surrogate FEM
#' models), vectorized evaluation can significantly improve performance.
#'
#'
#' \strong{Parallel execution and future plan handling}
#'
#' When \code{backend = "future"} is used, the function respects the currently
#' active \pkg{future} plan.
#'
#' If no parallel plan has been set by the user (i.e., the active plan is
#' sequential), a temporary plan is created internally using:
#'
#' \itemize{
#'   \item \code{multicore} on Unix-like systems,
#'   \item \code{multisession} on Windows.
#' }
#'
#' The original plan is automatically restored after completion of the
#' simulation.
#'
#' If the user has already defined a parallel strategy via
#' \code{future::plan()}, it will not be modified.
#'
#' This design ensures CRAN compliance while preserving full flexibility
#' for advanced users.
#'
#'
#' \strong{Reproducibility}
#'
#' Reproducible results across different numbers of cores are ensured through:
#'
#' \itemize{
#'   \item L'Ecuyer-CMRG random number streams,
#'   \item deterministic master sampling,
#'   \item worker-invariant parallel chunking,
#'   \item controlled seed propagation.
#' }
#'
#' If \code{seed} is supplied, results are reproducible regardless of the
#' number of worker threads.
#'
#'
#' @param lsf A single limit-state function \code{function(x)} or a list of
#'   such functions for system reliability. Must follow the LSF interface
#'   specification described below.
#'
#' @param lDistr List of marginal distribution objects (funlists with \code{$d},
#'   \code{$p}, \code{$q}, \code{$r}) corresponding to input variables.
#'
#' @param cov_user Target coefficient of variation (CoV) for the Monte Carlo
#'   estimator. Simulation stops once this threshold is reached.
#'
#' @param n_batch Number of samples generated per iteration.
#'
#' @param n_max Maximum total number of samples (upper stopping limit).
#'
#' @param use_threads Number of worker threads (for parallel execution).
#'
#' @param backend Parallel backend, either \code{"future"} or \code{"parallel"}.
#'
#' @param sys_type Character string, either \code{"serial"} or \code{"parallel"},
#'   specifying system configuration.
#'
#' @param dataRecord Logical; if \code{TRUE}, intermediate results per iteration
#'   are stored and returned.
#'
#' @param beta_l Optional threshold; limit states with \eqn{\beta > beta_l}
#'   may be excluded in system analysis.
#'
#' @param densityType Sampling density type (currently \code{"norm"} supported).
#'
#' @param dps Optional vector of design points in physical space; if supplied,
#'   FORM analysis is skipped.
#'
#' @param debug.level Integer verbosity level (0 = silent, 1 = summary, 2 = detailed).
#'
#' @param seed Optional integer seed for reproducible random numbers.
#'
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
#' @return MC_IS returns an object containing the following elements:
#' \itemize{
#'   \item \code{method} - method identifier
#'   \item \code{pf} - estimated failure probability
#'   \item \code{beta} - reliability index \eqn{\beta = -\Phi^{-1}(P_f)}
#'   \item \code{var} - variance estimate
#'   \item \code{cov_mc} - coefficient of variation
#'   \item \code{n_mc} - number of samples used
#'   \item \code{ESS} - effective sample size
#'   \item \code{ESS_modes} - effective sample size per mixture component (system case)
#'   \item \code{data} - optional iteration history (if \code{dataRecord = TRUE})
#' }
#'
#' @section Limit-State Function (LSF) Interface:
#'
#' The argument \code{lsf} must define a valid limit-state function
#' using one of the following supported signatures:
#'
#' \describe{
#'   \item{\code{function(x)}}{
#'     A single numeric vector \code{x} containing all basic variables
#'     in the order defined in \code{lDistr}.
#'
#'     This form is fully compatible with vectorized evaluation
#'     in Monte Carlo and importance sampling algorithms.
#'   }
#'
#'   \item{\code{function(Z, Fy, M, ...)}}{
#'     Explicitly named scalar arguments corresponding to the variable
#'     names defined in the probabilistic model.
#'   }
#' }
#'
#' The use of \code{function(...)} as the sole argument is not supported.
#' Such usage may cause inconsistent behaviour between FORM-based and
#' Monte Carlo-based reliability methods.
#'
#' The vector form \code{function(x)} is recommended for
#' best performance and consistency.
#'
#' @section Windows users:
#'
#' On Windows systems, it is recommended to use
#' \code{backend = "future"} for multi-core parallelization.
#' The \code{"parallel"} backend does not support forking on Windows.
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
                  use_threads = 4,
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
  stability_mode <- match.arg(stability_mode)

  ## -------------------------------------------------------------------------
  ## Backend selection with fallback and validation
  ## -------------------------------------------------------------------------
  valid_backends <- c("parallel", "future")

  ## message("Debug.level 0 %d", debug.level)

  if (is.null(backend)) {
    message("[TesiproV::MC_IS] No backend specified -- using 'future' by default.")
    backend <- "future"
  } else if (!(backend %in% valid_backends)) {
    stop(sprintf(
      "[TesiproV::MC_IS] Invalid backend '%s'. Valid options are: %s.",
      backend, paste(valid_backends, collapse = ", ")
    ), call. = FALSE)
  }

  ## Warning if parallel was chosen under Windows
  if (.Platform$OS.type == "windows" && backend == "parallel") {
    warning(
      "For Windows systems, backend = 'future' is recommended ",
      "for robust multi-core parallelization."
    )
  }

  ## Warning if no future plan was set by the user
  if (backend == "future") {
    current_plan <- future::plan()

    if (inherits(current_plan, "sequential")) {
      message(
        "[TesiproV::MC_IS] No active future plan detected.\n",
        "Computation will run sequentially.\n",
        "To enable parallel execution, set for example:\n",
        if (.Platform$OS.type == "windows") {
          "  future::plan(multisession, workers = 4)"
        } else {
          "  future::plan(multicore, workers = 4)"
        }
      )
    }
  }

  # ## Initialize Future plan
  # if (backend == "future") {
  #   current_plan <- future::plan()

  #   is_multicore <- inherits(current_plan, "multicore")
  #   is_multisession <- inherits(current_plan, "multisession")

  #   if (.Platform$OS.type == "windows") {
  #     if (!is_multisession) {
  #       future::plan(future::multisession, workers = use_threads)
  #     }
  #   } else {
  #     if (!is_multicore) {
  #       future::plan(future::multicore, workers = use_threads)
  #     }
  #   }
  # }

  # -------------------------------------------------------------------------
  ## RNG initialization depending on backend
  ## -------------------------------------------------------------------------
  if (backend == "future") {
    if (is.null(seed)) {
      # No fixed seed provided -> create and set a random global seed
      message("[TesiproV::MC_IS] No seed provided -- using random streams for 'future' backend.")

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
      message("[TesiproV::MC_IS] Fixed seed provided -- using reproducible streams for 'future' backend.")

      rng_state <- init_rng_master(seed, debug.level = debug.level)

      # The first integer element of the returned stream is sufficient as base seed
      # seed_used <- as.integer(rng_state[1])
      seed_used <- seed
    }
  } else if (backend == "parallel") {
    if (!is.null(seed)) {
      message("[TesiproV::MC_IS] Fixed seed provided -- using reproducible streams for 'parallel' backend.")
      rng_state <- init_rng_master(seed, debug.level = debug.level)
      seed_used <- rng_state
    } else {
      message("[TesiproV::MC_IS] No seed provided -- using random streams for 'parallel' backend.")
      rng_state <- init_rng_master(seed, debug.level = debug.level)
      seed_used <- NULL
    }
  }

  ## -------------------------------------------------------------------------
  ## -------------------------------------------------------------------------
  ## Windows fallback for threads
  ## -------------------------------------------------------------------------
  if (.Platform$OS.type == "windows") {
    message("[MC_IS] Windows detected - PSOCK cluster will be used (if backend == 'parallel').")
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
        use_threads = use_threads, sys_type = sys_type, dataRecord = dataRecord,
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
        n_max = n_max, use_threads = use_threads, dataRecord = dataRecord,
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
