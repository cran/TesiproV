#' @name MC_CRUDE
#' @title Crude Monte Carlo Simulation
#' @description
#' Crude Monte Carlo simulation (MC_CRUDE) method for the estimation of failure
#' probabilities in structural reliability analysis.
#'
#' The failure probability is estimated as
#' \deqn{
#'   \hat{p}_f = \frac{1}{N}
#'   \sum_{i=1}^{N}
#'   \mathbf{1}_{\{g(\mathbf{X}_i) < 0\}}.
#' }
#'
#' Failure is defined as \eqn{g(x) < 0}. Samples with \eqn{g(x) = 0} are
#' counted as safe.
#'
#' Although conceptually simple and unbiased, the method becomes
#' computationally inefficient for very small failure probabilities,
#' as a large number of samples is required to obtain stable estimates.
#'
#' @details
#' \strong{Stopping criterion}
#'
#' The simulation stops when the coefficient of variation falls below
#' \code{cov_user} or when \code{n_max} samples are reached. The last batch is
#' reduced so that the total does not exceed \code{n_max}.
#'
#' \strong{Parallel execution}
#'
#' The samples of each batch are generated in the master process and split
#' into chunks, which are evaluated in parallel. The workers do not generate
#' random numbers.
#' \itemize{
#'   \item \code{backend = "future"}: the chunks are evaluated with
#'     \code{future.apply::future_lapply()} under the active \pkg{future}
#'     plan. The function never changes the plan. The number of workers is
#'     determined by the plan, \code{use_threads} has no effect and issues a
#'     warning if supplied.
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
#' L'Ecuyer-CMRG stream in the master process. The result therefore does not
#' depend on the number of workers. The function sets
#' \code{RNGkind("L'Ecuyer-CMRG")} and changes \code{.Random.seed} in the
#' global environment.
#'
#' @param lsf Limit-state function. Must follow the LSF interface
#'   specification described in the "Limit-State Function (LSF) Interface"
#'   section below.
#' @param lDistr List of distribution objects as returned by
#'   \code{PROB_BASEVAR$getlDistr()}.
#' @param cov_user Target coefficient of variation to be achieved.
#' @param n_batch Number of samples per batch.
#' @param n_max Maximum number of Monte Carlo samples.
#' @param use_threads Number of workers for \code{backend = "parallel"}
#'   (default 1). Without effect for \code{backend = "future"}.
#' @param backend Parallel backend, \code{"future"} or \code{"parallel"}. If
#'   \code{NULL} (default), \code{"parallel"} with one worker is used, i.e.
#'   the calculation runs sequentially on every system.
#' @param dataRecord Logical; if \code{TRUE}, intermediate results are recorded.
#' @param debug.level Integer controlling verbosity (0 = silent, 2 = detailed
#'   output).
#' @param seed Optional integer value for reproducible random numbers.
#'
#' @return MC_CRUDE returns a list with the following elements:
#' \itemize{
#'   \item \code{method}: Character string "MC_CRUDE".
#'   \item \code{beta}: Estimated reliability index \eqn{\beta = -\Phi^{-1}(p_f)}.
#'   \item \code{pf}: Estimated probability of failure.
#'   \item \code{var}: Estimated variance of the Monte Carlo estimator.
#'   \item \code{cov_mc}: Estimated coefficient of variation of \code{pf}.
#'   \item \code{cov_user}: Target coefficient of variation.
#'   \item \code{n_mc}: Total number of Monte Carlo samples generated.
#'   \item \code{n_max}: Maximum allowed number of samples.
#'   \item \code{n_batch}: Batch size.
#'   \item \code{n_threads}: Number of workers used.
#'   \item \code{runtime}: CPU time returned by \code{proc.time()}.
#'   \item \code{data}: Data frame with intermediate results
#'         (only if \code{dataRecord = TRUE}).
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
#'     in the order defined in \code{lDistr}. Recommended.
#'   }
#'   \item{\code{function(Z, Fy, M, ...)}}{
#'     Explicitly named scalar arguments corresponding to the variable
#'     names defined in the probabilistic model.
#'   }
#' }
#'
#' The use of \code{function(...)} as the sole argument is not supported.
#'
#' @section Windows users:
#'
#' On Windows, \code{backend = "future"} with
#' \code{future::plan(future::multisession, workers = ...)} is recommended.
#' \code{backend = "parallel"} creates a new PSOCK cluster for every call.
#'
#' @references
#' Spaethe, G. (1991).
#' \emph{Die Sicherheit tragender Baukonstruktionen}.
#' Springer.
#'
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#'
#' @export
MC_CRUDE <- function(lsf,
                     lDistr,
                     cov_user = 0.05,
                     n_batch = 10000,
                     n_max = 1e7,
                     use_threads = NULL,
                     backend = NULL,
                     dataRecord = TRUE,
                     debug.level = 0,
                     seed = NULL) {
  # without a backend the calculation runs sequentially on every system
  if (is.null(backend)) backend <- "parallel"
  backend <- match.arg(backend, c("future", "parallel"))

  if (!is.numeric(n_batch) || length(n_batch) != 1L ||
      !is.finite(n_batch) || n_batch < 1) {
    stop("`n_batch` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(n_max) || length(n_max) != 1L ||
      !is.finite(n_max) || n_max < 1) {
    stop("`n_max` must be a positive integer.", call. = FALSE)
  }

  # ceck if lDistr is properly defined
  .check_lDistr(lDistr, "MC_CRUDE")

  # worker count: plan for "future", use_threads for "parallel"
  n_threads <- .resolve_backend_workers(
    backend = backend,
    use_threads = use_threads,
    default_parallel = 1L,
    call = "[TesiproV::MC_CRUDE]"
  )

  # normalize the limit-state function
  if (is.function(lsf)) {
    lsf_fun <- lsf
  } else if (!is.null(lsf$func) && is.function(lsf$func)) {
    lsf_fun <- lsf$func
  } else {
    stop("Parameter 'lsf' must be a function or SYS_LSF object with field $func.",
         call. = FALSE)
  }

  debug.TAG <- "MC_Crude"
  debug.print(debug.level, debug.TAG, c(TRUE), msg = "Crude Monte-Carlo Simulation started...")

  # PSOCK cluster for the Windows parallel backend. Created before the
  # RNG initialisation, clusterSetRNGStream() restores the master seed.
  cl <- NULL
  if (backend == "parallel" && n_threads > 1L &&
      .Platform$OS.type == "windows") {
    cl <- make_parallel_cluster(use_threads = n_threads, seed = NULL)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }

  # random numbers are generated in the master process only
  RNGkind(kind = "L'Ecuyer-CMRG")
  if (!is.null(seed)) {
    set.seed(seed)
    if (debug.level >= 1) message(sprintf("Fixed RNG seed set to %d", as.integer(seed)))
  }

  tic <- proc.time()

  pf <- 0
  I_n <- 0
  n_sim <- 0
  cov_mc <- Inf
  var_mc <- NA_real_
  k <- 0L

  n_vars <- length(lDistr)

  records <- list()
  saveEvery <- max(10L, as.integer(n_batch %/% 250))

  while (n_sim < n_max) {
    # last batch is reduced so that n_max is not exceeded
    total_batch <- as.integer(min(floor(n_batch), floor(n_max - n_sim)))
    if (total_batch < 1L) break

    v <- matrix(NA_real_, nrow = total_batch, ncol = n_vars)
    for (i in seq_len(n_vars)) {
      v[, i] <- lDistr[[i]][[1]]$r(total_batch)
    }

    if (n_threads == 1L) {
      chunk_indices <- list(seq_len(total_batch))
    } else {
      chunk_indices <- split(
        seq_len(total_batch),
        cut(seq_len(total_batch), breaks = n_threads, labels = FALSE)
      )
    }
    v_chunks <- lapply(chunk_indices, function(idx) v[idx, , drop = FALSE])

    I_list <- run_parallel(
      X = v_chunks,
      FUN = .mc_crude_worker,
      backend = backend,
      use_threads = n_threads,
      cl = cl,
      lsf = lsf_fun
    )

    I_n <- I_n + sum(unlist(I_list))
    n_sim <- n_sim + total_batch
    k <- k + 1L

    p_stamp <- proc.time() - tic

    debug.print(debug.level, debug.TAG, I_n, "I_n: ")

    if (I_n > 0) {
      pf <- I_n / n_sim
      var_mc <- pf * (1 - pf) / n_sim
      cov_mc <- sqrt((1 - pf) / (pf * n_sim))

      info.print(debug.TAG, debug.level, c("I_n", "pf", "cov", "nsim"), c(I_n, pf, cov_mc, n_sim))
      debug.print(debug.level, debug.TAG, pf, "Pf: ")
      debug.print(debug.level, debug.TAG, var_mc, "var_mc: ")
      debug.print(debug.level, debug.TAG, cov_mc, "CoV_mc: ")
      debug.print(debug.level, debug.TAG, n_sim, "n_sim: ")
      debug.print(debug.level, debug.TAG, n_max, "n_max: ")

      if (cov_mc < cov_user) break
    } else {
      info.print(debug.TAG, debug.level, c("I_n", "pf", "cov", "nsim"), c(I_n, "NA", "NA", n_sim))
    }

    if (dataRecord && k %% saveEvery == 0) {
      records[[length(records) + 1]] <- list(
        n_sim = n_sim,
        pf = pf,
        var = var_mc,
        cov = cov_mc,
        time.user = p_stamp[1],
        time.sys = p_stamp[2],
        time.elapsed = p_stamp[3],
        time.user.child = p_stamp[4],
        time.sys.child = p_stamp[5]
      )
    }
  }

  if (dataRecord && length(records) > 0) {
    df <- as.data.frame(do.call(rbind, lapply(records, as.data.frame)))
  } else {
    df <- data.frame()
  }

  duration <- proc.time() - tic

  beta_val <- if (pf > 0 && !is.na(pf)) -stats::qnorm(pf) else Inf

  output <- list(
    method = "MC_CRUDE",
    beta = beta_val,
    pf = pf,
    var = var_mc,
    cov_mc = cov_mc,
    cov_user = cov_user,
    n_mc = n_sim,
    n_max = n_max,
    n_batch = n_batch,
    n_threads = n_threads,
    runtime = duration[1:5],
    data = df
  )

  debug.print(debug.level, debug.TAG, c(duration), msg = "Crude Monte-Carlo Simulation finished in [s]:  ")

  return(output)
}


#' Evaluate one chunk in crude Monte Carlo simulation
#'
#' Counts the samples of a chunk that lie in the failure domain. Defined at
#' package level so that the worker does not capture the environment of
#' MC_CRUDE().
#'
#' @param v_chunk Numeric matrix, one sample per row, physical space.
#' @param lsf Limit-state function.
#' @return Number of samples with g(x) < 0.
#' @keywords internal
#' @noRd
.mc_crude_worker <- function(v_chunk, lsf) {
  sum(apply(v_chunk, 1, lsf) < 0)
}
