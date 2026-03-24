#' @name MC_CRUDE
#' @title Crude Monte Carlo Simulation
#' @description
#' Crude Monte Carlo simulation (MC_CRUDE) method for the estimation of failure
#' probabilities in structural reliability analysis.
#'
#' The crude Monte Carlo method estimates the probability of failure
#' by direct sampling of the basic random variables and evaluation
#' of the limit-state function.
#'
#' The failure probability is estimated as
#' \deqn{
#'   \hat{p}_f = \frac{1}{N}
#'   \sum_{i=1}^{N}
#'   \mathbf{1}_{\{g(\mathbf{X}_i) \le 0\}}.
#' }
#'
#' Although conceptually simple and unbiased, the method becomes
#' computationally inefficient for very small failure probabilities,
#' as a large number of samples is required to obtain stable estimates.
#'
#' @param lsf Limit-state function. Must follow the LSF interface
#'   specification described in the "Limit-State Function (LSF) Interface"
#'   section below.
#'
#' @param lDistr List of distribution objects as returned by \code{PROB_BASEVAR$getlDistr()}.
#'
#' @param cov_user Target coefficient of variation to be achieved.
#'
#' @param n_batch Batch size per iteration (used for parallel computation).
#'
#' @param n_max Maximum number of Monte Carlo samples (stopping criterion).
#'
#' @param use_threads Number of parallel threads. Set to 1 for single-core
#'   execution. Parallel execution is not supported on Windows.
#'
#' @param backend Parallel backend, either \code{"future"} (default) or \code{"parallel"}.
#'
#' @param dataRecord Logical; if \code{TRUE}, intermediate results are recorded.
#'
#' @param debug.level Integer controlling verbosity (0 = silent, 2 = detailed output).
#'
#' @param seed Optional integer value for reproducible random numbers.
#'
#' @return MC_CRUDE returns an object containing the following elements:
#' \itemize{
#'   \item \code{method}: Character string identifying the method ("MCC").
#'   \item \code{beta}: Estimated reliability index \eqn{\beta = -\Phi^{-1}(p_f)}.
#'   \item \code{pf}: Estimated probability of failure.
#'   \item \code{var}: Estimated variance of the Monte Carlo estimator.
#'   \item \code{cov_mc}: Estimated coefficient of variation of \code{pf}.
#'   \item \code{cov_user}: Target coefficient of variation specified by the user.
#'   \item \code{n_mc}: Total number of Monte Carlo samples generated.
#'   \item \code{n_max}: Maximum allowed number of samples.
#'   \item \code{n_batch}: Batch size per iteration.
#'   \item \code{n_threads}: Number of threads used.
#'   \item \code{runtime}: CPU time returned by \code{proc.time()}.
#'   \item \code{data}: Optional data frame containing intermediate results
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
#'     in the order defined in \code{lDistr}.
#'
#'     Example:
#'     \preformatted{
#'     lsf <- function(x) {
#'       x[1] - sum((x[-1]^2) / seq_along(x)[-1])
#'     }
#'     }
#'   }
#'
#'   \item{\code{function(Z, Fy, M, ...)}}{
#'     Explicitly named scalar arguments corresponding to the variable
#'     names defined in the probabilistic model.
#'
#'     Example:
#'     \preformatted{
#'     lsf <- function(Z, Fy, M) {
#'       Z * Fy - M
#'     }
#'     }
#'   }
#' }
#'
#' The following form is \strong{not supported}:
#'
#' \preformatted{
#' function(...)
#' }
#'
#' Using \code{...} as the sole argument may lead to inconsistent
#' behaviour across different reliability algorithms and is therefore
#' prohibited.
#'
#' For maximum robustness and cross-method compatibility,
#' the vector form \code{function(x)} is recommended.
#'
#' @section Windows users:
#'
#' On Windows systems, it is recommended to use
#' \code{backend = "future"} for multi-core parallelization.
#' The \code{"parallel"} backend does not support forking on Windows.
#'
#' @references
#' Spaethe, G. (1991).
#' \emph{Die Sicherheit tragender Baukonstruktionen}.
#' Springer.'
#'
#' @import parallel
#'
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#'
#' @export
#'

MC_CRUDE <- function(lsf,
                     lDistr,
                     cov_user = 0.05,
                     n_batch = 10000,
                     n_max = 1e7,
                     use_threads = parallel::detectCores(),
                     backend = c("future", "parallel"),
                     dataRecord = TRUE,
                     debug.level = 0,
                     seed = NULL) {
  backend <- match.arg(backend)

  # Warning if parallel was chosen under Windows
  if (.Platform$OS.type == "windows" && backend == "parallel") {
    warning(
      "For Windows systems, backend = 'future' is recommended ",
      "for robust multi-core parallelization."
    )
  }

  ## Initialize Future plan
  if (backend == "future") {
    current_plan <- future::plan()

    is_multicore <- inherits(current_plan, "multicore")
    is_multisession <- inherits(current_plan, "multisession")

    if (.Platform$OS.type == "windows") {
      if (!is_multisession) {
        future::plan(future::multisession, workers = use_threads)
      }
    } else {
      if (!is_multicore) {
        future::plan(future::multicore, workers = use_threads)
      }
    }
  }

  debug.TAG <- "MC_Crude"
  debug.print(debug.level, debug.TAG, c(TRUE), msg = "Crude Monte-Carlo Simulation started...")

  tic <- proc.time()

  pf <- 1
  I_n <- 0L
  # pf_i  <- 0
  # var_i <- 0
  n_sim <- 0L
  cov_mc <- +Inf
  var_mc <- NA_real_
  k <- 0L

  n_vars <- length(lDistr)

  v <- matrix(nrow = n_batch, ncol = n_vars)

  # Data record
  # If dataRecord = TRUE, record the following data (frequency saveEvery)

  if (dataRecord) {
    # create empty DataFraame
    df <- data.frame(
      n_sim = numeric(),
      pf = numeric(),
      var = numeric(),
      cov = numeric(),
      time.user = numeric(),
      time.sys = numeric(),
      time.elapsed = numeric(),
      time.user.child = numeric(),
      time.sys.child = numeric()
    )
  } else {
    df <- data.frame()
  }

  # Create empty list
  records <- list()

  # Definition of update frequency of dataRecord
  saveEvery <- max(10L, (use_threads * n_batch) / 250L)

  # Generate random numbers
  # Reproducible streams for parallel random numbers
  RNGkind(kind = "L'Ecuyer-CMRG")

  if (!is.null(seed)) {
    set.seed(seed)
    message(sprintf("Fixed RNG seed set to %d", seed))
  }

  # # mc_local is used for parallelisation purpose
  mc_local <- function(x) {
    v_local <- matrix(nrow = n_batch, ncol = n_vars)

    for (i in seq_len(n_vars)) {
      v_local[, i] <- lDistr[[i]][[1]]$r(n_batch)
    }
    if (is.function(lsf)) {
      I <- sum(apply(v_local, 1, lsf) < 0)
    } else if (!is.null(lsf$func) && is.function(lsf$func)) {
      I <- sum(apply(v_local, 1, lsf$func) < 0)
    } else {
      stop("Parameter 'lsf' must be a function or SYS_LSF object with field $func.")
    }

    return(I)
  }

  while (TRUE) {
    # # create Realisations
    # if (Sys.info()[[1]]=="Windows"){ #Windows is not able to fork - parallelisation under windows os not efficient
    #   I_n <- I_n + sum(unlist(parallel::mclapply(seq(1,use_threads),mc_local, mc.cores=1)))
    # }else{ # parallelisation for unix based platforms (macOS, Linux etc.)
    #   I_n <- I_n + sum(unlist(parallel::mclapply(seq(1,use_threads),mc_local, mc.set.seed = TRUE, mc.cores=use_threads)))
    #   stats::runif(1)
    # }
    # n_sim <- n_sim + use_threads*n_batch

    # k <- k+1
    # For only one thread, mclapply is slower than lapply
    if (use_threads == 1) {
      I_list <- lapply(seq_len(use_threads), mc_local)
    } else {
      if (backend == "parallel") {
        I_list <- parallel::mclapply(
          seq_len(use_threads),
          mc_local,
          mc.set.seed = TRUE,
          mc.cores = use_threads
        )
      } else if (backend == "future") {
        I_list <- future.apply::future_lapply(
          seq_len(use_threads),
          mc_local,
          future.seed = TRUE
        )
      }
    }

    I_n <- I_n + sum(unlist(I_list))
    n_sim <- n_sim + use_threads * n_batch

    k <- k + 1

    # time measurement
    p_stamp <- proc.time() - tic

    # Calculate stochastics
    debug.print(debug.level, debug.TAG, I_n, "I_n: ")

    if (I_n > 0) {
      pf <- I_n / n_sim
      # add more robust equation for the CoV
      # var_mc <- 1/(n_sim-1)*((1/n_sim)*I_n - pf^2)
      var_mc <- pf * (1 - pf) / n_sim
      # direct calculation of CoV is faster
      cov_mc <- sqrt((1 - pf) / (pf * n_sim)) # sqrt(var_mc)/pf

      # additional information if option debug = 2
      info.print(debug.TAG, debug.level, c("I_n", "pf", "cov", "nsim"), c(I_n, pf, cov_mc, n_sim))
      debug.print(debug.level, debug.TAG, pf, "Pf: ")
      debug.print(debug.level, debug.TAG, var_mc, "var_mc: ")
      debug.print(debug.level, debug.TAG, cov_mc, "CoV_mc: ")
      debug.print(debug.level, debug.TAG, n_sim, "n_sim: ")
      debug.print(debug.level, debug.TAG, n_max, "n_max: ")

      if (cov_mc < cov_user || n_sim > n_max) break # break loop, if target CoV is reached
    } else {
      info.print(debug.TAG, debug.level, c("I_n", "pf", "cov", "nsim"), c(I_n, "NA", "NA", n_sim))

      if (n_sim > n_max) break # break loop, if n_sim > n_max
    }


    # Save data only every 10th step and use a list to increase performance
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

    # # Allocated RAM
    #   if (debug.level > 0 && k %% saveEvery == 0) {
    #   mem_now <- pryr::mem_used()
    #   msg <- sprintf("[%s] [Iter %d] RAM: %.2f MB\n",
    #                  format(Sys.time(), "%H:%M:%S"),
    #                  n_sim,
    #                  mem_now / (1024^2))
    #   cat(msg, file="MC_CRUDE_log.txt", append=TRUE)
    # }
  } # end of while loop


  # Convert list to DataFrame
  if (dataRecord && length(records) > 0) {
    df <- as.data.frame(do.call(rbind, lapply(records, as.data.frame)))
  } else {
    df <- data.frame()
  }

  cat("\n")
  duration <- proc.time() - tic

  # calculate reliabilty index beta
  beta_val <- ifelse(pf > 0 && !is.na(pf), -stats::qnorm(pf), Inf)

  output <- list(
    method = "MC_CRUDE",
    beta = beta_val, #-stats::qnorm(pf),
    pf = pf,
    var = var_mc,
    cov_mc = cov_mc,
    cov_user = cov_user,
    n_mc = n_sim,
    n_max = n_max,
    n_batch = n_batch,
    n_threads = use_threads,
    runtime = duration[1:5],
    data = df
  )

  if (dataRecord) {
    output$data <- df
  } else {
    output$data <- data.frame()
  }

  debug.print(debug.level, debug.TAG, c(duration), msg = "Crude Monte-Carlo Simulation finished in [s]:  ")

  return(output)
}
