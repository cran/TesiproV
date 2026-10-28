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


# Package-level state (e.g. to show a hint only once per R session)
.mc_is_state <- new.env(parent = emptyenv())

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
#' @keywords internal
init_rng_master <- function(seed = NULL, debug.level = 0) {
  RNGkind("L'Ecuyer-CMRG")

  if (!is.null(seed)) {
    # fixed user seed -> reproducible
    set.seed(as.integer(seed))
    if (debug.level >= 1L) message(sprintf("[init_rng_master] User seed set to %d", as.integer(seed)))
  } else {
    base_seed <- sample.int(.Machine$integer.max, 1L)
    set.seed(base_seed)
    if (debug.level >= 1L) message(sprintf("[init_rng_master] Random RNG seed drawn: %d", base_seed))
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
  ## No new seed needs to be generated - the same base seed that
  ## `init_rng_master()` has already set is used.
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
  rec$cov[idx] <- if (isTRUE(pf_mc != 0)) sqrt(var_mc) / pf_mc else NA_real_
  p_stamp <- proc.time() - tic
  rec$time$user[idx] <- p_stamp[1]
  rec$time$sys[idx] <- p_stamp[2]
  rec$time$elapsed[idx] <- p_stamp[3]
  rec$time$child_user[idx] <- p_stamp[4]
  rec$time$child_sys[idx] <- p_stamp[5]
  rec
}


#' @title Convert recorder lists into tidy data frame after completion
#' @description
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
    ## seed is a complete CMRG stream (length 7)
    parallel::clusterSetRNGStream(cl, iseed = seed)
    if (debug.level >= 1) {
      message("[make_ready_cluster] Fixed CMRG stream set on workers (backend parallel)")
    }
  } else {
    ## generate random sub-streams
    parallel::clusterSetRNGStream(cl) # automatically creates new streams
    if (debug.level >= 1) {
      message("[make_ready_cluster] Random sub-streams generated (backend parallel)")
    }
  }

  #   ## 2) Export of all required objects
  #   all_export <- unique(c(export_objs, "libPaths_local"))
  #   parallel::clusterExport(cl,
  #     varlist = all_export,
  #     envir = parent.frame()
  #   )

  #   ## 3) Load packages on the workers
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

    for (pkg in c("TesiproV", "evd")) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' not found on worker.", pkg), call. = FALSE)
      }
    }

    NULL
  })

  # Remove automatic stop handler (cluster valid)
  on.exit(NULL)

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
