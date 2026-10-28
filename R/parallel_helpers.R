# File: R/parallel_helpers.R

#' @title Create a parallel cluster (fork on *nix, PSOCK on Windows)
#' @description Returns a `parallel` cluster object that already has the
#'   RNG streams (L'Ecuyer-CMRG) and required packages loaded.
#' @param use_threads Integer \>=1 - number of workers.
#' @param seed Optional integer or L'Ecuyer-CMRG stream (length 7) for reproducibility.
#' @return A cluster object (or NULL when the backend is "future").
#' @keywords internal
make_parallel_cluster <- function(use_threads, seed = NULL) {
    ## -------------------------------------------------
    ## 1) Determine cluster type (fork vs. PSOCK)
    ## -------------------------------------------------
    cl_type <- if (.Platform$OS.type == "windows") "PSOCK" else "FORK"

    ## -------------------------------------------------
    ## 2) Create the cluster
    ## -------------------------------------------------
    cl <- parallel::makeCluster(use_threads, type = cl_type)

    ## -------------------------------------------------
    ## 3) Initialise RNG streams
    ## -------------------------------------------------
    if (!is.null(seed)) {
        ## User supplied seed (integer or full CMRG stream)
        parallel::clusterSetRNGStream(cl, iseed = seed)
    } else {
        ## Random sub-streams - reproducible later via future.seed = TRUE
        parallel::clusterSetRNGStream(cl)
    }

    ## -------------------------------------------------
    ## 4) Load required packages on the workers
    ## -------------------------------------------------
    parallel::clusterEvalQ(cl, {
        suppressMessages({
            if (!requireNamespace("TesiproV", quietly = TRUE)) {
                stop("Package 'TesiproV' not found on worker.", call. = FALSE)
            }
            ## Add additional packages that LS-functions may need, e.g.:
            ## library(evd); library(EnvStats); ...
            invisible(NULL)
        })
    })

    return(cl)
}


#' @title Deterministic RNG Block Manager (L'Ecuyer-CMRG)
#' @description
#' Provides deterministic generation of random number blocks that are
#' independent of the number of worker threads.
#'
#' The RNG state is advanced only by the requested block size.
#'
#' @param seed Optional integer for reproducibility.
#' @return An RNG manager object (environment).
#' @keywords internal
create_rng_manager <- function(seed = NULL) {
    RNGkind("L'Ecuyer-CMRG")

    if (!is.null(seed)) {
        set.seed(seed)
    }

    rng_env <- new.env(parent = emptyenv())
    rng_env$current_state <- .Random.seed

    rng_env$generate_norm_block <- function(n) {
        assign(".Random.seed", rng_env$current_state, envir = .GlobalEnv)
        values <- rnorm(n)
        rng_env$current_state <- .Random.seed
        return(values)
    }

    return(rng_env)
}


# File: R/parallel_helpers.R
#' @title Run a parallel computation using either "parallel" or "future"
#' @description Small wrapper that hides the branching between the two back-ends.
#'   It always returns a **named list** whose names are identical to the names
#'   of the input split list. This guarantees that `unlist()` produces the same
#'   named vector for both back-ends, which is required by the unit tests.
#' @param X List (or vector) of inputs that will be split among the workers.
#' @param FUN Function that will be executed on each chunk.
#' @param backend Either "parallel" or "future".
#' @param use_threads Number of workers (relevant for both back-ends).
#' @param cl Optional cluster object created with \code{parallel::makeCluster()}.
#'   Only required when using the \code{"parallel"} backend on Windows, where
#'   \code{parallel::parLapply()} is used. Ignored on Unix-like systems where
#'   \code{parallel::mclapply()} is available.
#' @return List of results (identical to the output of `parLapply` /
#'   `future_lapply`). The list is **named** - the names are taken from the
#'   split index (`split_idx`).
#' @keywords internal
#'
run_parallel <- function(X,
                         FUN,
                         backend = c("parallel", "future"),
                         use_threads = NULL,
                         cl = NULL,
                         ...) {
    backend <- match.arg(backend)

    if (backend == "parallel") {
        if (use_threads <= 1L) {
            return(lapply(X, FUN, ...))
        }

        ## ---- Linux / macOS: use fork directly ----
        if (.Platform$OS.type != "windows") {
            return(parallel::mclapply(
                X,
                FUN,
                ...,
                mc.cores = use_threads,
                mc.preschedule = TRUE
            ))
        }

        ## ---- Windows: use PSOCK cluster ----
        if (is.null(cl)) {
            stop("Cluster object 'cl' must be provided for Windows parallel backend.",
                call. = FALSE
            )
        }

        return(parallel::parLapply(cl, X, FUN, ...))
    }

    ## ---- future backend ----
    return(future.apply::future_lapply(
        X = X,
        FUN = FUN,
        ...,
        future.seed = FALSE, # deterministic master RNG
        future.globals = FALSE # prevent unnecessary global export
    ))
}


#' Resolve the worker count for a parallel backend
#'
#' Resolves the effective number of workers according to the selected
#' parallel backend.
#'
#' For \code{backend = "parallel"}, the worker count is taken from
#' \code{use_threads}. If it is \code{NULL}, \code{default_parallel}
#' is used.
#'
#' For \code{backend = "future"}, the active \pkg{future} plan determines
#' the worker count. The function does not change the active plan. If
#' \code{use_threads} is explicitly supplied, a warning explains that the
#' value does not control the Future backend; when it differs from the
#' plan's worker count, the warning reports both values.
#'
#' @param backend Character scalar. Backend name: \code{"parallel"} or
#'   \code{"future"}.
#' @param use_threads Optional positive integer specifying the number of
#'   workers for the \code{"parallel"} backend. For the \code{"future"}
#'   backend, the active Future plan determines the worker count; supplying
#'   this argument triggers a warning.
#' @param default_parallel Positive integer used as the worker-count
#'   default when \code{backend = "parallel"} and \code{use_threads} is
#'   \code{NULL}.
#' @param call Character scalar used to identify the calling function in
#'   warning and error messages.
#'
#' @return A positive integer giving the effective worker count.
#'
#' @details
#' A sequential Future plan corresponds to one effective worker. The
#' function reads the active plan's worker count using the Future package
#' and leaves the plan unchanged.
#'
#' @keywords internal
#' @noRd
.resolve_backend_workers <- function(backend,
                                     use_threads,
                                     default_parallel,
                                     call = "TesiproV") {
    backend <- match.arg(backend, c("parallel", "future"))

    if (backend == "parallel") {
        workers <- use_threads

        if (is.null(workers)) {
            workers <- default_parallel
        }

        if (!is.numeric(workers) ||
            length(workers) != 1L ||
            !is.finite(workers) ||
            workers < 1 ||
            workers != floor(workers)) {
          stop(
            "`use_threads` must be a single positive integer when supplied.",
            call. = FALSE
          )
        }

        return(as.integer(workers))
    }

    # For the future backend, the active Future plan is authoritative.
    plan <- future::plan()

    plan_workers <- tryCatch(
        future::nbrOfWorkers(plan),
        error = function(e) NA_integer_
    )

    if (length(plan_workers) != 1L ||
        is.na(plan_workers) ||
        !is.finite(plan_workers) ||
        plan_workers < 1) {
        stop(
            sprintf(
                "%s could not determine the worker count from the active future plan.",
                call
            ),
            call. = FALSE
        )
    }

    plan_workers <- as.integer(plan_workers)

    # Warn only if the caller explicitly supplied use_threads.
    supplied_use_threads <- !missing(use_threads)

    if (supplied_use_threads && !is.null(use_threads)) {
      # Warn whenever the caller supplied use_threads, as in MC_IS.
      if (!missing(use_threads) && !is.null(use_threads)) {
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

        warning(
          sprintf(
            paste0(
              "%s backend = 'future' ignores 'use_threads = %d' ",
              "and uses the active Future plan with %d worker%s."
            ),
            call,
            as.integer(use_threads),
            plan_workers,
            if (plan_workers == 1L) "" else "s"
          ),
          call. = FALSE
        )
      }
    }

    plan_workers
}
