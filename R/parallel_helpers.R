# File: R/parallel_helpers.R

#' @title Create a parallel cluster (fork on *nix, PSOCK on Windows)
#' @description Returns a `parallel` cluster object that already has the
#'   RNG streams (L'Ecuyer‑CMRG) and required packages loaded.
#' @param use_threads Integer \>=1 – number of workers.
#' @param seed Optional integer or L'Ecuyer‑CMRG stream (length 7) for reproducibility.
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
        ## Random sub‑streams – reproducible later via future.seed = TRUE
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
            ## Add additional packages that LS‑functions may need, e.g.:
            ## library(evd); library(EnvStats); …
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
#' @description Small wrapper that hides the branching between the two back‑ends.
#'   It always returns a **named list** whose names are identical to the names
#'   of the input split list. This guarantees that `unlist()` produces the same
#'   named vector for both back‑ends, which is required by the unit tests.
#' @param X List (or vector) of inputs that will be split among the workers.
#' @param FUN Function that will be executed on each chunk.
#' @param backend Either "parallel" or "future".
#' @param use_threads Number of workers (relevant for both back‑ends).
#' @param cl Optional cluster object created with \code{parallel::makeCluster()}.
#'   Only required when using the \code{"parallel"} backend on Windows, where
#'   \code{parallel::parLapply()} is used. Ignored on Unix-like systems where
#'   \code{parallel::mclapply()} is available.
#' @return List of results (identical to the output of `parLapply` /
#'   `future_lapply`). The list is **named** – the names are taken from the
#'   split index (`split_idx`).
#' @keywords internal
#'
run_parallel <- function(X,
                         FUN,
                         backend = c("parallel", "future"),
                         use_threads = NULL,
                         cl = NULL) {
    backend <- match.arg(backend)

    if (backend == "parallel") {
        if (use_threads <= 1L) {
            return(lapply(X, FUN))
        }

        ## ---- Linux / macOS: use fork directly ----
        if (.Platform$OS.type != "windows") {
            return(parallel::mclapply(
                X,
                FUN,
                mc.cores = use_threads,
                mc.preschedule = TRUE
            ))
        }

        ## ---- Windows: use PSOCK cluster ----
        if (is.null(cl)) {
            stop("Cluster object 'cl' must be provided for Windows parallel backend.")
        }

        return(parallel::parLapply(cl, X, FUN))
    }

    ## ---- future backend ----
    return(future.apply::future_lapply(
        X = X,
        FUN = FUN,
        future.seed = FALSE, # deterministic master RNG
        future.globals = FALSE # prevent unnecessary global export
    ))
}
