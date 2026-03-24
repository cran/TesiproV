#' @title Package startup - register default parallelisation preferences
#' @description
#' This internal function is executed automatically when the
#' **TesiproV** package is loaded (via `library(TesiproV)`).
#'
#' Setting package-specific default options that can later be used by TesiproV
#' functions to configure parallel
#' execution in a controlled and CRAN-compliant way.
#'
#' Two global options are defined:
#'
#' * **`TesiproV.future.default`** - the default parallel backend.
#'   The default is `"multisession"`, which works on all operating systems
#'   including Windows.
#'
#' * **`TesiproV.max_workers.default`** - the default upper limit for the
#'   number of parallel workers. The default is `4L`.
#'
#' These values serve as internal fallback settings and may be overridden
#' explicitly by user-level function arguments or environment variables
#' within TesiproV functions.
#'
#' The function returns `invisible(NULL)` and is **not exported**. It is
#' intended solely for internal package initialisation.
#'
#' @param libname Ignored. Path to the library where the package is installed
#'   (standard argument for `.onLoad`).
#' @param pkgname Ignored. Name of the package (standard argument for `.onLoad`).
#'
#' @details
#' The function performs two steps:
#'
#' 1. It checks whether the namespace of the package `future` is available
#'    and aborts with a clear error message if it is missing.
#' 2. It registers internal default options via `options()` without modifying
#'    the active `future` plan.
#'
#' By not calling `future::plan()` during package startup, TesiproV avoids
#' altering the user's global parallel configuration. This follows CRAN
#' recommendations that packages should not change global execution plans
#' on load.
#'
#' @examples
#' ## Load the package (initialises default options)
#' library(TesiproV)
#'
#' @keywords internal
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#' @importFrom future availableCores plan sequential multisession multicore
#' @import future.apply
#' @import digest
#' @export

.onLoad <- function(libname, pkgname) {

    if (!requireNamespace("future", quietly = TRUE)) {
      stop("Package \"future\" is required but not installed.",
           call. = FALSE)
    }

    # DO NOT set a plan here
    # Only store default preferences

    options(
      TesiproV.future.default = "multisession",
      TesiproV.max_workers.default = 4L
    )

  invisible()
}
