#' Run code with a temporary Future plan
#'
#' Temporarily sets the active Future plan for a test and restores the
#' previously active plan when the code finishes, including when it errors.
#'
#' @param strategy Future strategy to use while evaluating \code{code}.
#' @param code Code to evaluate using the temporary Future plan.
#' @param ... Additional arguments passed to \code{future::plan()}, such
#'   as \code{workers}.
#'
#' @return The result of evaluating \code{code}.
#'
#' @keywords internal
#' @noRd
with_test_future_plan <- function(strategy, code, ...) {
    old_plan <- future::plan()

    on.exit(
        future::plan(old_plan),
        add = TRUE
    )

    future::plan(strategy, ...)
    force(code)
}
