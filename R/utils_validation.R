#' Test whether an object is *empty*
#'
#' This lightweight helper mimics the behaviour of `pracma::isempty()` but
#' depends only on base R.  An object is considered *empty* when it contains
#' no data:
#'   * `NULL`
#'   * an atomic vector, matrix, array or factor of length 0
#'   * a list (including named lists) of length 0
#'   * a data frame with zero rows (`nrow == 0`)
#'
#' The function is deliberately strict - environments, functions, S4/R6
#' objects, etc. are **not** treated as empty.  If you need additional
#' classes to be regarded as empty you can extend the internal checks.
#'
#' @param x Any R object.
#' @return A logical scalar (`TRUE` if `x` is empty, `FALSE` otherwise).
#' @keywords internal
#' @examples
#' \dontrun{
#' is_empty(NULL) # TRUE
#' is_empty(numeric(0)) # TRUE
#' is_empty(list()) # TRUE
#' is_empty(data.frame()) # TRUE
#' is_empty(1:5) # FALSE
#' is_empty(list(a = 1)) # FALSE
#' is_empty(data.frame(x = 1)) # FALSE
#' }

is_empty <- function(x) {
  ## 1) NULL -------------------------------------------------------------
  if (is.null(x)) {
    return(TRUE)
  }

  ## 2) Character vectors - empty strings count as empty -----------------
  if (is.character(x)) {
    return(length(x) == 0L || all(nchar(x) == 0L))
  }

  ## 3) Atomic, factor, matrix, array ------------------------------------
  if (is.atomic(x) || is.factor(x) || is.matrix(x) || is.array(x)) {
    return(length(x) == 0L)
  }

  ## 4) data.frame -------------------------------------------------------
  if (inherits(x, "data.frame")) {
    return(nrow(x) == 0L)
  }

  ## 5) Lists (including named lists) ------------------------------------
  if (is.list(x)) {
    return(length(x) == 0L)
  }

  ## 6) Anything else (functions, environments, R6/S4 objects,...) ------
  FALSE
}

#' Helper: is a value missing (NULL, NA, zero-length, or empty string)
#' Used inside validateArgs() to avoid the missing value where TRUE/FALSE needed
#' problem.
#' @param x any R object
#' @return TRUE if x should be treated as missing, FALSE otherwise
#'
#' @keywords internal
is_missing <- function(x) {
  if (is.null(x)) {
    return(TRUE)
  }
  if (length(x) == 0L) {
    return(TRUE)
  }
  if (is.atomic(x) && any(is.na(x))) {
    return(TRUE)
  }
  if (is.character(x) && any(nchar(x) == 0L)) {
    return(TRUE)
  }
  FALSE
}

#' @title Check whether a value is NULL or an uninitialized RefClass field
#' @description Reference Class fields of type "ANY" that were not supplied
#'   at construction are NOT NULL, but an object of class "uninitializedField"
#'   (a sentinel from the methods package). A plain is.null() check silently
#'   fails to detect this, so any optional "ANY" field must be tested with
#'   this helper instead.
#' @keywords internal
.is_unset <- function(x) {
  is.null(x) || methods::is(x, "uninitializedField")
}

#' @title Default value for NULL
#' @description Returns \code{y} if \code{x} is \code{NULL}, otherwise
#'   \code{x}. Internal replacement for the operator \code{\%||\%} of
#'   base R (available only from R 4.4.0), so the package also works with
#'   older R versions.
#' @param x Value to test.
#' @param y Default returned if \code{x} is \code{NULL}.
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
