#'
#' @title internal  Helper function to debug more easy
#' @param infoLevel If 0, no Output (just Errors), if 1 little output, if 2 bigger output
#'
#' @param flag Parse additonal info
#' @param values If you check variables then post this into values
#' @param msg here add some extra msg
#' @param type Type can be "INFO" or "ERROR"
#' @param step_interval Optional integer specifying
#'   the interval at which debug messages are printed.
#' @param iter Current iteration index used for
#'   conditional debug output.
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)


debug.print <- function(infoLevel,
                        flag = "",
                        values = character(0),
                        msg = "",
                        type = "INFO",
                        step_interval = 10,
                        iter = NA) {

  if (infoLevel == 0 && type != "ERROR") {
    return(invisible(NULL))
  }

  # nur alle n Schritte loggen bei Level=1
  if (infoLevel == 1 && is.numeric(iter) && !is.na(iter) && (iter %% step_interval != 0)) {
    return(invisible(NULL))
  }

  # Zeitstempel fuer jede Ausgabe
  t <- format(Sys.time(), "%H:%M:%OS3")

  # Ausgabe vorbereiten
  out_msg <- paste0(
    "#", type,
    if (flag != "") paste(" FROM ", flag),
    " - ", t,
    "\t", msg,
    if (length(values) > 0) paste0(" ", paste(values, collapse = " ")) else "",
    "\n"
  )

  # Ausgabeart: INFO oder ERROR
  if (type == "ERROR") {
    warning(out_msg)
  } else {
    message(out_msg)
  }
}

debug.print.getDebugLevel <- function() {
  return(0)
}


info.print <- function(tag, verbose, varnames, values) {
  if (verbose == 2) {
    cat("\r#", tag, ": ", paste(varnames, ":", values, "\t"))
  }
}


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
#' is_empty(NULL) # TRUE
#' is_empty(numeric(0)) # TRUE
#' is_empty(list()) # TRUE
#' is_empty(data.frame()) # TRUE
#' is_empty(1:5) # FALSE
#' is_empty(list(a = 1)) # FALSE
#' is_empty(data.frame(x = 1)) # FALSE
#' @export

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
#' @export

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
