#' @title internal  Helper function to debug more easy
#' @param infoLevel If 0, no Output (just Errors), if 1 little output, if 2 bigger output
#'
#' @param flag Parse additional info
#' @param values If you check variables then post this into values
#' @param msg here add some extra msg
#' @param type Type can be "INFO" or "ERROR"
#' @param step_interval Optional integer specifying
#'   the interval at which debug messages are printed.
#' @param iter Current iteration index used for
#'   conditional debug output.
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#' @keywords internal
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

  # at level 1, log only every n-th step
  if (infoLevel == 1 && is.numeric(iter) && !is.na(iter) && (iter %% step_interval != 0)) {
    return(invisible(NULL))
  }

  # timestamp for every output
  t <- format(Sys.time(), "%H:%M:%OS3")

  # prepare output
  out_msg <- paste0(
    "#", type,
    if (flag != "") paste(" FROM ", flag),
    " - ", t,
    "\t", msg,
    if (length(values) > 0) paste0(" ", paste(values, collapse = " ")) else "",
    "\n"
  )

  # output type: INFO or ERROR
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
