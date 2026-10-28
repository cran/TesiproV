#' @title Parametric Random Variable (PARAM_BASEVAR)
#'
#' @description Extends `PROB_BASEVAR` so that one of the distribution parameters (Mean,
#' Sd or DistributionType) can be **swept** over a user-defined vector of values.
#'
#' @field ParamValues  Numeric vector holding the values that will be assigned
#'                     during the sweep.
#' @field ParamType    Which parameter is varied: `"Mean"`, `"Sd"` or
#'                     `"DistributionType"`.
#' @field pos          Current index inside `ParamValues` (used internally).
#'
#' @details
#' The method `nextParam()` advances `pos` (circularly) and updates the basic
#' fields (`Mean`, `Sd`, `DistributionType`, `DistributionParameters`) before
#' calling `prepare()` again.  `getCurrentParam()` returns the *last* value that
#' was used (convenient for readable run names).
#'
#' @examples
#' pvar <- PARAM_BASEVAR(
#'   Name = "E_mod",
#'   DistributionType = "norm",
#'   Sd = 3000,
#'   ParamType = "Mean",
#'   ParamValues = c(30e3, 35e3, 40e3)
#' )
#' pvar$getCurrentParam() # first value (30000), set by the constructor
#' pvar$nextParam()       # sets Mean = 35000 and prepares the distribution
#' pvar$getCurrentParam()
#'
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#'
#' @export PARAM_BASEVAR
#' @exportClass PARAM_BASEVAR

PARAM_BASEVAR <- setRefClass(
  Class = "PARAM_BASEVAR",
  contains = "PROB_BASEVAR",
  fields = list(
    ParamValues = "vector",
    ParamType = "character",
    pos = "numeric"
  )
)

PARAM_BASEVAR$methods(
  list(
    initialize = function(...) {
      ## ---------------------------------------------------------------
      ## 1) Take over fields and basic initialisation
      ## ---------------------------------------------------------------
      initFields(...)

      ## ---- Ensure numeric fields ---------------------------
      if (is_missing(Id)) Id <<- NA_real_
      if (is_missing(Mean)) Mean <<- NA_real_
      if (is_missing(Sd)) Sd <<- NA_real_
      if (is_missing(Cov)) Cov <<- NA_real_
      if (is_missing(x0)) x0 <<- NA_real_
      if (is_missing(Name)) Name <<- ""
      if (is_missing(Description)) Description <<- ""

      ## ---- DistributionType / Package with defaults ------------------
      if (is_missing(DistributionType) || is.na(DistributionType) || DistributionType == "") {
        DistributionType <<- "norm"
      }
      Package <<- .resolve_package(Package, DistributionType)

      ## ---- Parameter vector always has length >= 2 --------------------------
      if (is.null(DistributionParameters) || length(DistributionParameters) < 2) {
        DistributionParameters <<- c(NA_real_, NA_real_)
      }

      .cache <<- new.env(parent = emptyenv())

      ## ---------------------------------------------------------------
      ## 2) Check integrity of parametric variable before validation
      ## ---------------------------------------------------------------

      # Ensure position counter exists
      if (is.null(pos) || length(pos) == 0L) {
        pos <<- 1L
      }

      # --- Early consistency check ------------------------------------
      # Case 1: ParamValues completely missing or empty
      if (is.null(ParamValues) || length(ParamValues) == 0L) {
        stop(sprintf(
          "PARAM_BASEVAR '%s' must define a non-empty ParamValues vector.",
          ifelse(is_empty(Name), "<unnamed>", Name)
        ), call. = FALSE)
      }

      # Case 2: Neither Mean nor Sd provided and no ParamValues -> abort
      if (all(is.na(c(Mean, Sd))) && length(ParamValues) == 0L) {
        stop(sprintf(
          "PARAM_BASEVAR '%s' must define either Mean/Sd or valid ParamValues.",
          ifelse(is_empty(Name), "<unnamed>", Name)
        ), call. = FALSE)
      }

      # Case 3: ParamType not specified -> meaningless sweep
      if (is_missing(ParamType) || is.na(ParamType) || ParamType == "") {
        stop(sprintf(
          "PARAM_BASEVAR '%s' requires a valid 'ParamType' (e.g., 'Mean', 'Sd', 'DistributionType').",
          ifelse(is_empty(Name), "<unnamed>", Name)
        ), call. = FALSE)
      }

      ## ---------------------------------------------------------------
      ## Proceed with validation and initial parameter setup
      ## ---------------------------------------------------------------
      validateArgs() # safe now - all fields exist and make sense

      if (length(ParamValues) == 0L || all(is.na(ParamValues))) {
        stop(sprintf("PARAM_BASEVAR '%s' requires non-empty ParamValues.", Name),
             call. = FALSE
        )
      }

      # Set first parameter value immediately to initialise Mean/Sd etc.
      nextParam()

      invisible(.self)
    },
    nextParam = function() {
      ## --- Make sure ParamValues exist -----------------
      if (length(ParamValues) == 0L || all(is.na(ParamValues))) {
        stop(sprintf("No parameter values defined for PARAM_BASEVAR '%s'.", Name),
             call. = FALSE
        )
      }

      ## --- Prepare numeric fields -------------------------------
      if (is_missing(Mean)) Mean <<- NA_real_
      if (is_missing(Sd)) Sd <<- NA_real_
      if (is_missing(Cov)) Cov <<- NA_real_
      if (is_missing(x0)) x0 <<- NA_real_

      ## ensure DistributionParameters always numeric vector of length 2
      if (is.null(DistributionParameters) || length(DistributionParameters) < 2) {
        DistributionParameters <<- c(NA_real_, NA_real_)
      }

      ## --- Wrap-around logic -----------------------------------------
      pos <<- (pos %% length(ParamValues)) + 1L
      i <- if (pos == 1L) length(ParamValues) else pos - 1L

      DistributionParameters <<- c(NA_real_, NA_real_) # Reset

      ## --- Set new parameter depending on type ------------------------
      if (tolower(trimws(ParamType)) == "mean") {
        Mean <<- ParamValues[i]
        if (!is_missing(Sd) && !is_missing(Mean) && Mean != 0) {
          Cov <<- as.numeric(Sd / Mean)
        } else if (!is_missing(Mean) && Mean == 0) {
          Cov <<- 0
        } else {
          Cov <<- NA_real_
        }
      } else if (tolower(trimws(ParamType)) == "sd") {
        Sd <<- ParamValues[i]
        if (!is_missing(Sd) && !is_missing(Mean) && Mean != 0) {
          Cov <<- as.numeric(Sd / Mean)
        } else if (!is_missing(Mean) && Mean == 0) {
          Cov <<- 0
        } else {
          Cov <<- NA_real_
        }
      } else if (tolower(trimws(ParamType)) == "distributiontype") {
        DistributionType <<- as.character(ParamValues[i])
        Package <<- .resolve_package(Package, DistributionType)
      } else {
        warning(sprintf(
          "Unknown ParamType '%s' in PARAM_BASEVAR '%s'.",
          ParamType, Name
        ))
      }

      # Safety check before prepare()
      stopifnot(is.numeric(Cov), is.numeric(Sd), is.numeric(Mean))

      prepare()

      # ## Clear cache - safeguard against stale functions ----------
      # if (!is.null(.cache)) {
      #   rm(list = ls(envir = .cache), envir = .cache)
      # }

      invisible(.self)
    },
    getCurrentParam = function() {
      idx <- if (pos == 1L) length(ParamValues) else pos - 1L
      ParamValues[idx]
    }
  )
)


#' @title Parametric Deterministic Variable (PARAM_DETVAR)
#'
#' @description Like `PROB_DETVAR` but the deterministic value itself can be swept.
#' The class inherits from `PROB_BASEVAR` and therefore re-uses the same
#' transformation machinery.
#'
#' @field ParamValues  Vector of deterministic values that will be assigned
#'                     sequentially in a parametric sweep.
#' @field pos          Current index (used internally).
#'
#' @examples
#' pdet <- PARAM_DETVAR(
#'   Name = "E_mod",
#'   ParamValues = c(30e3, 35e3, 40e3)
#' )
#' pdet$getCurrentParam() # first value (30000), set by the constructor
#' pdet$nextParam()       # sets Mean = 35000 and prepares the distribution
#' pdet$getCurrentParam()
#'
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#'
#' @export PARAM_DETVAR
#' @exportClass PARAM_DETVAR

PARAM_DETVAR <- setRefClass(
  Class = "PARAM_DETVAR",
  contains = "PROB_BASEVAR",
  fields = list(
    ParamValues = "vector",
    pos = "numeric"
  )
)

PARAM_DETVAR$methods(
  list(
    initialize = function(...) {
      callSuper(...)
      DistributionType <<- "norm" # Not strictly necessary, since PROB_BASEVAR sets "norm" for DistributionType = ""
      pos <<- 1L
      if (!is_empty(ParamValues)) nextParam()
    },
    nextParam = function() {
      ## 1) Determine next index - wrap-around
      pos <<- (pos %% length(ParamValues)) + 1L # <- compact code

      ## 2) The index actually used is now pos-1 (the previous value)
      i <- if (pos == 1L) length(ParamValues) else pos - 1L

      ## 3) Reset and insert the parameters
      DistributionParameters <<- c(0, 0)


      Mean <<- ParamValues[i]
      Sd <<- max(abs(Mean) / 1e7, 1e-12)
      Cov <<- ifelse(Mean == 0, 0, Sd / abs(Mean))

      ## 4) Finish - internal preparation
      prepare()
      if (!is.null(.cache)) rm(list = ls(envir = .cache), envir = .cache)
    },
    getCurrentParam = function() {
      idx <- if (pos == 1L) length(ParamValues) else pos - 1L
      ParamValues[idx]
    }
  )
)
