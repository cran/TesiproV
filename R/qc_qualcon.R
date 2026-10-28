# =============================================================================
#  PROB_QUALCON
# =============================================================================

# Default defect fractions at which the OC curve is evaluated
.oc_default_theta_grid <- c(0.001, 0.005, 0.007, 0.01,
                            0.012, 0.015, 0.02,  0.03,
                            0.05,  0.07,  0.09,  0.12,
                            0.15,  0.2,   0.3,   0.4,
                            0.5,   0.9)

#' @title Sampling plan and OC curve (PROB_QUALCON)
#'
#' @description
#' Defines a quality-control sampling plan and computes the corresponding
#' Operating Characteristic (OC) curve, i.e. the acceptance probability
#' \eqn{P_a(\theta)} as a function of the defect fraction \eqn{\theta}.
#'
#' A \code{PROB_QUALCON} object refers to a \code{PROB_BASEVAR} via the field
#' \code{var}, because the quality control relates to one specific structural
#' variable. The OC curve is used by \code{PROB_BASEVAR$updateBasevar()} to
#' derive the posterior distribution of that variable after acceptance.
#'
#' Two plan types are supported (field \code{ConformityType}):
#' \itemize{
#'   \item \strong{Variable plan} (\code{"fck_EN-206"}): conformity criterion of
#'         EN 206 for the concrete compressive strength. A lot is accepted if
#'         the mean of \code{n} samples is at least \eqn{f_{ck} + 1.48 s} and
#'         the smallest value at least \eqn{f_{ck} - 4}. \code{f_ck} is given
#'         as \code{DefectLimitLower}, \code{ConformityParam = c(n)}. The OC
#'         curve is simulated, the sample size per defect fraction is set
#'         with the field \code{NSim}.
#'   \item \strong{Attribute plan} (\code{"dim_k||n"}): binomial model, a lot is
#'         accepted if at most \code{c_accept} of \code{n} units are defective.
#'         \code{ConformityParam = c(n, c_accept)}. The plan can be used for
#'         any characteristic whose compliance with tolerance values is checked.
#'         A unit is defective if the variable lies below
#'         \code{DefectLimitLower} or above \code{DefectLimitUpper}. Both
#'         limits are freely selectable. They do not influence the OC curve,
#'         only the Bayesian update.
#' }
#'
#' The fields \code{DefectLimitLower} and \code{DefectLimitUpper} define when a
#' single unit counts as defective and thereby the defect fraction
#' \eqn{\theta}. They are not the acceptance criterion. How a lot is accepted
#' is determined by the sampling plan alone.
#'
#' The OC curve can only be computed for a variable \code{var} with
#' \code{DistributionType} \code{"st"} or \code{"lt"} (Student-t and
#' log-Student-t with hyperparameters \eqn{c(m', s', n', v')}).
#'
#' @field Id Numeric identifier.
#' @field Name Short name of the plan (character).
#' @field Description Optional description (character).
#' @field ConformityType Type of the sampling plan, \code{"fck_EN-206"} or
#'   \code{"dim_k||n"} (character).
#' @field ConformityParam Parameters of the sampling plan (numeric). For
#'   \code{"fck_EN-206"}: \code{c(n)} with the sample size \code{n}. For
#'   \code{"dim_k||n"}: \code{c(n, c_accept)}. It is not needed if the OC
#'   curve is specified manually.
#' @field NSim Monte-Carlo sample size per defect fraction for the simulation
#'   of the OC curve (numeric, default 100000). It is only used by
#'   \code{"fck_EN-206"}.
#' @field DefectLimitLower Lower limit (numeric, default \code{-Inf}). A unit
#'   below it is defective. For \code{"fck_EN-206"} it is the characteristic
#'   compressive strength \eqn{f_{ck}} in the unit of the variable and must be
#'   set, also for a manually specified OC curve.
#' @field DefectLimitUpper Upper limit (numeric, default \code{Inf}). A unit
#'   above it is defective. For \code{"fck_EN-206"} it must stay \code{Inf}.
#'   For \code{"dim_k||n"} at least one of the two limits must be finite.
#' @field var The \code{PROB_BASEVAR} the plan applies to (any). It is
#'   \code{NULL} as long as no variable is attached. \code{getOC()} requires it.
#' @field theta_grid Defect fractions at which the OC curve is evaluated
#'   (numeric). A default grid from 0.001 to 0.9 is set if nothing is given.
#'   \code{getOC(force = TRUE)} resets it to this default grid.
#' @field Pa_grid Acceptance probabilities belonging to \code{theta_grid}
#'   (numeric). It is \code{NA} until \code{getOC()} has been run. A manually
#'   specified OC curve can be given as \code{Pa_grid} together with
#'   \code{theta_grid} of the same length.
#' @field .prior_snapshot Internal. Digest of the distribution of \code{var} at
#'   the time of the last \code{getOC()} call, used by \code{isStale()}.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{getOC(force = FALSE)}}{Computes the OC curve and stores it in
#'     \code{Pa_grid}. An existing curve is only recomputed with
#'     \code{force = TRUE}. With \code{force = TRUE} the default
#'     \code{theta_grid} is used again, a grid or curve given by the user is
#'     discarded.}
#'   \item{\code{plot(log_x = TRUE, ...)}}{Plots \code{Pa_grid} against
#'     \code{theta_grid}. Runs \code{getOC()} first if necessary.}
#'   \item{\code{isStale()}}{\code{TRUE} if the distribution of \code{var} has
#'     changed since \code{getOC()}, \code{NA} if no snapshot exists (manual
#'     curve or no variable).}
#'   \item{\code{show()}, \code{summary()}}{Short printed representations.}
#' }
#'
#' @examples
#' \dontrun{
#' v <- PROB_BASEVAR$new(Name = "fc", DistributionType = "lt",
#'                       DistributionParameters = c(3.75, 0.105, 3, 10))
#' plan <- PROB_QUALCON$new(Name = "EN206_n15",
#'                          ConformityType = "fck_EN-206",
#'                          ConformityParam = 15, NSim = 2000,
#'                          DefectLimitLower = 30, var = v)
#' plan$getOC()
#' plan$plot()
#' post <- v$updateBasevar(plan)
#'
#' # Attribute plan: at most 1 of 10 units outside the limits 27 and 33
#' d <- PROB_BASEVAR$new(Name = "cover", DistributionType = "st",
#'                       DistributionParameters = c(30, 2, 3, 10))
#' plan2 <- PROB_QUALCON$new(Name = "cover_plan", ConformityType = "dim_k||n",
#'                           ConformityParam = c(10, 1),
#'                           DefectLimitLower = 27, DefectLimitUpper = 33, var = d)
#' plan2$getOC()
#' post2 <- d$updateBasevar(plan2)
#' }
#' @author (C) 2026 T. Lux, M. Ricker -- TU Dortmund University - Chair of Structural Concrete
#' @export PROB_QUALCON
#' @exportClass PROB_QUALCON
#'
PROB_QUALCON <- setRefClass(
  Class = "PROB_QUALCON",
  fields = list(
    Id              = "numeric",
    Name            = "character",
    Description     = "character",
    ConformityType  = "character",
    ConformityParam = "numeric",
    NSim            = "numeric",
    DefectLimitLower = "numeric",
    DefectLimitUpper = "numeric",
    var             = "ANY",   # PROB_BASEVAR
    theta_grid      = "numeric",
    Pa_grid         = "numeric",
    .prior_snapshot = "character"   # internal: digest of var's state at last getOC()
  ),
  methods = list(
    ## -----------------------------------------------------------------
    ## INITIALIZE - constructor: assign fields - defaults - validation
    ## -----------------------------------------------------------------
    initialize = function(...) {
      initFields(...)

      pa_given    <- !.is_unset(Pa_grid)    && length(Pa_grid)    > 0L
      theta_given <- !.is_unset(theta_grid) && length(theta_grid) > 0L

      if (xor(pa_given, theta_given)) {
        warning(sprintf(
          "PROB_QUALCON '%s': Pa_grid and theta_grid should be supplied together for a manual OC curve - only one of them was given.",
          ifelse(is_empty(Name), "<unnamed>", Name)
        ), call. = FALSE)
      } else if (pa_given && theta_given && length(Pa_grid) != length(theta_grid)) {
        warning(sprintf(
          "PROB_QUALCON '%s': Pa_grid and theta_grid have different lengths (%d vs. %d).",
          ifelse(is_empty(Name), "<unnamed>", Name), length(Pa_grid), length(theta_grid)
        ), call. = FALSE)
      }

      if (!.is_unset(var) && !is(var, "PROB_BASEVAR"))
        stop("'var' must be a PROB_BASEVAR object.", call. = FALSE)

      if (.is_unset(var)) var <<- NULL

      if (is_missing(Id)) Id <<- NA_real_
      if (is_missing(Pa_grid)) Pa_grid <<- NA_real_
      if (is_missing(Name)) Name <<- ""
      if (is_missing(Description)) Description <<- ""
      if (is_missing(DefectLimitLower)) DefectLimitLower <<- -Inf
      if (is_missing(DefectLimitUpper)) DefectLimitUpper <<- Inf

      if (length(DefectLimitLower) != 1L || length(DefectLimitUpper) != 1L ||
          is.na(DefectLimitLower) || is.na(DefectLimitUpper))
        stop("'DefectLimitLower' and 'DefectLimitUpper' must be single numbers (-Inf / Inf for an open side).",
             call. = FALSE)
      if (DefectLimitLower >= DefectLimitUpper)
        stop("'DefectLimitLower' must be smaller than 'DefectLimitUpper'.", call. = FALSE)
      if (is_missing(.prior_snapshot)) .prior_snapshot <<- NA_character_

      if (is_missing(NSim)) NSim <<- 100000
      if (length(NSim) != 1L || !is.finite(NSim) || NSim < 1)
        stop("'NSim' must be a single positive number.", call. = FALSE)
      if (identical(ConformityType, "fck_EN-206") && length(ConformityParam) > 1L)
        stop("'ConformityParam' for 'fck_EN-206' must be c(n). ",
             "The Monte-Carlo sample size is set with the field 'NSim'.",
             call. = FALSE)

      if (is_missing(theta_grid)) theta_grid <<- .oc_default_theta_grid

      invisible(.self)
    },
    ## -----------------------------------------------------------------
    ## Internal helper: fingerprint of var's current distribution state
    ## -----------------------------------------------------------------
    .snapshotVar = function() {
      if (is.null(var)) return(NA_character_)
      digest::digest(
        list(var$DistributionType, var$DistributionParameters),
        algo = "xxhash64"
      )
    },

    # ------------------------------------------------------------------
    #  getOC()  -  compute OC curve
    # ------------------------------------------------------------------
    getOC = function(force = FALSE) {
      "Computes the OC curve and stores it in $Pa_grid. ..."

      if (!all(is.na(Pa_grid)) && !force) {
        message("Pa_grid already available. Use getOC(force = TRUE) to recompute.")
        return(invisible(.self))
      }

      # force = TRUE: start again with the default grid of defect fractions
      if (force) theta_grid <<- .oc_default_theta_grid

      if (is.null(var))
        stop("No PROB_BASEVAR attached - set $var first.", call. = FALSE)
      if (!var$DistributionType %in% c("st", "lt"))
        stop("$var must have DistributionType 'st' or 'lt'.", call. = FALSE)
      if (is.null(ConformityType) || length(ConformityType) == 0L || nchar(ConformityType) == 0L)
        stop("ConformityType not set - define a sampling plan before calling getOC().", call. = FALSE)
      if (is.null(ConformityParam) || length(ConformityParam) == 0L || anyNA(ConformityParam))
        stop("ConformityParam not set - define a sampling plan before calling getOC().", call. = FALSE)


      Pa_grid <<- .oc_compute_pa(
        theta      = theta_grid,
        conf.type  = ConformityType,
        conf.param = ConformityParam,
        dist.type  = var$DistributionType,
        dist.param = var$DistributionParameters,
        lower      = DefectLimitLower,
        upper      = DefectLimitUpper,
        N          = NSim
      )

      .prior_snapshot <<- .snapshotVar()

      message("OC curve computed  |  ConformityType = '", ConformityType,
              "'  |  n = ", ConformityParam[1])
      invisible(.self)
    },

    # ------------------------------------------------------------------
    #  plot()  -  plots OC curve for visual checking
    # ------------------------------------------------------------------
    plot = function(log_x = TRUE, ...) {
      "Plots the OC curve (Pa_grid vs. theta_grid) using base R graphics.
   Runs $getOC() first if Pa_grid has not been computed yet."

      if (all(is.na(Pa_grid))) {
        message("Pa_grid not yet available - running getOC().")
        getOC()
      }

      use_log <- log_x && all(theta_grid > 0)
      if (log_x && !use_log) {
        warning("theta_grid contains non-positive values - using a linear x axis.",
                call. = FALSE)
      }

      title_str <- sprintf("OC curve%s", ifelse(is_empty(Name), "", paste0(" - ", Name)))
      sub_str <- paste0(
        "ConformityType = ", ConformityType,
        if (length(ConformityParam) > 0) paste0(" | n = ", ConformityParam[1]) else "",
        if (!is.null(var)) paste0(" | var: ", var$Name) else ""
      )

      graphics::plot(
        theta_grid, Pa_grid, type = "b", pch = 16,
        log = if (use_log) "x" else "",
        xlab = expression(paste("Defect fraction ", theta)),
        ylab = expression(paste("Acceptance probability ", P[a](theta))),
        main = title_str, sub = sub_str,
        ylim = c(0, 1),
        ...
      )
      invisible(.self)
    },
    # ------------------------------------------------------------------
    #  isStale()  -  check whether $var has changed since getOC()
    # ------------------------------------------------------------------
    isStale = function() {
      "Returns TRUE if $var's distribution has changed since Pa_grid was
       last computed via getOC() (e.g. by later re-assigning Mean/Sd on
       the same object). Returns NA if Pa_grid was set manually or $var
       is missing, since no snapshot exists to compare against."

      if (is.null(var) || is.na(.prior_snapshot)) return(NA)
      !identical(.snapshotVar(), .prior_snapshot)
    },

    # ------------------------------------------------------------------
    #  show()  -  compact print representation
    # ------------------------------------------------------------------
    show = function() {
      cat("Reference class object of class \"PROB_QUALCON\"\n")
      cat("  Id             :", Id, "\n")
      cat("  Name           :", Name, "\n")
      cat("  Description    :", Description, "\n")
      cat("  ConformityType :", ConformityType, "\n")
      cat("  ConformityParam:", paste(ConformityParam, collapse = ", "), "\n")
      if (identical(ConformityType, "fck_EN-206"))
        cat("  NSim           :", NSim, "\n")
      if (ConformityType %in% c("dim_k||n", "fck_EN-206"))
        cat("  Limits         : [", DefectLimitLower, ",", DefectLimitUpper, "]\n")
      cat("  var            :", if (!is.null(var)) var$Name else "<not attached>", "\n")
      cat("  theta_grid     :", theta_grid, "\n")
      cat("  Pa_grid        :", Pa_grid, "\n")
      stale <- isStale()
      if (isTRUE(stale)) cat("  [warning: $var has changed since getOC() was last run]\n")
      invisible(NULL)
    },

    # ------------------------------------------------------------------
    #  summary()
    # ------------------------------------------------------------------
    summary = function() {
      "Prints a short summary of the PROB_QUALCON object."
      cat("PROB_QUALCON\n")
      cat("  Id             :", Id,   "\n")
      cat("  Name           :", Name, "\n")
      cat("  ConformityType :", ConformityType, "\n")
      cat("  ConformityParam:", paste(ConformityParam, collapse = ", "), "\n")
      if (identical(ConformityType, "fck_EN-206"))
        cat("  NSim           :", NSim, "\n")
      if (ConformityType %in% c("dim_k||n", "fck_EN-206"))
        cat("  Limits         : [", DefectLimitLower, ",", DefectLimitUpper, "]\n")
      if (!is.null(var))
        cat("  var            :", var$Name, "(", var$DistributionType, ")\n")
      else
        cat("  var            : not attached\n")
      if (!all(is.na(Pa_grid)))
        cat("  Pa_grid        : available (", length(Pa_grid), "points)\n")
      else
        cat("  Pa_grid        : not yet computed - call $getOC()\n")
      invisible(NULL)
    }
  )
)
