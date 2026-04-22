#' @title System Probability Solution Object (SYS_PROB)
#'
#' @description
#' The \code{SYS_PROB} class represents a probabilistic system
#' consisting of one or more limit-state functions (\code{SYS_LSF})
#' and a set of reliability algorithms (\code{PROB_MACHINE}).
#'
#' It allows the execution of reliability analyses for
#' serial or parallel systems, aggregation of single
#' limit-state results, and system-level probability evaluation.
#'
#' Each limit-state function may include multiple random variables
#' (\code{PROB_BASEVAR}) with different probability distributions.
#' The transformation between mean/standard deviation and the
#' native distribution parameters is handled automatically.
#'
#' \itemize{
#'
#'   \item Normal ("norm"):
#'   mean = \eqn{\mu}, standard deviation = \eqn{\sigma}.
#'   Parameters: mean, sd.
#'
#'   \item Lognormal ("lnorm"):
#'   defined by mean \eqn{m} and sd \eqn{s}.
#'   Parameters:
#'   \eqn{\mu = \log\left( \frac{m}{\sqrt{1 + (s/m)^2}} \right)},
#'   \eqn{\sigma = \sqrt{ \log\left( 1 + (s/m)^2 \right) }}.
#'
#'   \item shifted Lognormal ("slnorm"):
#'   defined by Mean \eqn{m}, Sd \eqn{s} and \eqn{x0}.
#'   Parameters:
#'   \eqn{\mu = \log\left( \frac{m-x0}{\sqrt{1 + (s/(m-x0))^2}} \right)},
#'   \eqn{\sigma = \sqrt{ \log\left( 1 + (s/(m-x0))^2 \right) }}.
#'
#'   \item Gumbel ("gumbel", package "evd"):
#'   location-scale type I extreme value distribution.
#'   Parameters:
#'   location = \eqn{\mu - \gamma \cdot \beta},
#'   where \eqn{\gamma} is the Euler-Mascheroni constant,
#'   and
#'   \eqn{\beta = \frac{\sigma \sqrt{6}}{\pi}}.
#'
#'   \item Gamma ("gamma"):
#'   shape-scale parameterization.
#'   Parameters:
#'   \eqn{k = \frac{\mu^2}{\sigma^2}},
#'   \eqn{\theta = \frac{\sigma^2}{\mu}}.
#'
#'   \item Exponential ("exp"):
#'   rate \eqn{\lambda} or scale \eqn{\theta = 1/\lambda}.
#'   Parameter:
#'   \eqn{\lambda = 1/\mu}.
#'
#'   \item Beta ("beta"):
#'   bounded on \eqn{[0,1]}.
#'   Parameters:
#'   \eqn{\alpha = \mu t},
#'   \eqn{\beta = (1 - \mu) t},
#'   with
#'   \eqn{t = \frac{\mu (1-\mu)}{\sigma^2} - 1}.
#'
#'   \item Weibull ("weibull"):
#'   shape \eqn{k} and scale \eqn{\lambda}.
#'   Parameters:
#'   empirical approximation
#'   \eqn{k \approx \left( \frac{\sigma}{\mu} \right)^{-1.086}},
#'   \eqn{\lambda = \frac{\mu}{\Gamma\left(1 + \frac{1}{k}\right)}}.
#'
#'   \item Empirical ("emp"):
#'   non-parametric distribution defined by observed sample.
#'   Parameter: \code{obs} = numeric vector of observations.
#'   The empirical CDF is constructed from the ordered sample.
#'   Duplicate values are allowed but may trigger a warning.
#'
#'   \item Student t ("st"):
#'   marginal density of a normal-gamma model describing uncertainty
#'   in \eqn{\mu} and \eqn{\sigma}.
#'   Hyperparameters:
#'   \eqn{(m, s, n, \nu)}.
#'
#'   \item Log-Student t ("lt"):
#'   marginal density of a log-normal-gamma model describing uncertainty
#'   in \eqn{\mu} and \eqn{\sigma}.
#'   Hyperparameters:
#'   \eqn{(m, s, n, \nu)}.
#' }
#'
#' The types `"binom"` and `"unif"` are reserved for future extensions.
#' The package automatically checks consistency between mean, standard deviation and coefficient of variation.
#'
#' Methods available in this class:
#'
#' * `$runMachines()` - executes all solution algorithms for each limit state function.
#' * `$calculateSystemProbability()` - computes system reliability using bounds or Monte Carlo based methods (`MC_IS`, `MC_CRUDE`, `MC_SubSam`).
#' * `$printResults()` - prints a detailed report to console or file.
#' * `$saveProject(level)` - saves results at different detail levels.
#' * `$plotGraph(plotType)` - experimental plotting routine for simulation performance.
#'
#'
#' @field sys_input      List of `SYS_LSF` objects - the individual limit-state functions.
#' @field sys_type       Character string that defines whether the system
#'                       is *serial* or *parallel* (not yet implemented).
#' @field probMachines   List of `PROB_MACHINE` objects - the analysis methods.
#' @field res_single     List of results per machine for each limit-state function.
#' @field res_sys        List of results for the system-level reliability calculation.
#' @field beta_single    Matrix of \eqn{\beta}-values obtained from the single-problem analyses.
#' @field beta_sys       Matrix of \eqn{\beta}-values obtained from the system reliability calculations.
#' @field params_sys     List that stores the parameter sets used in Monte-Carlo runs.
#' @field debug.level    Numeric verbosity level (0 = silent, 1 = basic, 2 = detailed).
#'
#' @details
#' The class is built on \code{setRefClass} and therefore uses mutable fields.
#' All objects (`PROB_BASEVAR`, `SYS_LSF`, ...) are expected to be already
#' **checked** (via the `check()` method) before they are passed to `SYS_PROB`.
#'
#' @examples
#' ps <- SYS_PROB(
#'   sys_input = list(SYS_LSF(), SYS_LSF()),
#'   probMachines = list(PROB_MACHINE()),
#'   sys_type = "serial"
#' )
#' \dontrun{
#' ps$runMachines()
#' ps$beta_sys
#' ps$res_sys
#' ps$printResults("example_1")
#' ps$saveProject(4, "example_1")
#' }
#'
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#'
#' @import gridExtra
#' @import ggplot2
#' @importFrom methods "new"
#' @importFrom stats dnorm pnorm qnorm rnorm
#' @importFrom utils str
#' @export SYS_PROB
#' @exportClass SYS_PROB


SYS_PROB <- setRefClass(
  Class = "SYS_PROB",
  fields = list(
    sys_input = "list", # List of SYS_LSFs
    sys_type = "character", # determining serial or parallel system, not implemented yet
    probMachines = "list", # list of PROB_MACHINES
    res_single = "list",
    res_sys = "list",
    beta_single = "matrix",
    beta_sys = "matrix",
    params_sys = "list",
    sys_structure = "ANY",
    smooth_kappa = "numeric",
    debug.level = "numeric"
  )
)


SYS_PROB$methods(
  list(
    initialize = function(...) {
      initFields(...)
      beta_single <<- matrix(, 0, 0)
      beta_sys <<- matrix(, 0, 0)

      # Default smooth parameter
      if (is_empty(smooth_kappa)) {
        smooth_kappa <<- 20
      }

      # Default system structure
      if (is_empty(sys_structure)) {
        sys_structure <<- sys_type # fallback to existing behaviour
      }
      # Only set default if not provided by user
      if (is_empty(debug.level)) {
        debug.level <<- 0
      }
    },
    runMachines = function() {
      "Starts solving all given problems (sys_input) with all given algorithms (probMachines).
       First checks that all analysis methods are loaded before validating limit-state functions."

      # --- Step 1: check existence of reliability machines first ---
      for (j in seq_along(probMachines)) {
        probMachines[[j]]$checkMethodExists()
      }

      # --- Step 2: validate each limit-state function afterwards ---
      for (i in 1:length(sys_input)) {
        sys_input[[i]]$check()
      }

      if (is_empty(debug.level)) {
        debug.level <<- 0
      }

      # --- Safety check: all limit-state functions must be defined ---
      for (i in seq_along(sys_input)) {
        lsf_obj <- sys_input[[i]]

        # Case 1: no function object at all
        if (is_empty(lsf_obj$func)) {
          stop(sprintf(
            "Limit-state function for SYS_LSF[%d] ('%s') is not defined.\n",
            i,
            ifelse(is_empty(lsf_obj$name), "<unnamed>", lsf_obj$name)
          ), call. = FALSE)
        }

        # Case 2: function exists but has no arguments
        f_args <- formals(lsf_obj$func)
        if (is.null(f_args) || length(f_args) == 0L) {
          stop(sprintf(
            "Limit-state function for SYS_LSF[%d] ('%s') has no arguments.\n",
            i,
            ifelse(is_empty(lsf_obj$name), "<unnamed>", lsf_obj$name)
          ), call. = FALSE)
        }
      }

      # Helper: normalize getlDistr() output to a funlist (list with $d/$p/$q/$r)
      normalize_to_funlist <- function(dd) {
        # case 1: already funlist
        if (is.list(dd) && !is.null(dd$p) && is.function(dd$p)) {
          return(dd)
        }
        # case 2: wrapper -> first element is funlist
        if (is.list(dd) && length(dd) >= 1 && is.list(dd[[1]]) && !is.null(dd[[1]]$p) && is.function(dd[[1]]$p)) {
          return(dd[[1]])
        }
        # fallback: find any element that looks like a funlist
        if (is.list(dd)) {
          for (el in dd) {
            if (is.list(el) && !is.null(el$p) && is.function(el$p)) {
              return(el)
            }
          }
        }
        stop(sprintf("SYS_PROB$runMachines: unexpected getlDistr() format for variable (type=%s).", typeof(dd)), call. = FALSE)
      }

      beta_single <<- matrix(nrow = length(probMachines), ncol = length(sys_input))
      colnames(beta_single) <<- vector("character", length(sys_input))
      rownames(beta_single) <<- vector("character", length(probMachines))
      # Loop through each Problem in the System
      for (i in 1:length(sys_input)) {
        if (!is_empty(sys_input[[i]]$expr)) {
          lsf_expr <- sys_input[[i]]$expr
        }

        lsf <- sys_input[[i]]$getLSF()
        # lsf<-sys_input[[i]]$func
        if (is_empty(sys_input[[i]]$name)) {
          sys_input[[i]]$name <<- "Unknown problem name"
        }
        colnames(beta_single)[i] <<- sys_input[[i]]$name

        # Build basic distribution list for probMachines (normalized funlists)
        distr <- lapply(sys_input[[i]]$vars, function(v) {
          normalize_to_funlist(v$getlDistr())
        })

        lDistr_funlist <- distr

        # lDistr_for_FORM: wrapper list(funlist, params) (used only when calling FORM)
        lDistr_for_FORM <- lapply(lDistr_funlist, function(fun) {
          params <- c(NA_real_, NA_real_)
          if (!is.null(fun$mean)) params[1] <- as.numeric(fun$mean)
          if (!is.null(fun$Sd)) params[2] <- as.numeric(fun$Sd)
          list(fun, params)
        })

        # Loop through each probMachine
        res_machine <- list()
        for (j in 1:length(probMachines)) {
          if (probMachines[[j]]$fCall == "MVFOSM") {
            options <- probMachines[[j]]$options
            if (options$isExpression == "TRUE") {
              params <- list("lsf" = lsf_expr, "lDistr" = distr, "debug.level" = debug.level)
            } else {
              params <- list("lsf" = lsf, "lDistr" = distr, "debug.level" = debug.level)
            }
          } else {
            params <- list("lsf" = lsf, "lDistr" = distr, "debug.level" = debug.level)
          }

          if (!is_empty(probMachines[[j]]$options)) {
            params <- c(params, probMachines[[j]]$options)
          }

          # Later, when building params for the machine call:
          # decide which lDistr variant to pass
          # to be clarified if necessary
          if (probMachines[[j]]$fCall %in% c("FORM", "SORM") ||
            grepl("FORM", probMachines[[j]]$fCall, ignore.case = TRUE)) {
            params[["lDistr"]] <- lDistr_for_FORM
          } else {
            params[["lDistr"]] <- lDistr_for_FORM
          }

          probMachines[[j]]$checkMethodExists()

          safe_params <- params[names(params) %in% names(formals(get(probMachines[[j]]$fCall)))]
          res <- do.call(probMachines[[j]]$fCall, safe_params)

          # --- Debug: check structure of returned result -------------------
          if (debug.level >= 2) {
            message(sprintf("[SYS_PROB/runMachines] Structure of result from '%s':", probMachines[[j]]$fCall))
            str(res)
            str(res$beta)
          }
          # -----------------------------------------------------------------
          rownames(beta_single)[j] <<- probMachines[[j]]$name
          beta_single[j, i] <<- res$beta
          res_machine[[j]] <- res
        }
        res_single[[i]] <<- res_machine
      }
    },
    calculateSystemProbability = function(calcType = "simpleBounds", params = list()) {
      "Calculates the system probablity if more than one lsf is given and a system_type (serial or parallel) is set.
      If calcType is empty (or simpleBounds), only simpleBounds are applied to further calculation of single soultions.
      If calcType is MCIS, than a Monte Carlo Importance Sampling Method is used (only for parallel systems available).
      If calcType is MCC, than a Crude Monte Carlo Simulation is used.
      If calcType is MCSUS, than the Subset Sampling Algorithm ll be used.
      You can pass arguments to methods via the params field, while the argument has to be a named list (for example check the vignette)."
      params_sys[[length(params_sys) + 1]] <<- list("calcType" = calcType, params)

      ## ----- 1. ensure a debug level exists ---------------------------------
      if (is_empty(debug.level)) {
        debug.level <<- 0
      }

      ## make sure the debug level that you set globally is also forwarded
      ## (the MC_IS routine expects an argument called `debug.level`)
      # if (!is.null(debug.level)) {
      #   params$debug.level <- debug.level
      # }

      # debug-Message
      if (debug.level >= 2) {
        message("calcType = ", calcType)
        str(params)
      }

      # -------------------- Read and set parameters
      if (exists("cov_user", where = params)) {
        cov_user <- params$cov_user
      } else {
        cov_user <- 0.05
      }
      if (exists("n_batch", where = params)) {
        n_batch <- params$n_batch
      } else {
        n_batch <- 64
      }
      if (exists("n_max", where = params)) {
        n_max <- params$n_max
      } else {
        n_max <- 5e6
      }
      if (exists("use_threads", where = params)) {
        use_threads <- params$use_threads
      } else {
        use_threads <- 8
      }
      if (exists("beta_l", where = params)) {
        beta_l <- params$beta_l
      } else {
        beta_l <- 8
      }

      if (exists("Nsubset", where = params)) {
        Nsubset <- params$Nsubset
      } else {
        Nsubset <- 1e5
      }

      if (exists("variance", where = params)) {
        variance <- params$variance
      } else {
        variance <- "uniform"
      }

      if (exists("p0", where = params)) {
        p0 <- params$p0
      } else {
        p0 <- 0.1
      }

      if (exists("densityType", where = params)) {
        densityType <- params$densityType
      } else {
        densityType <- "norm"
      }

      if (exists("dps", where = params)) {
        dps <- params$dps
      } else {
        dps <- NULL
      }


      if (calcType == "simpleBounds") {
        min_beta <- vector("numeric", nrow(beta_single))
        max_beta <- vector("numeric", nrow(beta_single))

        if (sys_type == "parallel") {
          for (m in 1:nrow(beta_single)) { # m are the machines used
            min_beta[m] <- max(beta_single[m, ])
            max_beta[m] <- -1 * qnorm(pnorm(sum(-beta_single[m, ])))
          }
        } else if (sys_type == "serial") {
          for (m in 1:nrow(beta_single)) { # m are the machines used
            min_beta[m] <- -1 * qnorm(1 - prod(pnorm(beta_single[m, ])))
            max_beta[m] <- min(beta_single[m, ])
          }
        }

        if (nrow(beta_sys) < 1) {
          beta_sys <<- matrix(nrow = 2, ncol = 1)
          beta_sys[1, 1] <<- min_beta
          beta_sys[2, 1] <<- max_beta
          rownames(beta_sys) <<- c("SB min", "SB max")
        } else {
          beta_sys <<- rbind(beta_sys, "SB min" = min_beta)
          beta_sys <<- rbind(beta_sys, "SB max" = max_beta)
        }

        res <- list(
          "method" = "Simple Bounds",
          "system_mode" = sys_type,
          "single_betas" = beta_single,
          "min_beta" = min_beta,
          "max_beta" = max_beta
        )

        res_sys[[length(res_sys) + 1]] <<- res
        info.print("SimpleBounds", debug.level, c("system_mode", "beta_min", "beta_max"), c(sys_type, min_beta, max_beta))
      } else if (calcType == "MCIS") {
        # -------------------------------------------------
        # 1) Collect limit-state functions and distributions
        # -------------------------------------------------

        lsfs <- list()
        distr <- list()

        for (i in seq_along(sys_input)) {
          sys_input[[i]]$check()

          lsfs[[i]] <- sys_input[[i]]$getLSF()

          distr[[i]] <- lapply(sys_input[[i]]$vars, function(v) {
            d <- v$getlDistr()

            if (is.list(d) && length(d) >= 1 &&
              is.list(d[[1]]) && !is.null(d[[1]]$p)) {
              d[[1]]
            } else if (is.list(d) && !is.null(d$p)) {
              d
            } else {
              stop("Unexpected distribution object format returned by getlDistr()")
            }
          })
        }

        is_system <- length(lsfs) > 1

        # -------------------------------------------------
        # 2) Determine system method
        # -------------------------------------------------

        system_method <- if (!is.null(params$system_method)) {
          params$system_method
        } else {
          "multimodal"
        }

        # -------------------------------------------------
        # 3) Flatten distributions (needed for smooth)
        # -------------------------------------------------

        distr_flat <- unlist(distr, recursive = FALSE)

        # -------------------------------------------------
        # 4) Build global smooth LSF if requested
        # -------------------------------------------------

        if (is_system && system_method == "smooth") {
          g_sys <- function(x) {
            g_vals <- numeric(length(lsfs))

            k <- 1

            for (i in seq_along(lsfs)) {
              n_var <- length(sys_input[[i]]$vars)

              g_vals[i] <- lsfs[[i]](x[k:(k + n_var - 1)])

              k <- k + n_var
            }

            if (.self$sys_structure == "parallel") {
              smooth_max(g_vals, .self$smooth_kappa)
            } else {
              smooth_min(g_vals, .self$smooth_kappa)
            }
          }
        }

        # -------------------------------------------------
        # 5) Select MC_IS input
        # -------------------------------------------------

        if (!is_system) {
          lsf_arg <- lsfs[[1]]
          lDistr_arg <- distr[[1]]
        } else if (system_method == "multimodal") {
          lsf_arg <- lsfs
          lDistr_arg <- distr
        } else if (system_method == "smooth") {
          lsf_arg <- g_sys
          lDistr_arg <- distr_flat
        } else {
          stop("system_method must be 'multimodal' or 'smooth'", call. = FALSE)
        }

        # -------------------------------------------------
        # 6) Prepare MC_IS arguments
        # -------------------------------------------------

        seed_val <- params$seed %||% NULL
        backend_val <- params$backend %||% "future"
        dataRecord <- params$dataRecord %||% FALSE

        mc_args <- list(
          lsf = lsf_arg,
          lDistr = lDistr_arg,
          cov_user = cov_user,
          n_batch = n_batch,
          n_max = n_max,
          use_threads = use_threads,
          beta_l = beta_l,
          densityType = densityType,
          dps = dps,
          debug.level = debug.level,
          seed = seed_val,
          backend = backend_val,
          dataRecord = dataRecord,
          sys_type = sys_type
        )

        # -------------------------------------------------
        # 7) Run MC_IS
        # -------------------------------------------------

        res <- do.call(TesiproV::MC_IS, mc_args)

        # -------------------------------------------------
        # 8) Store result
        # -------------------------------------------------

        if (nrow(beta_sys) < 1) {
          beta_sys <<- matrix(nrow = 1, ncol = 1)

          rownames(beta_sys) <<- "MC IS"
        } else {
          beta_sys <<- rbind(beta_sys, "MC IS" = NA)
        }

        beta_sys[nrow(beta_sys), 1] <<- res$beta

        res_sys[[length(res_sys) + 1]] <<- res


        # } else if (calcType == "MCIS") {
        #   ## -------------------------------------------------
        #   ## 1) Build the list of limit-state functions and
        #   ##    marginal distribution objects (unchanged)
        #   ## -------------------------------------------------
        #   lsfs <- list()
        #   distr <- list()

        #   for (i in 1:length(sys_input)) {
        #     sys_input[[i]]$check()
        #     lsfs[[i]] <- sys_input[[i]]$getLSF()
        #     # subdistr <- lapply(sys_input[[i]]$vars, function(v) v$getlDistr())
        #     subdistr <- lapply(sys_input[[i]]$vars, function(v) {
        #       d_wrapped <- v$getlDistr()
        #       # If getlDistr() returned the wrapper list (funlist, params), extract funlist
        #       if (is.list(d_wrapped) && length(d_wrapped) >= 1 && is.list(d_wrapped[[1]]) && !is.null(d_wrapped[[1]]$p)) {
        #         return(d_wrapped[[1]])
        #       }
        #       # otherwise assume d_wrapped already is the funlist
        #       if (is.list(d_wrapped) && !is.null(d_wrapped$p)) {
        #         return(d_wrapped)
        #       }
        #       stop("Unexpected distribution object format returned by getlDistr()")
        #     })

        #     # subdistr <- list()
        #     # for (k in 1:length(sys_input[[i]]$vars)) {
        #     #   v <- sys_input[[i]]$vars[[k]]
        #     #   subdistr[[k]] <- v$getlDistr()
        #     # }
        #     distr[[i]] <- subdistr
        #   }

        #   is_system <- length(lsfs) > 1

        #   system_method <- if (!is.null(params$system_method)) {
        #     params$system_method
        #   } else {
        #     "multimodal"
        #   }

        #   # --- DEBUG: print structure to verify what is passed to MC_IS_system
        #   if (debug.level >= 2) {
        #     message("DEBUG before calling MC_IS_system: lsfs and distr structures:")
        #     str(lsfs)
        #     str(distr)
        #   }


        #   # optional sanity check - will stop early if something is wrong
        #   for (i in seq_along(distr)) {
        #     for (j in seq_along(distr[[i]])) {
        #       if (!is.list(distr[[i]][[j]]) || is.null(distr[[i]][[j]]$p)) {
        #         stop(sprintf("Sanity check failed: distr[[%d]][[%d]] not a funlist", i, j), call. = FALSE)
        #       }
        #     }
        #   }

        #   # -------------------------------------------------
        #   # Build global system LSF g_sys
        #   # -------------------------------------------------
        #   if (is_system) {
        #     g_sys <- function(x) {
        #       # -------------------------------------------------
        #       # Evaluate all individual limit-state functions
        #       # -------------------------------------------------
        #       g_vals <- numeric(length(lsfs))

        #       k <- 1
        #       for (i in seq_along(lsfs)) {
        #         n_var_i <- length(sys_input[[i]]$vars)
        #         g_vals[i] <- lsfs[[i]](x[k:(k + n_var_i - 1)])
        #         k <- k + n_var_i
        #       }

        #       # -------------------------------------------------
        #       # Determine system aggregation structure
        #       # -------------------------------------------------

        #       structure_type <- .self$sys_structure

        #       # -------------------------------------------------
        #       # Case 1: Custom function
        #       # -------------------------------------------------
        #       if (is.function(structure_type)) {
        #         return(structure_type(g_vals))
        #       }

        #       # -------------------------------------------------
        #       # Case 2: Character system structure
        #       # -------------------------------------------------
        #       if (is.character(structure_type) && length(structure_type) == 1L) {
        #         if (structure_type == "serial") {
        #           return(smooth_min(g_vals, .self$smooth_kappa))
        #         }

        #         if (structure_type == "parallel") {
        #           return(smooth_max(g_vals, .self$smooth_kappa))
        #         }

        #         stop(sprintf(
        #           "Unsupported character system structure: '%s'",
        #           structure_type
        #         ))
        #       }

        #       # -------------------------------------------------
        #       # Case 3: Fallback to sys_type
        #       # -------------------------------------------------
        #       sys_type_local <- .self$sys_type

        #       if (is.character(sys_type_local) && length(sys_type_local) == 1L) {
        #         if (sys_type_local == "serial") {
        #           return(smooth_min(g_vals, .self$smooth_kappa))
        #         }

        #         if (sys_type_local == "parallel") {
        #           return(smooth_max(g_vals, .self$smooth_kappa))
        #         }
        #       }

        #       # -------------------------------------------------
        #       # Final fallback (legacy default)
        #       # -------------------------------------------------
        #       return(smooth_min(g_vals, .self$smooth_kappa))
        #     }
        #   }

        #   # -------------------------------------------------
        #   # Flatten distribution list for global LSF
        #   # -------------------------------------------------

        #   distr_flat <- list()

        #   for (i in seq_along(distr)) {
        #     for (j in seq_along(distr[[i]])) {
        #       distr_flat[[length(distr_flat) + 1]] <- distr[[i]][[j]]
        #     }
        #   }

        #   # ---------------------------------------------------------------
        #   # 2) Pull the optional entries from `params` (if they exist)
        #   # ---------------------------------------------------------------
        #   seed_val <- if (!is.null(params$seed)) params$seed else NULL
        #   backend_val <- if (!is.null(params$backend)) params$backend else "future"
        #   dataRecord_val <- if (!is.null(params$dataRecord)) params$dataRecord else FALSE

        #   ## -------------------------------------------------
        #   ## 3) Gather *all* arguments for MC_IS
        #   ## -------------------------------------------------
        #   ## a) The arguments that are always required (hard-coded)
        #   # base_args <- list(
        #   #   lsf = lsfs,
        #   #   lDistr = distr,
        #   #   sys_type = sys_type,
        #   #   cov_user = cov_user,
        #   #   n_batch = n_batch,
        #   #   n_max = n_max,
        #   #   use_threads = use_threads,
        #   #   beta_l = beta_l,
        #   #   densityType = densityType,
        #   #   dps = dps,
        #   #   debug.level = debug.level,
        #   #   seed = seed_val,
        #   #   backend = backend_val, # fallback, wenn nicht angegeben
        #   #   dataRecord = dataRecord_val
        #   # )

        #   # if (is_system) {
        #   #   lsf_arg <- g_sys
        #   #   lDistr_arg <- distr_flat
        #   # } else {
        #   #   lsf_arg <- lsfs[[1]]
        #   #   lDistr_arg <- distr[[1]]
        #   # }

        #   if (is_system) {
        #     if (system_method == "multimodal") {
        #       # true system reliability
        #       lsf_arg <- lsfs
        #       lDistr_arg <- distr
        #     } else if (system_method == "smooth") {
        #       # smooth aggregation
        #       lsf_arg <- g_sys
        #       lDistr_arg <- distr_flat
        #     } else {
        #       stop("Unknown system_method. Use 'multimodal' or 'smooth'.", call. = FALSE)
        #     }
        #   } else {
        #     lsf_arg <- lsfs[[1]]
        #     lDistr_arg <- distr[[1]]
        #   }

        #   base_args <- list(
        #     lsf = lsf_arg,
        #     lDistr = lDistr_arg,
        #     cov_user = cov_user,
        #     n_batch = n_batch,
        #     n_max = n_max,
        #     use_threads = use_threads,
        #     beta_l = beta_l,
        #     densityType = densityType,
        #     dps = dps,
        #     debug.level = .self$debug.level,
        #     seed = seed_val,
        #     backend = backend_val,
        #     dataRecord = dataRecord_val
        #   )

        #   ## b) Append the user-supplied options (the content of `params`)
        #   ##    - but only keep those names that really exist in MC_IS
        #   ##    (otherwise `do.call` would fail with unused argument).
        #   ##    `formals(MC_IS)` returns the vector of formal argument names.
        #   mc_formals <- names(formals(TesiproV::MC_IS))

        #   ## keep only entries whose name appears in the formal list
        #   user_args <- params[
        #     (names(params) %in% mc_formals) & # a) valid for MC_IS
        #       !(names(params) %in% names(base_args)) # b) not duplicate
        #   ]

        #   ## c) Final argument list -> base args + (filtered) user args
        #   mc_args <- c(base_args, user_args)

        #   ## -------------------------------------------------
        #   ## 3) Call the Monte-Carlo routine
        #   ## -------------------------------------------------
        #   if (debug.level >= 2) {
        #     message("Calling MC_IS_system with:")
        #     print(names(mc_args))
        #   }

        #   # Ensure debug.level is always forwarded explicitly
        #   mc_args$debug.level <- debug.level

        #   res <- do.call(TesiproV::MC_IS, mc_args)

        #   ## -------------------------------------------------
        #   ## 4) Store the result (unchanged)
        #   ## -------------------------------------------------
        #   if (nrow(beta_sys) < 1) {
        #     beta_sys <<- matrix(nrow = 1, ncol = 1)
        #     beta_sys[1, 1] <<- res$beta
        #     rownames(beta_sys) <<- "MC IS"
        #   } else {
        #     beta_sys <<- rbind(beta_sys, "MC IS" = res$beta)
        #   }
        #   res_sys[[length(res_sys) + 1]] <<- res
        #   ## MC_CRUDE
      } else if (calcType == "MCC" | calcType == "MCSUS") {
        lsfs <- list()
        distr <- list()
        n <- 1
        n_vars <- vector("numeric", length(sys_input))
        for (i in 1:length(sys_input)) {
          sys_input[[i]]$check()
          lsfs[[i]] <- sys_input[[i]]$getLSF()
          n_vars[i] <- length(sys_input[[i]]$vars)
          for (k in 1:n_vars[i]) {
            v <- sys_input[[i]]$vars[[k]]
            distr[[n]] <- v$getlDistr()
            n <- n + 1
          }
        }
        if (sys_type == "parallel") {
          sys_lsf <- function(x) {
            n_lsfs <- length(lsfs)
            res <- vector("numeric", n_lsfs)
            k <- 1
            for (i in 1:n_lsfs) {
              res[i] <- lsfs[[i]](x[k:(k + n_vars[i] - 1)])
              k <- k + n_vars[i]
            }
            return(max(res))
          }
        } else if (sys_type == "serial") {
          sys_lsf <- function(x) {
            n_lsfs <- length(lsfs)
            res <- vector("numeric", n_lsfs)
            k <- 1
            for (i in 1:n_lsfs) {
              res[i] <- lsfs[[i]](x[k:(k + n_vars[i] - 1)])
              k <- k + n_vars[i]
            }
            return(min(res))
          }
        }
        if (calcType == "MCC") {
          res <- TesiproV::MC_CRUDE(
            lsf = sys_lsf,
            lDistr = distr,
            cov_user = cov_user,
            n_batch = n_batch,
            n_max = n_max,
            use_threads = use_threads,
            debug.level = debug.level
          )
          if (nrow(beta_sys) < 1) {
            beta_sys <<- matrix(nrow = 1, ncol = 1)
            beta_sys[1, 1] <<- res$beta
            rownames(beta_sys) <<- "MC Crude"
          } else {
            beta_sys <<- rbind(beta_sys, "MC Crude" = res$beta)
          }
          res_sys[[length(res_sys) + 1]] <<- res
        } else {
          res <- TesiproV::MC_SubSam(
            lsf = sys_lsf,
            lDistr = distr,
            Nsubset = Nsubset,
            variance = variance,
            p0 = p0,
            debug.level = debug.level
          )
          if (nrow(beta_sys) < 1) {
            beta_sys <<- matrix(nrow = 1, ncol = 1)
            beta_sys[1, 1] <<- res$beta
            rownames(beta_sys) <<- "MC SubSam"
          } else {
            beta_sys <<- rbind(beta_sys, "MC SubSam" = res$beta)
          }
          res_sys[[length(res_sys) + 1]] <<- res
        }
      }
    },
    printResults = function(path = "") {
      "TesiproV can create a report file with all the necessary data for you. If you provide a path (or filename, without ending) it will store
      the data there, otherwise it will report to the console. Set the path via setwd() or check it via getwd()."
      if (!path == "") {
        sink(paste(path, ".txt", sep = ""), type = "output")
        cat("Ergebnisausdruck TesiproV Berechnung\n")
        cat(date())
        cat("\n")
      }

      n_sys <- length(sys_input)
      n_machines <- length(probMachines)
      cat("\n -----------  1. Berechnungsergebnisse (Zusammenfassung)    ---------\n\n")

      cat("1.1 Versagensindicies:\n")
      if (nrow(beta_single) > 0) {
        cat("\nSicherheitsindices bei Einzelversagen:")
        print(beta_single)
      }

      if (nrow(beta_sys) > 0) {
        cat("\nSicherheitsindices bei Systemversagen:")
        print(beta_sys)
      }

      cat("\n1.2.1 Ergebnisse je Problemstellung und Algorithmus\n")

      for (i in 1:n_sys) {
        cat("\n____________________________________________\n")
        cat(sprintf("1.2.1.%d Ergebnisse fuer Problem Nr: %d \t-%s\n", i, i, sys_input[[i]]$name))
        for (j in 1:n_machines) {
          cat(sprintf("Beta: %.4f\tPf: %f\t", res_single[[i]][[j]]$beta, res_single[[i]][[j]]$pf))
          cat(sprintf("%s (%s)\n", probMachines[[j]]$name, probMachines[[j]]$fCall))
          print(res_single[[i]][[j]]$runtime)
        }
      }

      cat("\n1.2.2 Ergebnisse je Systemberechnung und Simulationsmethode\n")
      if (nrow(beta_sys) > 0) {
        for (i in 1:length(res_sys)) {
          cat("\n____________________________________________\n")
          cat(sprintf("1.2.2.%d Ergebnisse Systemsimulation %d \t %s\n", i, i, res_sys[[i]]$method))
          if (res_sys[[i]]$method == "Simple Bounds") {
            cat(sprintf("Min Beta: %.4f\t Max Beta: %f\t", res_sys[[i]]$min_beta, res_sys[[i]]$max_beta))
          } else {
            cat(sprintf("Beta: %.4f\tPf: %f\t", res_sys[[i]]$beta, res_sys[[i]]$pf))
          }
        }
      }

      cat("\n\n-------------------------  2. SYSTEM BESCHREIBUNG ------------------\n")

      cat(sprintf("2.1 Das System umfasste %d Gleichungen die in %s Beziehung stehen.\n\n", n_sys, sys_type))
      for (i in 1:n_sys) {
        cat("\n___________________________________\n")
        cat(sprintf("2.1.%d Gleichung Nr: %d \t-%s\n", i, i, sys_input[[i]]$name))
        if (!is_empty(sys_input[[i]]$expr)) {
          cat("Folgender symbolischer Ausdruck wurde definiert:\n")
          print(sys_input[[i]]$expr)
        }
        if (!is_empty(sys_input[[i]]$func)) {
          cat("Folgende Funktion wurde definiert:\n")
          print(sys_input[[i]]$func)
        }

        n_vars <- length(sys_input[[i]]$vars)
        cat(sprintf("\nBasisvariablen (%d gegeben):\n", n_vars))
        for (k in 1:n_vars) {
          v <- sys_input[[i]]$vars[[k]]
          # cat(sprintf("Name (ID): %s (%d)\tBeschreibung",v$Name,v$Id))

          cat(sprintf(
            "%d.\tName (ID): %s (%d)\tPackage::Verteilungstyp: %s::%s\tMean: %.3f\tSd: %.3f\tCov: %.3f\tx0: %.3f\tVerteilungsparameter: %.5f\t%.5f\n",
            k, v$Name, v$Id, v$Package, v$DistributionType, v$Mean, v$Sd, v$Cov, v$x0, v$DistributionParameters[1], v$DistributionParameters[2]
          ))
        }
      }
      cat("\n\n ")

      cat("\n ------------------------  3. METHODEN BESCHREIBUNG ------------------\n")
      if (length(n_machines) > 0) {
        cat(sprintf("\n\n3.1 Die Berechnung wurde von %d Methoden analysiert.\nIm folgenden werden die Eingangsparameter der Methoden beschrieben.", n_sys))
        for (i in 1:n_machines) {
          cat("\n______________________________________\n")
          cat(sprintf("3.%d Machine Name (Type): %s (%s)\n", i, probMachines[[i]]$name, probMachines[[i]]$fCall))
          if (!is.null(probMachines[[i]]$options)) {
            cat("Die folgenden Optionen waren von den Standartwerten abweichend:\n")
            for (j in 1:length(probMachines[[i]]$options)) {
              cat(sprintf("\t-%s\t=\t%s\n", names(probMachines[[i]]$options[j]), probMachines[[i]]$options[j]))
            }
          } else {
            cat("Ausschliesslich Standartoptionen verwendet!\n")
          }
        }
      }

      if (length(params_sys) > 0) {
        cat(sprintf("\n\n3.2 Parameter der Systemsimulationen."))
        for (i in 1:length(params_sys)) {
          cat("\n______________________________________\n")
          cat(sprintf("3.%d Simulation Method: %s\n", i, params_sys[[i]]$calcType))
          print(params_sys)
        }
      }


      cat("\n -------------------------  AUSFUEHRLICHE ERGEBNISSE ------------------\n")
      cat("\n4.1 Einzelprobleme\n")
      if (length(n_sys) > 0) {
        for (i in 1:n_sys) {
          cat("\n_____________________________________\n")
          cat(sprintf("Ergebnisse fuer Problem Nr: %d \t-%s\n", i, sys_input[[i]]$name))
          for (j in 1:n_machines) {
            cat(sprintf("Machine Name (Type): %s (%s)\n", probMachines[[j]]$name, probMachines[[j]]$fCall))

            n_erg <- length(res_single[[i]][[j]])
            for (k in 1:n_erg) {
              if (names(res_single[[i]][[j]][k]) == "data") {} else if (names(res_single[[i]][[j]][k]) == "runtime") {
                cat("\t -runtime:")
                print(res_single[[i]][[j]][k])
              } else {
                cat(sprintf("\t-%s\t=\t%s\n", names(res_single[[i]][[j]][k]), res_single[[i]][[j]][k]))
              }
            }
          }
        }
      }

      if (length(res_sys) > 0) {
        cat("\n\n4.2 Systemprobleme\n")
        for (i in 1:length(res_sys)) {
          cat("\n_____________________________________\n")
          cat(sprintf("Simulationmethode: %s\n", res_sys[[i]]$method))

          n_erg <- length(res_sys[[i]])
          for (k in 1:n_erg) {
            if (names(res_sys[[i]][k]) == "data") {} else if (names(res_sys[[i]][k]) == "runtime") {
              cat("\t -runtime:")
              print(res_sys[[i]][k])
            } else {
              cat(sprintf("\t-%s\t=\t%s\n", names(res_sys[[i]][k]), res_sys[[i]][k]))
            }
          }
        }
      }


      cat("--------------------- ENDE -----------------")


      if (!path == "") {
        sink()
      }
    },
    saveProject = function(level, filename = "tesiprov_project") {
      "You can save your calculation project with saveProject().
      There are four different levels of detail to save
      1st Level: Only the beta values
      2nd Level: The result Objects of single or systemcalculation
      3th Level: All The Probablity System Object, including limit state functions, machines and solutions
      4th Level: An image of your entire workspace"
      if (level == 1) {
        beta <- list(
          "beta_single" = beta_single,
          "beta_sys" = beta_sys
        )
        filename_loc <- paste(filename, "_beta.rds", sep = "")
        saveRDS(beta, filename_loc)
        cat(paste("Beta values successfully stored in\n", getwd(), "/", filename_loc, ".\nOpen with readRDS().", sep = ""))
      } else if (level == 2) {
        res <- list(
          "res_single" = res_single,
          "res_sys" = res_sys
        )
        filename_loc <- paste(filename, "_res.rds", sep = "")
        saveRDS(res, filename_loc)
        cat(paste("Results successfully stored in\n", getwd(), "/", filename_loc, ".\nOpen with readRDS().", sep = ""))
      } else if (level == 3) {
        ps <- list(
          "sys_input" = sys_input,
          "sys_type" = sys_type,
          "probMachines" = probMachines,
          "res_single" = res_single,
          "res_sys" = res_sys,
          "beta_single" = beta_single,
          "beta_sys" = beta_sys
        )
        filename_loc <- paste(filename, "_ps.rds", sep = "")
        saveRDS(ps, filename_loc)
        cat(paste("Prob. Sys. object successfully saved in\n", getwd(), "/", filename_loc, ".\nOpen with readRDS().", sep = ""))
      } else if (level == 4) {
        filename_loc <- paste(filename, ".RData", sep = "")
        save.image(filename_loc)
        cat(paste("Entire workspace successfully saved in\n", getwd(), "/", filename_loc, ".\nOpen with load().", sep = ""))
      } else {
        cat(paste("Level needs to be between 1 (less output) and 4 (entire workspace).", sep = ""))
      }
    },
    plotGraph = function(plotType = "sim.performance") {
      "not finally implemented. Do not use."
      for (i in 1:length(sys_input)) {
        sysName <- sys_input[[i]]$name
        for (j in 1:length(probMachines)) {
          pmName <- probMachines[[j]]$name
          pmfCall <- probMachines[[j]]$fCall
          if (pmfCall == "MC_IS") {
            if (plotType == "sim.performance") {
              df <- res_single[[i]][[j]]$data

              g1 <- ggplot(df, aes(x = n_sim)) +
                geom_line(aes(y = cov)) +
                labs(
                  x = "Simulationen (-)",
                  y = "CoV (-)",
                  title = paste(pmName, "fuer", sysName)
                ) +
                theme_bw()

              g2 <- ggplot(df, aes(x = n_sim)) +
                geom_line(aes(y = time)) +
                labs(
                  x = "Simulationen (-)",
                  y = "Zeit (s)"
                ) +
                theme_bw()

              ga <- grid.arrange(g1, g2, ncol = 1)
              return(ga)
            } else if (plotType == "sim.beta") {
              df <- res_single[[i]][[j]]$data

              g1 <- ggplot(df, aes(x = n_sim)) +
                geom_line(aes(y = -qnorm(pf)), size = 1) +
                geom_hline(yintercept = res_single[[i]][[j]]$beta, linetype = 2) +
                geom_text(x = 0, y = res_single[[i]][[j]]$beta, vjust = "bottom", hjust = "left", label = paste("Beta =", round(res_single[[i]][[j]]$beta, 4))) +
                labs(
                  x = "Simulationen (-)",
                  y = "Beta (-)",
                  title = paste(pmName, "fuer", sysName)
                ) +
                theme_bw()
              return(g1)
            }
          }
        }
      }

      if (plotType == "sim.hv") {
        hv_func <- res_single[[1]]$hv_func
        plots <- list()
        for (i in 1:length(hv_func)) {
          rval <- 0.0001
          samplerate <- 1e5
          x.min <- hv_func[[i]]$qfun(rval)
          x.max <- hv_func[[i]]$qfun(1 - rval)
          x <- seq(x.min, x.max, length.out = samplerate)
          dhv <- hv_func[[i]]$dfun(x)
          rhv <- hv_func[[i]]$rfun(samplerate)
          df <- data.frame(x = x, d = dhv)
          plots[[i]] <- ggplot(data = df, aes(x = x, y = d)) +
            geom_line() +
            labs(
              x = "Groesse",
              y = "Wahrscheinlichkeit (-)",
              title = paste("Verbunddichte hv_", i, sep = "")
            ) +
            theme_bw()
        }
        return(plots)
      }
    }
  )
)

#' @title System Limit-State Function (SYS_LSF)
#'
#' @description A reference-class that stores a user-defined limit-state function.
#' The function can be written either as an **expression** (`expr`) or
#' as a **regular R function** (`func`).  All random variables that appear
#' in the function must be supplied in the `vars` list as `PROB_BASEVAR`
#' objects.
#'
#' @field name  Optional human-readable name of the limit-state function.
#' @field expr  Symbolic expression (e.g. `expression(f_ck - d_nom)`) that is
#'               later converted to a callable function via `ExpressionToFunction()`.
#' @field func  Actual R function implementing the limit-state equation.
#'              May be `NULL` until defined by user or converted from `expr`.
#' @field vars  List of `PROB_BASEVAR` objects that provide the random variables.
#'
#' @details
#' The method `getLSF()` returns a **closure** that expects a numeric vector
#' `x` (ordered exactly like `vars`) and evaluates the limit-state function.
#'
#' @details
#' For consistency across FORM and Monte Carlo reliability methods,
#' limit-state functions must either accept a single numeric vector
#' argument or explicitly named scalar arguments.
#' The use of \code{function(...)} is not supported.
#'
#' @examples
#' vars <- list(
#'   PROB_BASEVAR(Name = "f_ck", DistributionType = "norm", Mean = 30, Sd = 1),
#'   PROB_BASEVAR(Name = "d_nom", DistributionType = "norm", Mean = 0.2, Sd = 0.01)
#' )
#' lsf <- SYS_LSF(name = "Bending resistance", vars = vars)
#' lsf$ExpressionToFunction()
#'
#' @section Limit-State Function (LSF) API:
#'
#' A limit-state function must follow one of the two supported
#' interface patterns:
#'
#' \describe{
#'   \item{\code{function(x)}}{
#'     A single numeric vector argument \code{x} containing all
#'     basic variables in the order specified in \code{vars}.
#'
#'     Example:
#'     \preformatted{
#'     lsf$func <- function(x) {
#'       x[1] - sum((x[-1]^2) / seq_along(x)[-1])
#'     }
#'     }
#'   }
#'
#'   \item{\code{function(Z, Fy, M, ...)}}{
#'     Explicitly named scalar arguments corresponding to the
#'     variable names defined in \code{vars}.
#'
#'     Example:
#'     \preformatted{
#'     lsf$func <- function(Z, Fy, M) {
#'       Z * Fy - M
#'     }
#'     }
#'   }
#' }
#'
#' The following form is \strong{not supported}:
#'
#' \preformatted{
#' function(...)
#' }
#'
#' Using \code{...} as the sole argument leads to inconsistent
#' behaviour between FORM-based and Monte Carlo-based reliability
#' methods and is therefore prohibited.
#'
#' For maximum robustness and compatibility across all reliability
#' algorithms, the use of the vector form \code{function(x)} is
#' recommended.
#'
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#'
#' @export SYS_LSF
#' @exportClass SYS_LSF
#'
SYS_LSF <- setRefClass(
  Class = "SYS_LSF",
  fields = list(
    name = "character",
    expr = "expression",
    # func = "function",
    func = "ANY",
    vars = "list"
  ),
  methods = list(
    initialize = function(...) {
      initFields(...)
      # Explicitly set func to NULL if not provided or invalid
      if (!is.function(func)) func <<- NULL
      invisible(.self)
    },
    ExpressionToFunction = function() {
      "Transforms a valid expression into a objective function. Need the set of Variables with correct spelled names and IDs"
      s <- deparse(expr)
      s <- substr(s, 12, nchar(s) - 1)
      for (i in 1:length(vars)) {
        symb <- vars[[i]]$Name
        vec <- paste(" x[", vars[[i]]$Id, "] ", sep = "")
        s <- gsub(symb, vec, s)
      }
      func <<- eval(parse(text = paste("f <- function(x) {", s, "}", sep = "")))
    },
    FunctionToExpression = function() {
      # Conversion back is not implemented yet
    },
    check = function() {
      "Checks all variables. You dont need to execute this, since the system object will do anyway."

      # --- Validate that a limit-state function exists and is callable -------

      # Case 1: completely missing or default dummy ---------------------------
      if (is.null(func)) {
        stop(sprintf(
          "Limit-state function '%s' is missing.",
          ifelse(is_empty(name), "<unnamed>", name)
        ), call. = FALSE)
      }

      if (is.function(func) &&
        length(formals(func)) == 0L &&
        identical(body(func), quote(NULL)) &&
        is_empty(expr)) {
        stop(sprintf(
          "Limit-state function '%s' is missing.",
          ifelse(is_empty(name), "<unnamed>", name)
        ), call. = FALSE)
      }

      # Case 2: not callable at all -------------------------------------------
      if (!is.function(func)) {
        stop(sprintf(
          "Limit-state function '%s' is not callable or invalid.",
          ifelse(is_empty(name), "<unnamed>", name)
        ), call. = FALSE)
      }

      # Case 3: formally defined but has no arguments -------------------------
      f_args <- tryCatch(formals(func), error = function(e) NULL)
      if (is.null(f_args) || length(f_args) == 0L) {
        stop(sprintf(
          "Limit-state function '%s' has no formal arguments.",
          ifelse(is_empty(name), "<unnamed>", name)
        ), call. = FALSE)
      }

      # Case 4: function defined only with '...' is not supported ----------
      # Using 'function(...)' leads to inconsistent behavior between
      # FORM-based and Monte Carlo-based reliability methods.
      # Therefore, limit-state functions must either accept:
      #   (a) a single vector argument 'function(x)', or
      #   (b) explicitly named scalar arguments.
      if (length(f_args) == 1L && identical(names(f_args), "...")) {
        stop(sprintf(
          "Limit-state function '%s' must not use '...' as sole argument.\n",
          ifelse(is_empty(name), "<unnamed>", name)
        ), call. = FALSE)
      }

      # --- Check all associated base variables -------------------------------
      for (i in seq_along(vars)) {
        vars[[i]]$prepare()
      }

      invisible(TRUE)
    },
    getLSF = function() {
      # Cache variable names once
      var_names <- vapply(vars, function(v) v$Name, character(1))
      n_vars <- length(var_names)

      lsf.local <- function(x) {
        fml <- formals(func)

        # ------------------------------------------------------------
        # Case 1: function(x)
        # ------------------------------------------------------------
        if (length(fml) == 1) {
          return(func(x))
        }

        # ------------------------------------------------------------
        # Case 2: function(Z, Fy, M, ...)
        # ------------------------------------------------------------

        # Vector case (FORM)
        if (is.vector(x)) {
          args <- as.list(x)
          names(args) <- var_names
          return(do.call(func, args))
        }

        # Matrix case (MC_IS)
        if (is.matrix(x)) {
          # Fast column extraction (no as.list, no do.call)
          cols <- lapply(seq_len(n_vars), function(i) x[, i])
          names(cols) <- var_names

          return(do.call(func, cols))
        }

        stop("Unsupported input type for LSF")
      }

      return(lsf.local)
    }
  )
)


#' @title Object to store the distribution model for basic variables (PROB_BASEVAR)
#'
#' @description Stores the distribution model for a basic random variable that can be used
#' in any limit-state function.  The class knows how to transform between
#' *mean / standard deviation* and the native distribution parameters.
#'
#' @field Id                     Integer identifier (position inside the LSF vector).
#' @field Name                   Human-readable name (e.g. `"f_ck"`).
#' @field Description            Short textual description.
#' @field DistributionType       Distribution identifier (`"norm"`, `"lnorm"`, `"slnorm"`, `"gumbel"`,
#'                               `"gamma"`, `"exp"`, `"beta"`, `"st"`, `"binom"`,`"emp"`,
#'                               `"unif"`, `"weibull"`, `"lt"`).
#' @field Mean                   Mean value of the variable.
#' @field Sd                     Standard deviation.
#' @field Cov                    Coefficient of variation (`Sd / Mean`).
#' @field x0                     Shifting parameter (used for shifted log-normal distribution `"slnorm"`).
#' @field Package                Name of the R package that provides the density
#'                               functions (`"stats"`, `"evd"`, `"EnvStats"` or `"brms"`).
#' @field DistributionParameters Numeric vector with the *native* distribution
#'                               parameters (e.g. `c(mu, sigma)` for a normal).
#' @field .cache                 Private environment that stores cached density
#'                               functions (used internally, **do not** access
#'                               from user code).
#'
#' @details
#' The method `prepare()` calculates missing parameters and checks for consistency.
#' `getlDistr()` returns a list containing *four* functions (`d`, `p`, `q`, `r`)
#' together with meta-information; the result is cached to avoid repeated
#' `loadNamespace()` calls.
#'
#' @examples
#' ## Normal random variable --------------------------------------------------
#' var1 <- PROB_BASEVAR(
#'   Name = "f_ck", DistributionType = "norm",
#'   Mean = 30, Sd = 1
#' )
#' var1$prepare()
#'
#' ## Gumbel random variable (package `evd`) ----------------------------------
#' var2 <- PROB_BASEVAR(
#'   Name = "M", DistributionType = "gumbel",
#'   Package = "evd", Mean = 2000, Sd = 200
#' )
#' var2$prepare()
#'
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#'
#' @export PROB_BASEVAR
#' @exportClass PROB_BASEVAR
#'
PROB_BASEVAR <- setRefClass(
  Class = "PROB_BASEVAR",
  fields = list(
    Id = "numeric",
    Name = "character",
    Description = "character",
    DistributionType = "character",
    Mean = "numeric",
    Sd = "numeric",
    Cov = "numeric",
    x0 = "numeric",
    Package = "character",
    DistributionParameters = "vector",
    .cache = "environment"
  ),
  methods = list(
    ## -----------------------------------------------------------------
    ##  Initialise the *private* cache (called once per object)
    ## -----------------------------------------------------------------
    .initCache = function() {
      .cache <<- new.env(parent = emptyenv())
    },

    ## -----------------------------------------------------------------
    ## INITIALIZE - constructor: assign fields - defaults - validation
    ## -----------------------------------------------------------------
    initialize = function(...) {
      ## Copy user-supplied arguments into object fields
      initFields(...)

      ## --- Early integrity checks ---------------------------------------
      # Case 1: completely empty object (no Mean, no Sd, no parameters)
      if ((is.null(Mean) || length(Mean) == 0L) &&
        (is.null(Sd) || length(Sd) == 0L) &&
        (is.null(DistributionParameters) ||
          all(is.na(DistributionParameters)))) {
        stop(sprintf(
          "PROB_BASEVAR '%s' must define at least Mean/Sd or valid distribution parameters.",
          ifelse(is_empty(Name), "<unnamed>", Name)
        ), call. = FALSE)
      }

      ## ---- Numeric fields: never leave them NULL --------------------
      if (is_missing(Id)) Id <<- NA_real_
      if (is_missing(Mean)) Mean <<- NA_real_
      if (is_missing(Sd)) Sd <<- NA_real_
      if (is_missing(Cov)) Cov <<- NA_real_
      if (is_missing(x0)) x0 <<- NA_real_

      ## 3) character slots - empty strings
      if (is_missing(Name)) Name <<- ""
      if (is_missing(Description)) Description <<- ""

      ## 4) distribution type / package defaults (the only change)
      if (is_missing(DistributionType) || DistributionType == "") {
        DistributionType <<- "norm"
      }
      if (is_missing(Package) || Package == "") {
        Package <<- "stats"
      }
      ## 5) DistributionParameters - always length 2
      if (is.null(DistributionParameters) ||
        length(DistributionParameters) < 2) {
        DistributionParameters <<- c(NA_real_, NA_real_)
      }

      ## ---- Create private cache environment --------------------------------
      .cache <<- new.env(parent = emptyenv())

      ## ---- Validate consistency --------------------------------------------
      validateArgs()

      prepare()

      invisible(.self)
    },


    ## -------------------------------------------------
    ##  validateArgs() - pure sanity checks (no calculations)
    ## -------------------------------------------------
    validateArgs = function() {
      ## ---- distribution type -----------------------------------------------
      allowed <- c(
        "norm", "lnorm", "slnorm", "gumbel", "gamma", "exp",
        "beta", "weibull", "st", "binom", "unif", "lt", "emp"
      )
      if (!(DistributionType %in% allowed)) {
        stop(sprintf("Distribution type '%s' is not known.", DistributionType),
          call. = FALSE
        )
      }

      ## ---- if Sd is missing but Mean & Cov are present ---------------------
      if (is_missing(Sd) && !is_missing(Mean) && !is_missing(Cov)) {
        Sd <<- Mean * Cov
      }

      ## ---- exponential - missing Sd - set equal to Mean --------------------
      if (tolower(trimws(DistributionType)) == "exp" &&
        is_missing(Sd) && !is.na(Mean)) {
        Sd <<- Mean
      }

      ## ---- Sd must be non-negative -----------------------------------------
      if (!is_missing(Sd) && Sd < 0) {
        stop("Standard deviation (Sd) must be non-negative.", call. = FALSE)
      }

      invisible(TRUE)
    },
    getlDistr = function() {
      ## 1)  make sure the object is prepared (sets defaults)
      prepare()

      if (is_missing(x0) || is.na(x0)) x0 <<- 0

      ## 2)  guarantee that Name and Package are not empty
      if (is_empty(Name)) Name <<- "Unknown Basevar Name"
      if (is_empty(Package)) Package <<- "stats"

      ## 3)  load the required package
      ensure_pkg(Package)

      ## 4)  create a private cache if it does not exist yet
      # if (!exists(".cache", inherits = FALSE)) .initCache()

      ## 5)  unique key based on distribution type + parameters
      if (!requireNamespace("digest", quietly = TRUE)) {
        stop("Package 'digest' is required for caching distribution functions.",
          call. = FALSE
        )
      }
      cache_key <- digest::digest(list(
        DistributionType,
        DistributionParameters
      ), algo = "xxhash64")

      ## 6)  return cached functions if they already exist
      if (exists(cache_key, envir = .cache, inherits = FALSE)) {
        cached <- get(cache_key, envir = .cache, inherits = FALSE)
        return(list(cached, DistributionParameters))
      }

      ## 7)  otherwise build the four functions
      if (DistributionType == "emp") {
        d_fun <- function(x) {
          do.call("demp",
            envir = loadNamespace("EnvStats"),
            args = list(x = x, obs = DistributionParameters)
          )
        }
        p_fun <- function(x) {
          do.call("pemp",
            envir = loadNamespace("EnvStats"),
            args = list(q = x, obs = DistributionParameters)
          )
        }
        q_fun <- function(x) {
          do.call("qemp",
            envir = loadNamespace("EnvStats"),
            args = list(p = x, obs = DistributionParameters)
          )
        }
        r_fun <- function(x) {
          do.call("remp",
            envir = loadNamespace("EnvStats"),
            args = list(n = x, obs = DistributionParameters)
          )
        }
      } else if (DistributionType == "slnorm") {
        d_fun <- function(x) {
          do.call("dshifted_lnorm",
            envir = loadNamespace("brms"),
            c(list(x), as.list(DistributionParameters))
          )
        }
        p_fun <- function(x) {
          do.call("pshifted_lnorm",
            envir = loadNamespace("brms"),
            c(list(x), as.list(DistributionParameters))
          )
        }
        q_fun <- function(x) {
          do.call("qshifted_lnorm",
            envir = loadNamespace("brms"),
            c(list(x), as.list(DistributionParameters))
          )
        }
        r_fun <- function(x) {
          do.call("rshifted_lnorm",
            envir = loadNamespace("brms"),
            c(list(x), as.list(DistributionParameters))
          )
        }
      } else if (DistributionType == "st" || DistributionType == "lt") {
        d_fun <- function(x) {
          do.call(paste0("d", DistributionType),
            envir = loadNamespace(Package),
            args = list(x = x, hyper.param = DistributionParameters)
          )
        }
        p_fun <- function(x) {
          do.call(paste0("p", DistributionType),
            envir = loadNamespace(Package),
            args = list(q = x, hyper.param = DistributionParameters)
          )
        }
        q_fun <- function(x) {
          do.call(paste0("q", DistributionType),
            envir = loadNamespace(Package),
            args = list(p = x, hyper.param = DistributionParameters)
          )
        }
        r_fun <- function(x) {
          do.call(paste0("r", DistributionType),
            envir = loadNamespace(Package),
            args = list(n_vals = x, hyper.param = DistributionParameters)
          )
        }
      } else {
        d_fun <- function(x) {
          do.call(paste0("d", DistributionType),
            envir = loadNamespace(Package),
            c(list(x), as.list(DistributionParameters))
          )
        }
        p_fun <- function(x) {
          do.call(paste0("p", DistributionType),
            envir = loadNamespace(Package),
            c(list(x), as.list(DistributionParameters))
          )
        }
        q_fun <- function(x) {
          do.call(paste0("q", DistributionType),
            envir = loadNamespace(Package),
            c(list(x), as.list(DistributionParameters))
          )
        }
        r_fun <- function(x) {
          do.call(paste0("r", DistributionType),
            envir = loadNamespace(Package),
            c(list(x), as.list(DistributionParameters))
          )
        }
      }

      ## 8)  store in cache
      funlist <- list(
        d = d_fun, p = p_fun, q = q_fun, r = r_fun,
        name = Name, mean = Mean, Sd = Sd, X0 = x0,
        Cov = Cov, DistributionType = DistributionType,
        DistributionPackage = Package
      )
      assign(cache_key, funlist, envir = .cache)

      ## 9)  return the structure expected by the rest of the code
      list(funlist, DistributionParameters)
    },

    ## -----------------------------------------------------------------------
    ## Performs transformations between mean/sd and distribution parameters.
    ## Ensures that Mean, Sd and Cov are consistent and valid.
    ## -----------------------------------------------------------------------
    prepare = function() {
      "Performs transformations between mean/sd and distribution parameters.
      Ensures that Mean, Sd and Cov are consistent and valid."

      ## ---- compute missing Sd from Mean & Cov -------------------------------
      if (is.na(Sd) && !is.na(Mean) && !is.na(Cov)) {
        Sd <<- Mean * Cov
      }

      ## ---- exponential shortcut --------------------------------------------
      if (tolower(DistributionType) == "exp" && is.na(Sd) && !is.na(Mean)) {
        Sd <<- Mean
      }
      ## ---------------------------------------------------------------
      ##  Forward: Mean & Sd (Cov darf NA sein) -> native Parameter
      ## ---------------------------------------------------------------
      if (!any(is.na(c(Mean, Sd)))) { # <- Cov is not required here
        switch(tolower(DistributionType),
          norm = {
            DistributionParameters[1] <<- Mean
            DistributionParameters[2] <<- Sd
          },
          lnorm = {
            sn <- sqrt(log(1 + (Sd / Mean)^2))
            mn <- log(Mean / sqrt(1 + (Sd / Mean)^2))
            DistributionParameters[1] <<- mn
            DistributionParameters[2] <<- sn
          },
          slnorm = {
            if (any(is.na(x0))) {
              x0 <<- 0
            } else if (x0 > Mean) {
              stop("For shifted lognormal distribution: Mean must be larger than x0.")
            }
            mx <- Mean
            sn <- sqrt(log(1 + (Sd / (mx - x0))^2))
            mn <- log((mx - x0) / sqrt(1 + (Sd / (mx - x0))^2))
            DistributionParameters[1] <<- mn
            DistributionParameters[2] <<- sn
            DistributionParameters[3] <<- x0
          },
          gumbel = {
            scale <- (Sd * sqrt(6)) / pi
            location <- Mean + digamma(1) * scale
            DistributionParameters[1] <<- location
            DistributionParameters[2] <<- scale
          },
          exp = {
            rate <- ifelse(is.na(Mean), NA_real_, 1 / Mean)
            scale <- Mean
            DistributionParameters[1] <<- rate
            DistributionParameters[2] <<- scale
          },
          beta = {
            var_val <- Sd^2
            t_val <- Mean * (1 - Mean) / var_val - 1
            alpha <- Mean * t_val
            beta_ <- (1 - Mean) * t_val
            DistributionParameters[1] <<- alpha
            DistributionParameters[2] <<- beta_
          },
          gamma = {
            shape <- Mean^2 / Sd^2
            scale <- Sd^2 / Mean
            DistributionParameters[1] <<- shape
            DistributionParameters[2] <<- scale
          },
          weibull = {
            ## 1) Versuche numerische Loesung fuer k (Shape)
            root_fun <- function(k) {
              if (k <= 0) {
                return(Inf)
              }
              cov_k <- sqrt(exp((base::gamma(1 + 2 / k) / base::gamma(1 + 1 / k)^2)) - 1)
              cov_target <- Sd / Mean
              cov_k - cov_target
            }
            ## 2) Intervall waehlen (typisch 0.1 - 100)
            k_est <- tryCatch(
              uniroot(root_fun, interval = c(0.1, 100))$root,
              error = function(e) {
                ## Falls die Root-Suche scheitert - gehe zurueck zur empirischen Formel
                (Sd / Mean)^(-1.086)
              }
            )
            lam_est <- Mean / base::gamma(1 + 1 / k_est)

            DistributionParameters[1] <<- k_est
            DistributionParameters[2] <<- lam_est
          },
          emp = {
            # not required
          },
          st = {
            # not required
          },
          lt = {
            # not required
          },
          {
            warning(sprintf("Forward transformation for '%s' not implemented", DistributionType))
          }
        )
      } else if ((all(is.na(c(Mean, Sd))) && any(!is.na(DistributionParameters)))) {
        if (DistributionType == "lt" && length(DistributionParameters) != 4) {
          stop("Log-Student-t requires 4 hyperparameters: m', s', n', v'")
        }
        switch(tolower(DistributionType),
          norm = {
            Mean <<- DistributionParameters[1]
            Sd <<- DistributionParameters[2]
            Cov <<- ifelse(DistributionParameters[1] == 0, 0, DistributionParameters[2] / DistributionParameters[1])
          },
          lnorm = {
            mn <- DistributionParameters[1]
            sn <- DistributionParameters[2]
            m <- exp(mn + (sn^2) / 2)
            s <- m * sqrt(exp(sn^2) - 1)
            ## use "<<-" for field assignment
            Cov <<- ifelse(m == 0, 0, s / m)
            Sd <<- s
            Mean <<- m
          },
          slnorm = {
            mn <- DistributionParameters[1]
            sn <- DistributionParameters[2]
            x0 <<- DistributionParameters[3]
            m <- x0 + exp(mn + (sn^2) / 2)
            s <- exp(mn + (sn^2) / 2) * sqrt(exp(sn^2) - 1)
            Cov <<- ifelse(m == 0, 0, s / m)
            Sd <<- s
            Mean <<- m
          },
          gumbel = {
            loc <- DistributionParameters[1]
            sc <- DistributionParameters[2]
            m <- loc - digamma(1) * sc
            s <- (pi / sqrt(6)) * sc
            Cov <<- ifelse(m == 0, 0, s / m)
            Sd <<- s
            Mean <<- m
          },
          exp = {
            rate <- DistributionParameters[1]
            m <- ifelse(rate == 0 || is.na(rate), NA_real_, 1 / rate)
            s <- m
            covv <- ifelse(m == 0 || !is.finite(s) || !is.finite(m),
              NA_real_, s / m
            )
            Cov <<- covv
            Sd <<- s
            Mean <<- m
          },
          gamma = {
            shape <- DistributionParameters[1]
            scale <- DistributionParameters[2]
            m <- shape * scale
            s <- sqrt(shape) * scale
            covv <- ifelse(m == 0 || !is.finite(s) || !is.finite(m),
              NA_real_, s / m
            )
            Cov <<- covv
            Sd <<- s
            Mean <<- m
          },
          beta = {
            ## Extract alpha and beta parameters from DistributionParameters -------------------
            alpha <- DistributionParameters[1]
            beta_ <- DistributionParameters[2]

            ## Validate parameters -----------------------------------------
            if (is.na(alpha) || is.na(beta_) || alpha <= 0 || beta_ <= 0) {
              Mean <<- NA_real_
              Sd <<- NA_real_
              Cov <<- NA_real_
            } else {
              ## Compute mean and standard deviation from parameters -------
              m <- alpha / (alpha + beta_)
              var_val <- (alpha * beta_) /
                ((alpha + beta_)^2 * (alpha + beta_ + 1))

              s <- sqrt(var_val)
              covv <- ifelse(m == 0, 0, s / m)

              ## Assign results to object fields ---------------------------
              Mean <<- m
              Sd <<- s
              Cov <<- covv
            }
          },
          weibull = {
            k <- DistributionParameters[1]
            lam <- DistributionParameters[2]
            Mean <<- lam * base::gamma(1 + 1 / k)
            Sd <<- sqrt(lam^2 *
              (base::gamma(1 + 2 / k) -
                base::gamma(1 + 1 / k)^2))
          },
          emp = {
            # https://en.wikipedia.org/wiki/Empirical_distribution_function
            # obs (sample) is given via DistributionParameters
            # https://rdrr.io/cran/EnvStats/man/Empirical.html

            Mean <<- mean(DistributionParameters)
            Sd <<- sd(DistributionParameters)
            Cov <<- Sd / Mean
          },
          st = {
            m <- DistributionParameters[1]
            s <- DistributionParameters[2]
            n <- DistributionParameters[3]
            nue <- DistributionParameters[4]

            Mean <<- if (nue > 1) m else NA
            Sd <<- if (nue > 2) s / sqrt(n / (n + 1)) * sqrt(nue / (nue - 2)) else if (nue > 1 && nue <= 2) Inf else NA
            # Sd <<- s/sqrt(n/(n+1))/sqrt(nue/(nue-2))
            Cov <<- if (nue > 2) Sd / Mean else NA
          },
          lt = {
            Mean <<- Inf # Median will be used as Start Value for FORM
            Sd <<- Inf
            Cov <<- NA_real_
          },
          {
            warning(sprintf("Backward transformation for '%s' not implemented", DistributionType))
          }
        )
      }

      # --- Recalculate CoV safely -------------------------------------------
      eps <- .Machine$double.eps^0.5 # ~1e-8

      if (is.finite(Mean) && is.finite(Sd) && abs(Mean) > eps) {
        Cov <<- Sd / Mean
      } else if (is.finite(Sd) && abs(Mean) <= eps) {
        Cov <<- Inf
      } else {
        Cov <<- NA_real_
      }

      invisible(TRUE)
    }
  )
)


#' @title Object to store a deterministic parameters (PROB_DETVAR)
#'
#' @description A convenience subclass of `PROB_BASEVAR` that represents a deterministic
#' value (i.e. a variable with practically zero variance).  The object is
#' automatically converted to a normal distribution with an infinitesimal
#' standard deviation (`Sd = Mean / 1e7`).
#'
#' @section Inherited fields:
#' \describe{
#' Inherits all fields of `PROB_BASEVAR` (see that documentation for a full list).
#' }
#'
#' @section Additional fields:
#' \describe{
#'   \item{\code{Value}}{
#'     The deterministic value used internally as \code{Mean}.
#'   }
#' }
#'
#' @examples
#' det <- PROB_DETVAR(Name = "E_mod", Value = 30e3)
#' det$prepare()
#'
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#'
#' @export PROB_DETVAR
#' @exportClass PROB_DETVAR
PROB_DETVAR <- setRefClass(
  Class = "PROB_DETVAR",
  contains = "PROB_BASEVAR",
  fields = c(
    Value = "numeric"
  ),
  methods = list(
    initialize = function(...) {
      initFields(...)
      DistributionType <<- "norm" # Eigentlich nicht notwendig, da PROB_BASE var bei DistributionTyp = "" auf "norm" setzt
      if (length(Value) != 1) {
        stop("PROB_DETVAR: 'Value' must be a scalar numeric.")
      }

      Value <<- as.numeric(Value)[1]

      Mean <<- Value
      Sd <<- Mean / 1e7

      Mean <<- as.numeric(Mean)[1]
      Sd <<- as.numeric(Sd)[1]
      Cov <<- if (Mean == 0) 0 else Sd / Mean

      prepare()
    }
  )
)


#' @title Reliability Algorithm (PROB_MACHINE)
#'
#' @description Stores the definition of a reliability method (FORM, SORM, Monte-Carlo,
#' etc.) that can be executed by `SYS_PROB` or `SYS_PARAM`.  The actual
#' algorithm is called via the character string `fCall` (e.g. `"FORM"` or
#' `"MC_IS"`).  Optional arguments are passed through the list `options`.
#'
#' @field name    Descriptive name of the machine (appears in reports).
#' @field fCall   Function name that implements the algorithm (must be
#'                available in the package namespace, e.g. `FORM`,
#'                `SORM`, `MC_CRUDE`, `MC_IS`, `MC_SubSam`).
#' @field options Additional options for the algorithm (list, e.g.
#'                `list(n_max = 1e6, cov_user = 0.05)`).
#'
#' @details
#' The helper method `getMethodLevel()` returns an integer that indicates
#' the *complexity* of the method (1 = first-order, 2 = second-order, 3 = Monte-Carlo).
#'
#' @examples
#' ## FORM machine ------------------------------------------------------------
#' form_rf <- PROB_MACHINE(
#'   name = "FORM - Rack-Fies",
#'   fCall = "FORM",
#'   options = list(n_optim = 20, loctol = 1e-3, optim_type = "rackfies")
#' )
#'
#' ## Monte-Carlo importance sampling -----------------------------------------
#' mcis <- PROB_MACHINE(
#'   name = "MC-IS",
#'   fCall = "MC_IS",
#'   options = list(cov_user = 0.05, n_max = 3e5, seed = 1234)
#' )
#'
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#'
#' @export PROB_MACHINE
#' @exportClass PROB_MACHINE

PROB_MACHINE <- setRefClass(
  Class = "PROB_MACHINE",
  fields = list(
    name    = "character",
    fCall   = "character",
    options = "list" # <-- changed from "vector"
  ),
  methods = list(
    initialize = function(name = "", fCall = "", options = list()) {
      .self$name <<- name
      .self$fCall <<- fCall
      .self$options <<- options

      ## ---- FORM / SORM need n_optim and loctol -------------------------
      ## If the user did not provide them we inject default values.
      if (grepl("^FORM$", fCall, ignore.case = TRUE) ||
        grepl("^SORM$", fCall, ignore.case = TRUE)) {
        if (is.null(.self$options$n_optim)) .self$options$n_optim <- 20L
        if (is.null(.self$options$loctol)) .self$options$loctol <- 1e-4
      }
    },
    getMethodLevel = function() {
      if (grepl("MC", fCall)) {
        return(3L)
      } else if (grepl("ORM", fCall)) {
        return(2L)
      } else if (identical(fCall, "MVFOSM")) {
        return(2L)
      } else {
        return(1L)
      }
    },

    ## -------------------------------------------------
    ##  Helper used by SYS_PROB$runMachines() and SYS_PARAM$runMachines()
    ##  to produce a **uniform** error message.
    ## -------------------------------------------------
    checkMethodExists = function() {
      ## `exists()` looks for a *function* in the global namespace.
      ## If it is not there we abort with the exact wording the tests expect.
      if (!exists(fCall, mode = "function")) {
        stop(sprintf("Method '%s' is not loaded", fCall), call. = FALSE)
      }
    }
  )
)

# utils.R (oder am Anfang einer deiner .R-Dateien)
#' Ensure that a required package is available
#'
#' @param pkg Character string with the package name
#' @return Invisible TRUE (used only for its side-effect)
#' @note The function stops with an informative error if the required
#'   package is not installed.
#' @keywords internal
ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf(
        "Package \"%s\" is required but not installed.\n",
        pkg
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}


#' Smooth approximation of the minimum function
#'
#' Computes a differentiable approximation of \code{min(x)}
#' using a log-sum-exp formulation:
#'
#' \deqn{
#'   \min(x) \approx -\frac{1}{\kappa} \log\left( \sum_i e^{-\kappa x_i} \right)
#' }
#'
#' The approximation becomes exact as \code{kappa → ∞}.
#'
#' @param x Numeric vector.
#' @param kappa Positive smoothing parameter controlling
#'   the sharpness of the approximation (default = 20).
#'
#' @return Numeric scalar.
#'
#' @details
#' This function is used internally for smooth system
#' limit-state aggregation (serial systems) to improve
#' numerical stability and differentiability.
#'
#' @keywords internal
#' @noRd
smooth_min <- function(x, kappa = 20) {
  -log(sum(exp(-kappa * x))) / kappa
}

#' Smooth approximation of the maximum function
#'
#' Computes a differentiable approximation of \code{max(x)}
#' using a log-sum-exp formulation:
#'
#' \deqn{
#'   \max(x) \approx \frac{1}{\kappa} \log\left( \sum_i e^{\kappa x_i} \right)
#' }
#'
#' The approximation becomes exact as \code{kappa → ∞}.
#'
#' @param x Numeric vector.
#' @param kappa Positive smoothing parameter controlling
#'   the sharpness of the approximation (default = 20).
#'
#' @return Numeric scalar.
#'
#' @details
#' This function is used internally for smooth system
#' limit-state aggregation (parallel systems) to improve
#' numerical stability and differentiability.
#'
#' @keywords internal
#' @noRd
smooth_max <- function(x, kappa = 20) {
  log(sum(exp(kappa * x))) / kappa
}
