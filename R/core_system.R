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
#'   \eqn{k} is determined numerically from the coefficient of variation,
#'   with the empirical approximation
#'   \eqn{k \approx \left( \frac{\sigma}{\mu} \right)^{-1.086}}
#'   as fallback,
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
#'
#'
#' @field sys_input      List of `SYS_LSF` objects - the individual limit-state functions.
#' @field sys_type       Character string, `"serial"` or `"parallel"`. Used by
#'                       `calculateSystemProbability()`.
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
#' @importFrom methods "new"
#' @importFrom stats dnorm pnorm qnorm rnorm
#' @importFrom utils str
#' @export SYS_PROB
#' @exportClass SYS_PROB


SYS_PROB <- setRefClass(
  Class = "SYS_PROB",
  fields = list(
    sys_input = "list", # List of SYS_LSFs
    sys_type = "character", # "serial" or "parallel"
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
      if (.is_unset(sys_structure) || is_empty(sys_structure)) {
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
            if (isTRUE(as.logical(options$isExpression))) {
              params <- list("lsf" = lsf_expr, "lDistr" = distr, "debug.level" = debug.level)
            } else {
              params <- list("lsf" = lsf, "lDistr" = distr, "debug.level" = debug.level)
            }
          } else {
            params <- list("lsf" = lsf, "lDistr" = distr, "debug.level" = debug.level)
          }

          if (!is_empty(probMachines[[j]]$options)) {
            # options of the machine take precedence over the system settings (e.g. debug.level)
            opts <- probMachines[[j]]$options
            params[names(opts)] <- opts
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
      "Calculates the system probability if more than one lsf is given and a system_type (serial or parallel) is set.
      If calcType is empty (or simpleBounds), only simple bounds are computed from the single results of runMachines() (the machines can be FORM, SORM, ...).
      If calcType is MCIS, an importance sampling method for systems is used (serial and parallel systems).
      If calcType is MCC, a crude Monte Carlo simulation of the system limit-state function is used.
      If calcType is MCSUS, subset simulation of the system limit-state function is used.
      A system analysis with FORM or SORM is not available.
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

      # debug message
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
      use_threads_was_supplied <- exists("use_threads", where = params)

      if (use_threads_was_supplied) {
        use_threads <- params$use_threads
      } else {
        use_threads <- NULL
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
          for (m in seq_len(nrow(beta_single))) { # m are the machines used
            min_beta[m] <- max(beta_single[m, ])
            # lower pf bound for independent components: product of the pf_i
            max_beta[m] <- -stats::qnorm(
              sum(stats::pnorm(-beta_single[m, ], log.p = TRUE)),
              log.p = TRUE
            )
          }
        } else if (sys_type == "serial") {
          for (m in seq_len(nrow(beta_single))) { # m are the machines used
            min_beta[m] <- -1 * qnorm(1 - prod(pnorm(beta_single[m, ])))
            max_beta[m] <- min(beta_single[m, ])
          }
        }

        sb_names <- rownames(beta_single)
        if (is.null(sb_names) || any(sb_names == "")) {
          sb_names <- paste0("machine", seq_len(nrow(beta_single)))
        }
        sb <- rbind(
          matrix(min_beta, ncol = 1, dimnames = list(paste0("SB min (", sb_names, ")"), NULL)),
          matrix(max_beta, ncol = 1, dimnames = list(paste0("SB max (", sb_names, ")"), NULL))
        )
        beta_sys <<- if (nrow(beta_sys) < 1) sb else rbind(beta_sys, sb)

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

        # Variables that appear in several limit-state functions are
        # identified by their name and enter the simulation only once.
        var_names <- character(0)
        distr_flat <- list()
        idx_map <- vector("list", length(sys_input))
        for (i in seq_along(sys_input)) {
          vars_i <- sys_input[[i]]$vars
          idx_map[[i]] <- integer(length(vars_i))
          for (k in seq_along(vars_i)) {
            nm <- vars_i[[k]]$Name
            pos <- if (is_empty(nm) || !nzchar(nm)) NA_integer_ else match(nm, var_names)
            if (is.na(pos)) {
              var_names <- c(var_names, if (is_empty(nm)) "" else nm)
              distr_flat[[length(var_names)]] <- distr[[i]][[k]]
              pos <- length(var_names)
            }
            idx_map[[i]][k] <- pos
          }
        }

        # -------------------------------------------------
        # 4) Build global smooth LSF if requested
        # -------------------------------------------------

        if (is_system && system_method == "smooth") {
          g_sys <- function(x) {
            g_vals <- numeric(length(lsfs))

            for (i in seq_along(lsfs)) {
              g_vals[i] <- lsfs[[i]](x[idx_map[[i]]])
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
          beta_l = beta_l,
          densityType = densityType,
          dps = dps,
          debug.level = debug.level,
          seed = seed_val,
          backend = backend_val,
          dataRecord = dataRecord,
          sys_type = sys_type
        )

        if (identical(backend_val, "parallel")) {
          mc_args$use_threads <- use_threads
        } else if (use_threads_was_supplied) {
          # Forward only a user-supplied value so MC_IS can issue its warning.
          mc_args$use_threads <- use_threads
        }

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
        #   #   backend = backend_val, # fallback if not specified
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
        # Variables that appear in several limit-state functions are
        # identified by their name and enter the simulation only once.
        lsfs <- list()
        distr <- list()
        var_names <- character(0)
        idx_map <- vector("list", length(sys_input))
        for (i in seq_along(sys_input)) {
          sys_input[[i]]$check()
          lsfs[[i]] <- sys_input[[i]]$getLSF()
          idx_map[[i]] <- integer(length(sys_input[[i]]$vars))
          for (k in seq_along(sys_input[[i]]$vars)) {
            v <- sys_input[[i]]$vars[[k]]
            pos <- if (is_empty(v$Name) || !nzchar(v$Name)) NA_integer_ else match(v$Name, var_names)
            if (is.na(pos)) {
              var_names <- c(var_names, if (is_empty(v$Name)) "" else v$Name)
              distr[[length(var_names)]] <- v$getlDistr()
              pos <- length(var_names)
            }
            idx_map[[i]][k] <- pos
          }
        }
        sys_lsf <- function(x) {
          g <- vapply(seq_along(lsfs), function(i) lsfs[[i]](x[idx_map[[i]]]), numeric(1))
          if (sys_type == "parallel") max(g) else min(g)
        }
        if (calcType == "MCC") {
          backend_mcc <- params$backend %||% "future"

          mcc_args <- list(
            lsf = sys_lsf,
            lDistr = distr,
            cov_user = cov_user,
            n_batch = n_batch,
            n_max = n_max,
            backend = backend_mcc,
            debug.level = debug.level
          )

          # Pass use_threads to the parallel backend, or when explicitly supplied.
          if (identical(backend_mcc, "parallel") || use_threads_was_supplied) {
            mcc_args$use_threads <- use_threads
          }

          res <- do.call(TesiproV::MC_CRUDE, mcc_args)
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
        on.exit(sink(), add = TRUE)
        cat("TesiproV calculation report\n")
        cat(date())
        cat("\n")
      }

      n_sys <- length(sys_input)
      n_machines <- length(probMachines)
      cat("\n -----------  1. SUMMARY OF RESULTS  ---------\n\n")

      cat("1.1 Reliability indices:\n")
      if (nrow(beta_single) > 0) {
        cat("\nReliability indices for single failure modes:")
        print(beta_single)
      }

      if (nrow(beta_sys) > 0) {
        cat("\nReliability indices for system failure:")
        print(beta_sys)
      }

      cat("\n1.2.1 Results per problem and algorithm\n")

      for (i in seq_len(n_sys)) {
        cat("\n____________________________________________\n")
        cat(sprintf("1.2.1.%d Results for problem no.: %d \t-%s\n", i, i, sys_input[[i]]$name))
        for (j in seq_len(n_machines)) {
          cat(sprintf("Beta: %.4f\tPf: %f\t", res_single[[i]][[j]]$beta, res_single[[i]][[j]]$pf))
          cat(sprintf("%s (%s)\n", probMachines[[j]]$name, probMachines[[j]]$fCall))
          print(res_single[[i]][[j]]$runtime)
        }
      }

      cat("\n1.2.2 Results per system calculation and simulation method\n")
      if (nrow(beta_sys) > 0) {
        for (i in 1:length(res_sys)) {
          cat("\n____________________________________________\n")
          cat(sprintf("1.2.2.%d Results of system simulation %d \t %s\n", i, i, res_sys[[i]]$method))
          if (res_sys[[i]]$method == "Simple Bounds") {
            cat(sprintf("Min Beta: %.4f\t Max Beta: %f\t", res_sys[[i]]$min_beta, res_sys[[i]]$max_beta))
          } else {
            cat(sprintf("Beta: %.4f\tPf: %f\t", res_sys[[i]]$beta, res_sys[[i]]$pf))
          }
        }
      }

      cat("\n\n-------------------------  2. SYSTEM DESCRIPTION ------------------\n")

      cat(sprintf("2.1 The system comprised %d equations in a %s relationship.\n\n", n_sys, sys_type))
      for (i in seq_len(n_sys)) {
        cat("\n___________________________________\n")
        cat(sprintf("2.1.%d Equation no.: %d \t-%s\n", i, i, sys_input[[i]]$name))
        if (!is_empty(sys_input[[i]]$expr)) {
          cat("The following symbolic expression was defined:\n")
          print(sys_input[[i]]$expr)
        }
        if (!is_empty(sys_input[[i]]$func)) {
          cat("The following function was defined:\n")
          print(sys_input[[i]]$func)
        }

        n_vars <- length(sys_input[[i]]$vars)
        cat(sprintf("\nBasic variables (%d given):\n", n_vars))
        for (k in 1:n_vars) {
          v <- sys_input[[i]]$vars[[k]]
          # cat(sprintf("Name (ID): %s (%d)\tBeschreibung",v$Name,v$Id))

          cat(sprintf(
            "%d.\tName (ID): %s (%d)\tPackage::DistributionType: %s::%s\tMean: %.3f\tSd: %.3f\tCov: %.3f\tx0: %.3f\tDistribution parameters: %.5f\t%.5f\n",
            k, v$Name, v$Id, v$Package, v$DistributionType, v$Mean, v$Sd, v$Cov, v$x0, v$DistributionParameters[1], v$DistributionParameters[2]
          ))
        }
      }
      cat("\n\n ")

      cat("\n ------------------------  3. METHOD DESCRIPTION ------------------\n")
      if (n_machines > 0) {
        cat(sprintf("\n\n3.1 The calculation was analysed by %d methods.\nThe input parameters of the methods are described below.", n_machines))
        for (i in seq_len(n_machines)) {
          cat("\n______________________________________\n")
          cat(sprintf("3.%d Machine Name (Type): %s (%s)\n", i, probMachines[[i]]$name, probMachines[[i]]$fCall))
          if (length(probMachines[[i]]$options) > 0) {
            cat("The following options differ from the default values:\n")
            for (j in seq_along(probMachines[[i]]$options)) {
              cat(sprintf("\t-%s\t=\t%s\n", names(probMachines[[i]]$options[j]), probMachines[[i]]$options[j]))
            }
          } else {
            cat("Only default options were used.\n")
          }
        }
      }

      if (length(params_sys) > 0) {
        cat("\n\n3.2 Parameters of the system simulations.")
        for (i in 1:length(params_sys)) {
          cat("\n______________________________________\n")
          cat(sprintf("3.%d Simulation Method: %s\n", i, params_sys[[i]]$calcType))
          print(params_sys)
        }
      }


      cat("\n -------------------------  DETAILED RESULTS ------------------\n")
      cat("\n4.1 Single problems\n")
      if (length(n_sys) > 0) {
        for (i in seq_len(n_sys)) {
          cat("\n_____________________________________\n")
          cat(sprintf("Results for problem no.: %d \t-%s\n", i, sys_input[[i]]$name))
          for (j in seq_len(n_machines)) {
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
        cat("\n\n4.2 System problems\n")
        for (i in 1:length(res_sys)) {
          cat("\n_____________________________________\n")
          cat(sprintf("Simulation method: %s\n", res_sys[[i]]$method))

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


      cat("--------------------- END -----------------")
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
#' Sign convention: \eqn{g(x) > 0} is the safe domain, \eqn{g(x) = 0} the
#' limit state and \eqn{g(x) < 0} failure. All reliability methods of
#' TesiproV use this convention.
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
#' lsf <- SYS_LSF(name = "Bending resistance", expr = expression(f_ck - d_nom), vars = vars)
#' lsf$ExpressionToFunction()
#' lsf$func(c(30, 0.2))
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
      if (length(expr) == 0) {
        stop("ExpressionToFunction(): no expression defined, set 'expr' first.", call. = FALSE)
      }
      e <- if (is.expression(expr)) expr[[1]] else expr
      # replace the variable symbols by x[i] (i = position in vars, same order as the input vector)
      mapping <- lapply(seq_along(vars), function(i) bquote(x[.(i)]))
      names(mapping) <- vapply(vars, function(v) v$Name, character(1))
      f <- function(x) NULL
      body(f) <- do.call(substitute, list(e, mapping))
      environment(f) <- globalenv()
      func <<- f
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
#' The limit-state function (`lsf`) and the basevariables (`lDistr`) are supplied by
#' the system (`SYS_PROB` or `SYS_PARAM`) when the machine is run. They must not be
#' part of `options`. `options` only holds the remaining arguments of the algorithm.
#' Arguments that are not set use the defaults of the algorithm.
#'
#' Arguments that can be set via `options` (defaults in brackets, see the help
#' page of the respective function for details):
#'
#' \tabular{ll}{
#'   \strong{fCall} \tab \strong{Options} \cr
#'   `FORM` \tab `n_optim` (100), `loctol` (1e-5), `optim_type` ("auto"), `lsf_list`, `debug.level` \cr
#'   `SORM` \tab `debug.level` \cr
#'   `MVFOSM` \tab `h` (1e-4), `isExpression` (FALSE), `debug.level` \cr
#'   `MC_CRUDE` \tab `cov_user` (0.05), `n_batch` (10000), `n_max` (1e7), `use_threads`,
#'     `backend`, `dataRecord` (TRUE), `debug.level`, `seed` \cr
#'   `MC_IS` \tab as `MC_CRUDE`, plus `sys_type`, `beta_l`, `densityType`, `dps`,
#'     `adaptive_alpha`, `alpha_update_rate`, `adaptive_batch`, `batch_control`,
#'     `min_adapt_samples`, `alpha_min`, `stability_mode` \cr
#'   `MC_SubSam` \tab `Nsubset` (1e5), `p0` (0.1), `MaxSubsets` (10), `Alpha` (0.05),
#'     `variance` ("uniform"), `debug.level`, `seed` \cr
#' }
#'
#' The helper method `getMethodLevel()` returns the level of the method:
#' 2 for MVFOSM, FORM and SORM, 3 for the Monte Carlo methods.
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
