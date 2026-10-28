#' @title Object for parametric Studies (SYS_PARAM)
#'
#' @description Extends `SYS_PROB` to perform a **parametric sweep** over one or more
#' `PARAM_BASEVAR` objects.  For each combination of parameter values a full
#' reliability analysis is executed and the results are stored in
#' `beta_params` and `res_params`.
#'
#' @field beta_params   List of \eqn{\beta}-matrices (one matrix per parameter run).
#' @field res_params    List of detailed result objects (one per run).
#' @field nParams       Numeric vector: number of parameter values for each
#'                       `PARAM_BASEVAR`.
#' @field nLsfs         Number of limit-state functions in the system.
#' @field nMachines     Number of analysis machines (methods) used.
#'
#' @details
#' The method `runMachines()` is almost identical to the one in `SYS_PROB`,
#' but it loops over **all** parameter values, creates a readable run name,
#' and stores the whole result hierarchy.
#'
#' @examples
#' # Simple 2-parameter sweep ------------------------------------------------
#' # two random variables that will be varied
#' pv1 <- PARAM_BASEVAR(
#'   Name = "E_mod",
#'   DistributionType = "norm",
#'   Sd = 3000,
#'   ParamType = "Mean",
#'   ParamValues = c(30e3, 35e3, 40e3)
#' )
#' pv2 <- PARAM_BASEVAR(
#'   Name = "f_ck",
#'   DistributionType = "norm",
#'   Sd = 5,
#'   ParamType = "Mean",
#'   ParamValues = c(30, 35, 40)
#' )
#'
#' lsf <- SYS_LSF(vars = list(pv1, pv2), name = "Bending")
#' ps <- SYS_PARAM(
#'   sys_input = list(lsf),
#'   probMachines = list(PROB_MACHINE(name = "FORM", fCall = "FORM")),
#'   sys_type = "serial"
#' )
#'
#' \dontrun{
#' ps$runMachines()
#' ps$beta_params # matrix of \eqn{\beta}-values for every run
#' }
#'
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#'
#' @import future.apply
#' @import future
#' @export SYS_PARAM
#' @exportClass SYS_PARAM

SYS_PARAM <- setRefClass(
  Class = "SYS_PARAM",
  contains = "SYS_PROB",
  fields = list(
    beta_params = "list", # list of beta-matrices (one per parameter run)
    res_params = "list", # list of full result objects (one per run)
    nParams = "vector", # vector: number of parameter values per PARAM_BASEVAR
    nLsfs = "numeric", # total number of limit-state functions in the system
    nMachines = "numeric" # total number of analysis machines (methods)
  ),
  methods = list(
    ## -------------------------------------------------
    ##  INITIALIZER - create empty containers
    ## -------------------------------------------------
    initialize = function(...) {
      # --- DEBUG: capture all user-supplied arguments ---
      args <- list(...)

      # --- Step 1: manually assign inherited fields from parent class ---
      # We do not use initFields(), because it fails during reference-class creation.
      if (!is.null(args$sys_input)) sys_input <<- args$sys_input
      if (!is.null(args$probMachines)) probMachines <<- args$probMachines
      if (!is.null(args$sys_type)) sys_type <<- args$sys_type

      # --- Step 2: call parent initializer to set up default values ---
      callSuper(...)

      # --- Step 3: initialize own fields for parametric study bookkeeping ---
      beta_params <<- list()
      res_params <<- list()
      nParams <<- numeric()
      nLsfs <<- length(sys_input)
      nMachines <<- length(probMachines)
    },

    ## -------------------------------------------------
    ##  runMachines - perform the parametric study
    ## -------------------------------------------------
    runMachines = function() {
      ## ---------------------------------------------------------------
      ## Guard checks - prevent "subscript out of bounds"
      ## ---------------------------------------------------------------
      if (length(sys_input) == 0L) {
        stop("No limit-state functions ('sys_input') defined", call. = FALSE)
      }
      if (length(probMachines) == 0L) {
        stop("No reliability machines ('probMachines') defined", call. = FALSE)
      }

      ## ---------------------------------------------------------------
      ## Additional integrity check: ensure PARAM_BASEVAR objects still have ParamValues
      lapply(sys_input, function(lsf) {
        lapply(lsf$vars, function(v) {
          if (inherits(v, "PARAM_BASEVAR") && length(v$ParamValues) == 0L) {
            stop(sprintf(
              "PARAM_BASEVAR '%s' must define a non-empty ParamValues vector.",
              v$Name
            ), call. = FALSE)
          }
        })
      })


      # #Check integrity of basic variables first
      # for(i in 1:length(sys_input)){
      #   sys_input[[i]]$check()
      # }

      ## 1. Verify the integrity of all basic variables (faster)
      lapply(sys_input, function(ls) ls$check())

      res_full_loc <- list()
      p <- 0
      paramVar_i <- 0
      # runNames <- vector()

      nLsfs <<- length(sys_input)
      nMachines <<- length(probMachines)

      ## -----------------------------------------------------------------
      ## 1)  Compute total number of runs (needed for pre-allocation)
      ## -----------------------------------------------------------------
      totalRuns <- 0L
      for (m in seq_along(sys_input)) {
        for (v in sys_input[[m]]$vars) {
          if (inherits(v, "PARAM_BASEVAR")) {
            if (length(v$ParamValues) == 0L) {
              stop(sprintf(
                "Parametric variable \"%s\" (ID=%s) has no `ParamValues`.",
                v$Name, as.character(v$Id)
              ), call. = FALSE)
            }
            totalRuns <- totalRuns + length(v$ParamValues)
          }
        }
      }

      runNames <- character(totalRuns) # pre-allocated vector
      p <- 0L

      ## Result containers (one beta matrix and one list object per run)
      beta_params <<- list()
      res_params <<- list()

      ## Helper vector storing the number of values of each parameter variable
      nParams <<- numeric()

      ## 3. Loop over every LSF and then over every variable belonging to it
      for (m in seq_along(sys_input)) {
        for (n in 1:length(sys_input[[m]]$vars)) {
          var <- sys_input[[m]]$vars[[n]]

          ## Only handle PARAM_BASEVAR objects
          if (inherits(var, "PARAM_BASEVAR")) {
            paramVar_i <- paramVar_i + 1L
            nParams[paramVar_i] <<- length(var$ParamValues)

            ## 4. Iterate over all values of the current parameter variable
            for (u in seq_along(var$ParamValues)) {
              p <- p + 1L
              runNames[p] <- sprintf(
                "Paramrun No. %d. Paramname %s with value(Type): %f (%s)",
                p, var$Name, var$getCurrentParam(), var$ParamType
              )
              ## -------------------------------------------------
              ## Containers for results of this single run
              ## -------------------------------------------------
              beta_loc <- matrix(
                NA_real_,
                nrow = length(probMachines),
                ncol = length(sys_input),
                dimnames = list(
                  NULL, # row names (set later)
                  sapply(sys_input, function(ls) {
                    if (is_empty(ls$name)) ls$name <- "Unknown Problem Name"
                    ls$name
                  })
                )
              )

              colnames(beta_loc) <- vector("character", length(sys_input))
              rownames(beta_loc) <- vector("character", length(probMachines))
              # Loop through each Problem in the System

              ## Store the results of each LSF for this run
              res_full_loc <- list()

              for (i in 1:length(sys_input)) {
                ## Optional: use the symbolic expression instead of the compiled function
                lsf_expr <- NULL
                if (!is_empty(sys_input[[i]]$expr)) {
                  lsf_expr <- sys_input[[i]]$expr
                }

                ## Compiled objective function (always available)
                lsf <- sys_input[[i]]$getLSF()
                # lsf<-sys_input[[i]]$func

                ## Ensure every LSF has a name (required for matrix column names)
                if (is_empty(sys_input[[i]]$name)) {
                  sys_input[[i]]$name <<- "Unknown Problem Name"
                }
                colnames(beta_loc)[i] <- sys_input[[i]]$name

                ## Build a list of distribution objects for all variables used in this LSF (Faster)
                distr <- lapply(sys_input[[i]]$vars, function(v) v$getlDistr())

                ## -------------------------------------------------
                ## Parallel execution of all machines for this LSF
                ## -------------------------------------------------
                ## library(future.apply)   # make sure it is loaded

                #  A) Before future_lapply:
                res_machine <- vector("list", length(probMachines))

                res_machine_raw <- future.apply::future_lapply(
                  seq_along(probMachines),
                  function(j) {
                    pm <- probMachines[[j]]

                    ## ---- Build the argument list (identical to the original code) ----
                    if (identical(pm$fCall, "MVFOSM") &&
                        !is_empty(pm$options) &&
                        isTRUE(as.logical(pm$options$isExpression))) {
                      args <- list(lsf = lsf_expr, lDistr = distr, debug.level = debug.level)
                    } else {
                      args <- list(lsf = lsf, lDistr = distr, debug.level = debug.level)
                    }

                    ## ---- Append user-defined options (if any) ------------------------
                    ## options of the machine take precedence over the system settings
                    if (!is_empty(pm$options)) args[names(pm$options)] <- pm$options

                    ## ---- Existence check - **must be done before** `do.call()` -----
                    pm$checkMethodExists()

                    ## Pass only arguments the algorithm knows (as in SYS_PROB)
                    args <- args[names(args) %in% names(formals(get(pm$fCall)))]

                    ## Run the algorithm
                    res <- do.call(pm$fCall, args)

                    ## ---- Optional debug output (only when requested) ---------------
                    if (debug.level >= 2) {
                      message(sprintf(
                        "[SYS_PARAM/runMachines] %s -> beta = %g",
                        pm$fCall, res$beta
                      ))
                    }

                    ## Return a small list that contains everything we need later
                    list(
                      name    = pm$name,
                      beta    = res$beta,
                      result  = res
                    )
                  },
                  future.seed = TRUE
                )

                ## -------------------------------------------------
                ## Copy results back into the beta matrix and the list
                ## -------------------------------------------------

                if (!exists("res_machine")) res_machine <- vector("list", length(probMachines))

                for (j in seq_along(res_machine_raw)) {
                  mech <- res_machine_raw[[j]]
                  beta_loc[j, i] <- mech$beta
                  res_machine[[j]] <- mech$result
                  rownames(beta_loc)[j] <- mech$name
                }

                ## Save the results of this LSF
                res_full_loc[[i]] <- res_machine
                # end of LSF loop
              }

              ## -------------------------------------------------
              ## Save the complete results of this parameter run
              ## -------------------------------------------------
              beta_params[[p]] <<- beta_loc
              res_params[[p]] <<- res_full_loc

              ## set the parametric basic variable one step further
              var$nextParam()

              ## check if there are deterministic vars that also need to be shifted
              for (z in 1:length(sys_input[[m]]$vars)) {
                if (class(sys_input[[m]]$vars[[z]])[1] == "PARAM_DETVAR") {
                  sys_input[[m]]$vars[[z]]$nextParam()
                }
              }
            } # end of u-loop (all values of the current parameter variable)
          }
        } # end of n-loop (variables inside a single LSF)
      } # end of m-loop (over all LSFs)

      ## -------------------------------------------------
      ## Attach the run names to the result lists for easier post-processing
      ## -------------------------------------------------
      names(beta_params) <<- runNames
      names(res_params) <<- runNames
    },

    ## -------------------------------------------------
    ##  printResults - adapted to the parametric study output
    ## -------------------------------------------------
    printResults = function(path = "") {
      if (!path == "") {
        sink(paste(path, ".txt", sep = ""), type = "output")
        on.exit(sink(), add = TRUE)
        cat("TesiproV calculation report\n")
        cat(date())
        cat("\n")
      }

      n_sys <- length(sys_input)
      n_machines <- length(probMachines)
      cat("\n -----------  1. Summary of Results  ---------\n\n")

      cat("\n1.1 Results per solution algorithm\n")
      for (i in seq_len(n_machines)) {
        cat("\n____________________________________________\n")
        for (j in seq_along(res_params)) {
          loc_res <- res_params[[j]][[1]][[i]]
          cat(sprintf("Run %d\tBeta: %.4f\tPf: %f\tMethod: %s\tRuntime: %s\n", j, loc_res$beta, loc_res$pf, loc_res$method, loc_res$runtime))
        }
      }

      cat("\n\n-------------------------  2. SYSTEM DESCRIPTION ------------------\n")

      cat("2.1 The system comprises one equation.\n\n")
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

          if (class(v)[1] == "PARAM_BASEVAR") {
            cat("\n____________________")
            cat("\nParameter variable:\n")
            cat(sprintf(
              "%d.\tName (ID): %s (%d)\tPackage::DistributionType: %s::%s\tMean: %.2f\tSd: %.2f\tCov: %.2f\tx0: %.2f\tDistribution parameters: %.5f\t%.5f\n",
              k, v$Name, v$Id, v$Package, v$DistributionType, v$Mean, v$Sd, v$Cov, v$x0, v$DistributionParameters[1], v$DistributionParameters[2]
            ))
            cat(sprintf("ParamType: %s\n", v$ParamType))
            cat("ParamValues:")
            cat(v$ParamValues)
            cat("\n____________________\n\n")
          } else if (inherits(v, "PARAM_DETVAR")) {
            cat("\n____________________")
            cat("\nDeterministic variable depending on a parameter variable\n")
            cat(sprintf("%d.\tName (ID): %s (%d)\n", k, v$Name, v$Id))
            cat("ParamValues:")
            cat(v$ParamValues)
            cat("\n____________________\n\n")
          } else {
            cat(sprintf(
              "%d.\tName (ID): %s (%d)\tPackage::DistributionType: %s::%s\tMean: %.3f\tSd: %.3f\tCov: %.3f\tx0: %.3f\tDistribution parameters: %.5f\t%.5f\n",
              k, v$Name, v$Id, v$Package, v$DistributionType, v$Mean, v$Sd, v$Cov, v$x0, v$DistributionParameters[1], v$DistributionParameters[2]
            ))
          }
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


      cat("--------------------- END -----------------")

    }
  )
)



#' @title Parametric Limit-State Function (PARAM_LSF)
#'
#' @description Identical to `SYS_LSF` but inherits from it to make the class hierarchy clear
#' when a parametric study is performed.  No additional fields or methods are
#' required.
#'
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#'
#' @export PARAM_LSF
#' @exportClass PARAM_LSF
#'
PARAM_LSF <- setRefClass(
  Class = "PARAM_LSF",
  contains = "SYS_LSF"
)
