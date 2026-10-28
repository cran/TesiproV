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
#'                               functions (`"stats"`, `"evd"`, `"EnvStats"` or `"TesiproV"`).
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
    .cache = "ANY"
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

      # # --- Early integrity checks ---------------------------------------
      # # Case 1: completely empty object (no Mean, no Sd, no parameters)
      # if ((is.null(Mean) || length(Mean) == 0L) &&
      #     (is.null(Sd) || length(Sd) == 0L) &&
      #     (is.null(DistributionParameters) ||
      #      all(is.na(DistributionParameters)))) {
      #   warning(sprintf(
      #     "PROB_BASEVAR '%s' only useable with valid Mean/Sd or valid distribution parameters.",
      #     ifelse(is_empty(Name), "<unnamed>", Name)
      #   ), call. = FALSE)
      # }

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

      Package <<- .resolve_package(Package, DistributionType)

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

      if ((is_missing(Package) || Package == "") && (DistributionType == "st" || DistributionType == "lt")) {
        Package <<- "TesiproV"
      } else if (is_missing(Package) || Package == "") {
        Package <<- "stats"
      }

      ## the shifted lognormal is implemented in TesiproV itself,
      ## older scripts may still pass Package = "brms"
      if (DistributionType == "slnorm") Package <<- "TesiproV"

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
        DistributionParameters,
        Package
      ), algo = "xxhash64")

      ## 6)  return cached functions if they already exist
      if (exists(cache_key, envir = .cache, inherits = FALSE)) {
        cached <- get(cache_key, envir = .cache, inherits = FALSE)
        return(list(cached, DistributionParameters))
      }

      .native_args <- function(type, par) {
        par <- as.list(unname(par))
        switch(type,
               gamma   = list(shape = par[[1]], scale = par[[2]]),
               exp     = list(rate  = par[[1]]),
               weibull = list(shape = par[[1]], scale = par[[2]]),
               par
        )
      }

      ## 7)  otherwise build the four functions
      par <- DistributionParameters   # snapshot, used by all branches
      if (DistributionType == "emp") {
        d_fun <- function(x) {
          do.call("demp",
                  envir = loadNamespace("EnvStats"),
                  args = list(x = x, obs = par)
          )
        }
        p_fun <- function(x) {
          do.call("pemp",
                  envir = loadNamespace("EnvStats"),
                  args = list(q = x, obs = par)
          )
        }
        q_fun <- function(x) {
          do.call("qemp",
                  envir = loadNamespace("EnvStats"),
                  args = list(p = x, obs = par)
          )
        }
        r_fun <- function(x) {
          do.call("remp",
                  envir = loadNamespace("EnvStats"),
                  args = list(n = x, obs = par)
          )
        }
      }
      else if (DistributionType == "slnorm") {
        d_fun <- function(x) {
          do.call("dshifted_lnorm",
                  envir = loadNamespace("TesiproV"),
                  c(list(x), as.list(par))
          )
        }
        p_fun <- function(x) {
          do.call("pshifted_lnorm",
                  envir = loadNamespace("TesiproV"),
                  c(list(x), as.list(par))
          )
        }
        q_fun <- function(x) {
          do.call("qshifted_lnorm",
                  envir = loadNamespace("TesiproV"),
                  c(list(x), as.list(par))
          )
        }
        r_fun <- function(x) {
          do.call("rshifted_lnorm",
                  envir = loadNamespace("TesiproV"),
                  c(list(x), as.list(par))
          )
        }
      }
      else if (DistributionType == "st" || DistributionType == "lt") {
        ns_t <- asNamespace("TesiproV")
        f_d <- get(paste0("d", DistributionType), envir = ns_t)
        f_p <- get(paste0("p", DistributionType), envir = ns_t)
        f_q <- get(paste0("q", DistributionType), envir = ns_t)
        f_r <- get(paste0("r", DistributionType), envir = ns_t)
        d_fun <- function(x) f_d(x, par)
        p_fun <- function(x) f_p(x, par)
        q_fun <- function(x) f_q(x, par)
        r_fun <- function(x) f_r(x, par)
      }
      else {
        args <- .native_args(DistributionType, DistributionParameters)  # snapshot
        ns <- loadNamespace(Package)
        d_fun <- function(x) do.call(paste0("d", DistributionType), c(list(x), args), envir = ns)
        p_fun <- function(x) do.call(paste0("p", DistributionType), c(list(x), args), envir = ns)
        q_fun <- function(x) do.call(paste0("q", DistributionType), c(list(x), args), envir = ns)
        r_fun <- function(x) do.call(paste0("r", DistributionType), c(list(x), args), envir = ns)
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
      ##  Forward: Mean & Sd (Cov may be NA) -> native parameters
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
                 ## 1) Try numerical solution for k (shape)
                 root_fun <- function(k) {
                   sqrt(base::gamma(1 + 2 / k) / base::gamma(1 + 1 / k)^2 - 1) - Sd / Mean
                 }
                 ## 2) Choose interval (typically 0.1 - 100)
                 k_est <- tryCatch(
                   uniroot(root_fun, interval = c(0.1, 1000), tol = 1e-10)$root,
                   error = function(e) {
                     ## If the root search fails - fall back to the empirical formula
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
        else if (DistributionType == "st" && length(DistributionParameters) != 4) {
          stop("Student-t requires 4 hyperparameters: m', s', n', v'")
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

                 Mean <<- if (nue > 1) m else NA_real_
                 Sd <<- if (nue > 2) s / sqrt(n / (n + 1)) * sqrt(nue / (nue - 2)) else if (nue > 1 && nue <= 2) Inf else NA_real_
                 # Sd <<- s/sqrt(n/(n+1))/sqrt(nue/(nue-2))
                 Cov <<- if (nue > 2) Sd / Mean else NA_real_
               },
               lt = {
                 # Mean and Sd not defined for log-Student-t distribution
                 # Median will be used as Start Value for FORM
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
    },

    fromJCSS = function( BaseVar,
                         NominalValue,
                         ComponentType = "ready-mixed",...) {
      "Sets the distribution according to the JCSS Probabilistic Model Code
      and modifies this object. BaseVar = 'fc': concrete compressive strength,
      NominalValue is the strength class (e.g. 'C30'), ComponentType is
      'ready-mixed' or 'pre-cast'. The result is a log-Student-t distribution
      ('lt'). BaseVar = 'fy': yield strength, NominalValue is the steel grade
      (e.g. 'B500'), ComponentType is 'reinforcing steel'. The result is a
      normal distribution. Not all combinations are tabulated."

      dist.type <- .JCSS_disttype(BaseVar,ComponentType)

      params <- .JCSS_distparam(
        basic.var = BaseVar,
        for.value = NominalValue,
        component.type = ComponentType
      )

      Mean <<- NA_real_
      Sd <<- NA_real_
      Cov <<- NA_real_
      x0 <<- NA_real_
      DistributionType <<- dist.type
      DistributionParameters <<- params
      Package <<- .default_package(dist.type)

      prepare()

      invisible(.self)
    },

    genHyperParam = function(n.gen.num = 15,
                             n.gen.samp= 10000, ...
    ) {
    "Returns a new PROB_BASEVAR with a Student-t ('st') or log-Student-t
    ('lt') distribution for this normal or lognormal variable. The
    hyperparameters c(m', s', n', v') are estimated by maximum likelihood
    (Rackwitz 1983) from n.gen.samp simulated samples of size n.gen.num.
    The result depends on the random numbers, use set.seed() for
    reproducible values."

      if (DistributionType == "norm") {
        dist.type <- "st"
      }
      else if (DistributionType == "lnorm") {
        dist.type <- "lt"
      }
      else {
        stop("Generation of Hyperparameter of joint density only implemented for normal and lognormal distribution.",
             call. = FALSE)
      }

      params <- .MLE_distparam(
        mue            = Mean,
        sig            = Sd,
        dist.type      = DistributionType,
        n.gen.num      = n.gen.num,
        n.gen.samp     = n.gen.samp
      )

      obj <- PROB_BASEVAR(
        Id = Id,
        Name = Name,
        Description = Description,
        DistributionType = dist.type,
        DistributionParameters = params
      )

      return(obj)
    },

    updateBasevar = function(PROB_QUALCON, debug.level = 0, min_r2 = 0.99) {
     "Bayesian update with the OC curve of a PROB_QUALCON object. This
     variable must have DistributionType 'st' or 'lt' with the
     hyperparameters c(m', s', n', v'). Returns a new PROB_BASEVAR with the
     distribution after acceptance of the lot, this object is not changed.
     The outgoing distribution is fitted by the same distribution type. If
     the fit gives R^2 < min_r2, an empirical distribution ('emp') is
     returned instead. debug.level = 1 prints progress messages."

      # --- Input validation -------------------------------------------
      if (!is(PROB_QUALCON, "PROB_QUALCON"))
        stop("'PROB_QUALCON' must be a PROB_QUALCON object.", call. = FALSE)
      if (!DistributionType %in% c("st", "lt"))
        stop("Bayesian updating requires DistributionType 'st' or 'lt'.",
             call. = FALSE)
      if (isTRUE(PROB_QUALCON$isStale()))
        warning("PROB_QUALCON$var has changed since getOC() was last run - Pa_grid may be outdated.",
                call. = FALSE)

      # --- Ensure OC curve exists ------------------------------------
      if (all(is.na(PROB_QUALCON$Pa_grid))) {
        if (debug.level >= 1) message("OC data not found - running getOC() ...")
        PROB_QUALCON$getOC()
      }

      # --- Compute outgoing PDF --------------------------------------
      if (debug.level >= 1) message("Computing outgoing PDF ...")
      f_o_result <- .oc_get_outgoing_pdf(
        theta      = PROB_QUALCON$theta_grid,
        P_a        = PROB_QUALCON$Pa_grid,
        conf.type  = PROB_QUALCON$ConformityType,
        conf.param = PROB_QUALCON$ConformityParam,
        dist.type  = DistributionType,
        dist.param = DistributionParameters,
        lower      = PROB_QUALCON$DefectLimitLower,
        upper      = PROB_QUALCON$DefectLimitUpper
      )

      if (is.null(f_o_result))
        stop("Failed to compute outgoing PDF - check PROB_QUALCON parameters.",
             call. = FALSE)

      # --- Build updated copy ----------------------------------------
      var_post        <- .self$copy()
      var_post$field(".cache", new.env(parent = emptyenv()))

      hp_fit <- .oc_fit_posterior_hyper(f_o_result, dist.type = DistributionType, min_r2 = min_r2)

      probs <- f_o_result[[2]]
      probs[is.na(probs) | is.nan(probs) | probs < 0] <- 0
      probs <- probs / sum(probs)
      empirical_sample <- sample(f_o_result[[1]], size = 100000, replace = TRUE, prob = probs)

      if (!is.null(hp_fit)) {
        hp_post   <- hp_fit
        dist_type <- DistributionType
        if (debug.level >= 1)
          message(sprintf("Parametric fit accepted (R^2 = %.4f).", attr(hp_fit, "r2")))
      } else {
        hp_post   <- empirical_sample
        dist_type <- "emp"
        # the reason is reported by .oc_fit_posterior_hyper()
      }

      var_post$DistributionParameters <- hp_post
      var_post$DistributionType   <- dist_type
      var_post$Mean <- f_o_result$m_o
      var_post$Sd   <- f_o_result$s_o
      var_post$Cov <- if (isTRUE(var_post$Mean != 0))
        var_post$Sd / var_post$Mean else 0

      # Cache for chained updates
      assign("hyper_posterior", hp_post, envir = var_post$.cache)
      assign("hyper_prior",     DistributionParameters, envir = var_post$.cache)
      assign("prior_dist_type", DistributionType, envir = var_post$.cache)
      assign("f_o_result",   f_o_result, envir = var_post$.cache)
      assign("empirical_sample", empirical_sample, envir = var_post$.cache)
      if (!is.null(hp_fit)) {
        assign("fit_hp",        hp_fit,                    envir = var_post$.cache)
        assign("fit_dist_type", DistributionType,           envir = var_post$.cache)
        assign("fit_r2",        attr(hp_fit, "r2"),         envir = var_post$.cache)
      }

      # Prior moments: if they are not available analytically (NA for
      # hyperparameter-defined st / lt), approximate them numerically
      prior_mom <- c(Mean = Mean, Sd = Sd)
      prior_approx <- FALSE
      if (!is.finite(prior_mom[["Mean"]]) || !is.finite(prior_mom[["Sd"]])) {
        hp <- DistributionParameters
        if (DistributionType == "st" && hp[4] > 2) {
          # analytical moments exist for v' > 2
          prior_mom <- c(Mean = hp[1],
                         Sd = hp[2] * sqrt((hp[3] + 1) / hp[3]) * sqrt(hp[4] / (hp[4] - 2)))
        } else {
          num_mom <- .numeric_moments_hyper(DistributionType, hp)
          if (!is.null(num_mom)) {
            prior_mom <- num_mom
            prior_approx <- TRUE
          }
        }
      }

      assign("prior_mean", prior_mom[["Mean"]], envir = var_post$.cache)

      message(sprintf(
        "updateBasevar complete.\n  Prior     : Mean = %.4f  Sd = %.4f%s\n  Posterior : Mean = %.4f  Sd = %.4f",
        prior_mom[["Mean"]], prior_mom[["Sd"]],
        if (prior_approx) "  (numerical approximation)" else "",
        var_post$Mean, var_post$Sd
      ))
      return(var_post)
    },

    plot = function(x_range = NULL, n = 500, add = FALSE, compare = FALSE,
                    prior = FALSE, unit = NULL, ...) {
      "Plots the probability density function of this basic variable.
      If prior = TRUE and this object was produced by updateBasevar(), the
      density of the prior is drawn into the same plot (prior vs. posterior).
      unit is an optional character string with the unit of the variable
      (e.g. \"mm\"). It is added to the x axis label as [unit].
      If compare = TRUE and this object was produced by updateBasevar(),
      the fitted parametric curve (if a fit was attempted) and the raw
      empirical resample density are overlaid for visual comparison."

      ld <- getlDistr()
      d_fun <- ld[[1]]$d
      q_fun <- ld[[1]]$q

      if (is.null(x_range)) x_range <- c(q_fun(0.001), q_fun(0.999))

      has_prior <- prior && !is.null(.cache) &&
        exists("hyper_prior", envir = .cache, inherits = FALSE) &&
        exists("prior_dist_type", envir = .cache, inherits = FALSE)

      if (prior && !has_prior) {
        message("No prior data found in the cache - this instance was not created by updateBasevar(). Argument prior is ignored.")
      }

      if (has_prior) {
        hp_prior   <- get("hyper_prior", envir = .cache, inherits = FALSE)
        type_prior <- get("prior_dist_type", envir = .cache, inherits = FALSE)
        d_prior <- if (type_prior == "st") dst else dlt
        q_prior <- if (type_prior == "st") qst else qlt
        x_range <- range(c(x_range, q_prior(c(0.001, 0.999), hp_prior)))
      }

      has_cache_data <- compare && !is.null(.cache) &&
        exists("empirical_sample", envir = .cache, inherits = FALSE)

      if (compare && !has_cache_data) {
        message("No comparison data found in the cache - this instance was probably not created by updateBasevar(). Argument compare is ignored.")
      }

      if (has_cache_data) {
        emp <- get("empirical_sample", envir = .cache, inherits = FALSE)
        dens_emp <- stats::density(emp, from = x_range[1], to = x_range[2])
        x_range <- range(c(x_range, dens_emp$x))
      }

      x_vals <- seq(x_range[1], x_range[2], length.out = n)
      y_vals <- d_fun(x_vals)
      if (has_prior) y_prior <- d_prior(x_vals, hp_prior)

      # Title: "Var. <Id>: <Name>", subtitles: description and distribution type.
      # The name is used as label of the x axis.
      name_str <- if (is_empty(Name)) "PROB_BASEVAR" else Name
      id_ok    <- length(Id) == 1L && !is.na(Id)
      title_str <- if (id_ok) sprintf("Var. %s: %s", format(Id), name_str) else name_str
      desc_str <- if (length(Description) == 1L && !is.na(Description) &&
                      nzchar(trimws(Description))) Description else NULL
      dist_str <- sprintf("Distribution type: %s", .distribution_label(DistributionType))

      # colours: with a prior the prior is black and the posterior coloured
      post_col <- "dodgerblue3"
      cur_col  <- if (has_prior) post_col else "black"

      if (add) {
        graphics::lines(x_vals, y_vals, ...)
      } else {
        y_max <- max(y_vals, if (has_cache_data) dens_emp$y, if (has_prior) y_prior, na.rm = TRUE)
        dots <- list(...)
        auto_title <- is.null(dots$main)
        if (auto_title) {
          # extra space on top for title and subtitles
          op <- graphics::par(mar = graphics::par("mar") + c(0, 0, if (is.null(desc_str)) 1.5 else 2.5, 0))
          on.exit(graphics::par(op), add = TRUE)
        }
        has_unit <- is.character(unit) && length(unit) == 1L && !is.na(unit) && nzchar(trimws(unit))
        x_lab <- if (is_empty(Name)) "x" else Name
        if (has_unit) x_lab <- sprintf("%s [%s]", x_lab, unit)
        plot_args <- list(x = x_vals, y = y_vals, type = "l", lwd = 2, col = cur_col,
                          xlab = x_lab,
                          ylab = "Density [-]",
                          ylim = c(0, y_max * 1.05), main = "")
        plot_args[names(dots)] <- dots
        do.call(graphics::plot, plot_args)
        if (auto_title) {
          if (is.null(desc_str)) {
            graphics::title(main = title_str, line = 2.2)
            graphics::mtext(dist_str, side = 3, line = 0.8, cex = 0.9, col = "grey30")
          } else {
            graphics::title(main = title_str, line = 3.2)
            graphics::mtext(desc_str, side = 3, line = 1.8, cex = 0.9, col = "grey30")
            graphics::mtext(dist_str, side = 3, line = 0.6, cex = 0.9, col = "grey30")
          }
        }
      }

      legend_labels <- character(0); legend_lty <- integer(0); legend_col <- character(0)

      if (!add) {
        legend_labels <- c(legend_labels,
                           sprintf("%s (%s)", if (has_prior) "Posterior" else "PDF", DistributionType))
        legend_lty <- c(legend_lty, 1); legend_col <- c(legend_col, cur_col)
      }

      if (has_prior) {
        graphics::lines(x_vals, y_prior, lty = 1, lwd = 2, col = "black")
        legend_labels <- c(legend_labels, sprintf("Prior (%s)", type_prior))
        legend_lty <- c(legend_lty, 1); legend_col <- c(legend_col, "black")
      }

      # vertical lines at the mean (current distribution and prior)
      mean_col <- if (has_prior) post_col else "darkorange3"
      if (is.finite(Mean) && Mean >= x_range[1] && Mean <= x_range[2]) {
        graphics::abline(v = Mean, lty = 4, lwd = 2, col = mean_col)
        legend_labels <- c(legend_labels, sprintf("%s = %s", if (has_prior) "Posterior mean" else "Mean", format(signif(Mean, 5))))
        legend_lty <- c(legend_lty, 4); legend_col <- c(legend_col, mean_col)
      }
      if (has_prior && exists("prior_mean", envir = .cache, inherits = FALSE)) {
        pm <- get("prior_mean", envir = .cache, inherits = FALSE)
        if (is.finite(pm) && pm >= x_range[1] && pm <= x_range[2]) {
          graphics::abline(v = pm, lty = 4, lwd = 2, col = "black")
          legend_labels <- c(legend_labels, sprintf("Prior mean = %s", format(signif(pm, 5))))
          legend_lty <- c(legend_lty, 4); legend_col <- c(legend_col, "black")
        }
      }

      if (has_cache_data) {
        graphics::lines(dens_emp$x, dens_emp$y, lty = 2, lwd = 2, col = "grey40")
        legend_labels <- c(legend_labels, "empirical (resample)")
        legend_lty <- c(legend_lty, 2); legend_col <- c(legend_col, "grey40")

        if (exists("fit_hp", envir = .cache, inherits = FALSE)) {
          fit_hp   <- get("fit_hp", envir = .cache, inherits = FALSE)
          fit_type <- get("fit_dist_type", envir = .cache, inherits = FALSE)
          fit_r2   <- get("fit_r2", envir = .cache, inherits = FALSE)
          fit_d    <- if (fit_type == "st") dst else dlt

          graphics::lines(x_vals, fit_d(x_vals, fit_hp), lty = 3, lwd = 2, col = "firebrick")
          legend_labels <- c(legend_labels, sprintf("Fit (R^2 = %.3f)", fit_r2))
          legend_lty <- c(legend_lty, 3); legend_col <- c(legend_col, "firebrick")
        }
      }

      if (!add && length(legend_labels) > 1) {
        graphics::legend("topright", legend = legend_labels, lty = legend_lty,
                         lwd = 2, col = legend_col, bty = "n", cex = 0.85)
      }

      invisible(.self)
    }
  )
)

#' Full name of a distribution type
#'
#' Translates the short distribution identifier used in `PROB_BASEVAR`
#' (e.g. "norm") into a readable name (e.g. "Normal") for plot titles.
#' Unknown identifiers are returned unchanged.
#'
#' @param type Character, distribution identifier
#' @return Character, readable name
#' @keywords internal
#' @noRd
.distribution_label <- function(type) {
  labels <- c(
    norm = "Normal", lnorm = "Log-normal", slnorm = "Shifted log-normal",
    gumbel = "Gumbel", gamma = "Gamma", exp = "Exponential", beta = "Beta",
    weibull = "Weibull", st = "Student-t", lt = "Log-Student-t",
    binom = "Binomial", unif = "Uniform", emp = "Empirical"
  )
  if (length(type) == 1L && type %in% names(labels)) unname(labels[[type]]) else type
}

#' Numerical mean and standard deviation of a st / lt distribution
#'
#' Approximates the first two moments of a Student-t ("st") or log-Student-t
#' ("lt") distribution given by its hyperparameters c(m', s', n', v') by
#' trapezoidal integration of the density between the 1e-6 and 1 - 1e-6
#' quantiles. The result is a truncated approximation, which is used where the
#' analytical moments are undefined (e.g. for "lt").
#'
#' @param dist.type "st" or "lt"
#' @param hyper numeric vector c(m', s', n', v')
#' @return named numeric vector c(Mean, Sd) or NULL if the approximation fails
#' @keywords internal
#' @noRd
.numeric_moments_hyper <- function(dist.type, hyper) {
  if (!dist.type %in% c("st", "lt")) return(NULL)
  d_fun <- if (dist.type == "st") dst else dlt
  q_fun <- if (dist.type == "st") qst else qlt
  res <- tryCatch({
    lims <- q_fun(c(1e-6, 1 - 1e-6), hyper)
    if (!all(is.finite(lims)) || lims[1] >= lims[2]) return(NULL)
    x <- seq(lims[1], lims[2], length.out = 20001)
    f <- d_fun(x, hyper)
    f[!is.finite(f) | f < 0] <- 0
    dx <- diff(x)
    trap <- function(y) sum((y[-1] + y[-length(y)]) / 2 * dx)
    mass <- trap(f)
    if (!is.finite(mass) || mass <= 0) return(NULL)
    m <- trap(x * f) / mass
    s <- sqrt(max(trap((x - m)^2 * f) / mass, 0))
    c(Mean = m, Sd = s)
  }, error = function(e) NULL)
  res
}

#' @title Object to store a deterministic parameters (PROB_DETVAR)
#'
#' @description A convenience subclass of `PROB_BASEVAR` that represents a deterministic
#' value (i.e. a variable with practically zero variance).  The object is
#' automatically converted to a normal distribution with an infinitesimal
#' standard deviation (`Sd = max(|Value| / 1e7, 1e-12)`).
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
  fields = list(
    Value = "numeric"
  ),
  methods = list(
    initialize = function(...) {
      callSuper(...)
      DistributionType <<- "norm" # Not strictly necessary, since PROB_BASEVAR sets "norm" for DistributionType = ""
      if (length(Value) != 1) {
        stop("PROB_DETVAR: 'Value' must be a scalar numeric.")
      }

      Value <<- as.numeric(Value)[1]
      Mean <<- Value
      # positive standard deviation, also for Mean <= 0
      Sd <<- max(abs(Mean) / 1e7, 1e-12)
      Cov <<- if (Mean == 0) 0 else Sd / abs(Mean)

      prepare()
    }
  )
)
