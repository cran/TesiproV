# -----------------------------------------------------------------------------
#' @title JCSS hyperparameters for a structural basic variable
#'
#' @description Returns the hyperparameter vector of the Normal-Gamma or
#' LogNormal-Gamma prior distribution for a given structural basic variable,
#' as tabulated in the JCSS Probabilistic Model Code.
#'
#' @param basic.var Character. Structural basic variable identifier.
#'   Currently `"fc"` (concrete compressive strength) and `"fy"` (steel yield
#'   strength) are implemented.
#' @param for.value Character. For `basic.var = "fc"`: concrete class string
#'   (e.g. `"C25"`). For `basic.var = "fy"`: steel grade string (e.g. `"B500"`).
#' @param component.type Character. Production type. For `basic.var = "fc"`:
#'   `"ready-mixed"` or `"pre-cast"`. For `basic.var = "fy"`: `"structural steel"`
#'   or `"reinforcing steel"` (sic, matches internal spelling). Defaults to
#'   `"ready-mixed"`.
#'
#' @return A named numeric vector of hyperparameters. For `"fc"`:
#'   `c(m', s', n', v')` (Normal-Gamma / LogNormal-Gamma). For `"fy"`:
#'   `c(mu, sigma)`. Returns `NULL` with a warning if `basic.var` is not
#'   implemented.
#'
#' @details Internal helper, not exported. Used by [PROB_BASEVAR]'s
#'   `fromJCSS()` method to populate `DistributionParameters` from tabulated
#'   JCSS values.
#'
#' @noRd
.JCSS_distparam <- function(basic.var,
                            for.value,
                            component.type = "ready-mixed"
) {

  if (basic.var == "fc") {

    values_JCSS <- c(
      # ready-mixed  [m', s', n', v']
      c(3.40,  0.14,  3.0, 10),   # C15
      c(3.525, 0.13,  3.0, 10),   # C20
      c(3.65,  0.12,  3.0, 10),   # C25
      c(3.75,  0.105, 3.0, 10),   # C30
      c(3.85,  0.09,  3.0, 10),   # C35
      c(3.915, 0.08,  3.0, 10),   # C40
      c(3.98,  0.07,  3.0, 10),   # C45
      c(NA,    NA,    NA,  NA),   # C55
      # pre-cast
      c(NA,    NA,    NA,  NA),   # C15
      c(3.715, 0.095, 3.0, 10),   # C20
      c(3.80,  0.09,  3.0, 10),   # C25
      c(3.875, 0.085, 3.0, 10),   # C30
      c(3.95,  0.08,  3.0, 10),   # C35
      c(4.015, 0.075, 3.5, 10),   # C40
      c(4.08,  0.07,  4.0, 10),   # C45
      c(4.15,  0.05,  4.0, 10)    # C55
    )

    hyper_JCSS <- array(
      values_JCSS,
      dim      = c(4, 8, 2),
      dimnames = list(
        param          = c("m'", "s'", "n'", "v'"),
        fck_cyl        = c("C15","C20","C25","C30","C35","C40","C45","C55"),
        component.type = c("ready-mixed", "pre-cast")
      )
    )
    if (!for.value %in% dimnames(hyper_JCSS)$fck_cyl) {
      stop("Unknown concrete class: ", for.value)
    }

    DistributionParameter <- hyper_JCSS[, for.value, component.type]
    #names(DistributionParameter) <- c("m'", "s'", "n'", "v'")
    return(as.numeric(DistributionParameter))

  }
  else if (basic.var == "fy") {
    if (component.type == "structural steel") {
      stop("Not implemented yet.")
    }
    else if (component.type == "reinforcing steel") {
      values_JCSS <- c(
        c(420+2*30,  30),   # Grade 60
        c(520+2*30,  30),   # Grade 75
        c(500+2*30,  30),   # B500
        c(550+2*30,  30)    # B550
      )

      param_JCSS <- array(
        values_JCSS,
        dim      = c(2, 4, 1),
        dimnames = list(
          param          = c("mu", "sigma"),
          fy             = c("Grade 60","Grade 75","B500","B550"),
          component.type = c("reinforcing steel")
        )
      )
      if (!for.value %in% dimnames(param_JCSS)$fy) {
        stop("Unknown reinforcing steel grade: ", for.value)
      }

      DistributionParameter <- param_JCSS[, for.value, component.type]
    }
    #names(DistributionParameter) <- c("mu", "sigma")
    return(as.numeric(DistributionParameter))
  }
  else {
    warning("basic.var '", basic.var, "' not available - returning NULL")
    return(NULL)
  }
}

# -----------------------------------------------------------------------------
#' @title Distribution type for a JCSS basic variable
#'
#' @description Maps a structural basic variable (and, for `"fy"`, its
#' component type) to the corresponding TesiproV distribution type identifier.
#'
#' @param basic.var Character. Structural basic variable identifier (`"fc"`,
#'   `"fy"`, or any other value, which falls back to `"norm"`).
#' @param component.type Character. Production/component type, only relevant
#'   for `basic.var = "fy"` (`"structural steel"` -> `"lnorm"`,
#'   `"reinforcing steel"` -> `"norm"`).
#'
#' @return Character scalar: one of `"lt"`, `"lnorm"`, or `"norm"`.
#'
#' @details Internal helper, not exported. Used by [PROB_BASEVAR]'s
#'   `fromJCSS()` method to determine `DistributionType`.
#'
#' @noRd

.JCSS_disttype <- function(basic.var,component.type) {
  if (basic.var == "fc") {"lt"}
  else if (basic.var == "fy") {
    if (component.type == "structural steel") {"lnorm"}
    else if (component.type == "reinforcing steel") {"norm"}
  }
  else {"norm"}
}


# -----------------------------------------------------------------------------
#' @title Maximum-likelihood hyperparameter estimation (Rackwitz 1983)
#'
#' @description Estimates the four hyperparameters `(m', s', n', v')` of the
#' Normal-Gamma (or LogNormal-Gamma) conjugate prior from simulated datasets,
#' following the maximum-likelihood procedure of Rackwitz (1983). Used to
#' generate a Student-t (`"st"`) or Log-Student-t (`"lt"`) representation of a
#' normal or lognormal basic variable for Bayesian updating.
#'
#' @param mue Numeric. Mean of the underlying normal or lognormal variable.
#' @param sig Numeric. Standard deviation of the underlying normal or
#'   lognormal variable.
#' @param dist.type Character. `"norm"` or `"lnorm"`.
#' @param n.gen.num Integer. Sample size per simulated dataset.
#' @param n.gen.samp Integer. Number of simulated datasets (Monte Carlo
#'   repetitions) over which the ML estimates are averaged.
#'
#' @return A rounded numeric vector `c(m', s', n', v')` of length 4. Returns
#'   `NULL` with a warning if `dist.type` is neither `"norm"` nor `"lnorm"`.
#'
#' @details The fourth hyperparameter `v'` is obtained by solving
#'   `digamma(v/2) - h_snake - log(v / (2 * h_quer)) = 0` via [stats::uniroot()],
#'   with a grid-search fallback if root-finding fails. For `dist.type =
#'   "lnorm"`, the first two hyperparameters are additionally transformed from
#'   the linear to the logarithmic scale.
#'
#'   Internal helper, not exported. Used by [PROB_BASEVAR]'s `genHyperParam()`
#'   method.
#'
#' @references Rackwitz, R. (1983). *Predictive distribution of strength under
#'   control.* Materials and Structures, 16(4).
#'
#' @noRd

.MLE_distparam <- function(mue,
                           sig,
                           dist.type,
                           n.gen.num,
                           n.gen.samp) {

  if (dist.type =="norm") {
    # Simulate n.gen.samp datasets of size n.gen.num (vectorised)
    myData <- matrix(
      round(stats::rnorm(n.gen.samp * n.gen.num, mue, sig), 2),
      nrow = n.gen.samp, ncol = n.gen.num
    )
  } else if (dist.type =="lnorm") {
    sn <- sqrt(log(1 + (sig / mue)^2))
    mn <- log(mue / sqrt(1 + (sig / mue)^2))
    myData <- matrix(
      round(stats::rlnorm(n.gen.samp * n.gen.num, mn, sn), 2),
      nrow = n.gen.samp, ncol = n.gen.num
    )
  }


  # Maximum-Likelihood estimation (Rackwitz 1983)
  m_i       <- rowMeans(myData)
  h_i       <- 1 / apply(myData, 1, stats::var)
  log_h_i   <- log(h_i)
  hxm_i     <- h_i * m_i
  hxsqm_i   <- h_i * m_i^2

  h_quer   <- mean(h_i)
  h_snake  <- mean(log_h_i)
  h_hat    <- mean(hxm_i)
  h_revhat <- mean(hxsqm_i)

  hyper_prior    <- numeric(4)
  hyper_prior[1] <- h_hat / h_quer
  hyper_prior[2] <- h_quer^(-0.5)
  hyper_prior[3] <- (h_revhat - (h_hat^2) / h_quer)^(-1)

  # Solve for v' (truncated series, Rackwitz 1983)
  # Exact: digamma() + uniroot() instead of grid search
  funA <- function(v) digamma(v / 2) - h_snake - log(v/(2*h_quer))

  # Choose upper bound such that digamma(v_max/2) > h_snake
  v_max <- max(500, exp(h_snake) * 4)   # safe upper margin

  hyper_prior[4] <- tryCatch(
    stats::uniroot(funA, interval = c(1e-6, v_max))$root,
    error = function(e) {
      # Fall back to grid search if uniroot fails
      v_grid    <- seq(1e-6, v_max, by = 0.1)
      residuals <- abs(sapply(v_grid, function(v) digamma(v/2) - h_snake - log(v/(2*h_quer))))
      v_grid[which.min(residuals)]
    }
  )

  #hyper_prior[5] <- (log(h_quer)-h_snake)^(-1)


  if (dist.type == "lnorm") {
    s_lt <- sqrt(log((hyper_prior[2] / hyper_prior[1])^2 + 1))
    m_lt <- log(hyper_prior[1]) - s_lt^2 / 2
    hyper_prior[1] <- m_lt
    hyper_prior[2] <- s_lt
  } else if (dist.type != "norm") {
    warning("dist.type '", dist.type, "' not available - returning NULL")
    return(NULL)
  }

  hyper_prior <- round(hyper_prior, 2)

  #names(hyper_prior) <- c("m'", "s'", "n'", "v'")
  return(hyper_prior)
}

.default_package <- function(type) {
  switch(tolower(type),
         st = , lt = "TesiproV",
         gumbel    = "evd",
         "stats")
}

.resolve_package <- function(package, type) {
  # "", "stats" and "TesiproV" count as set automatically and follow the type.
  # Anything else (e.g. "extraDistr") was chosen deliberately by the user.
  if (is_missing(package) || package %in% c("stats", "TesiproV")) {
    .default_package(type)
  } else {
    package
  }
}

# Stops if a basic variable in lDistr is not completely defined, i.e. its
# distribution cannot be evaluated. Accepts funlists ($d, $p, $q, $r) and the
# wrapper format list(funlist, params) used by SYS_PROB.
.check_lDistr <- function(lDistr, method) {
  for (i in seq_along(lDistr)) {
    d <- lDistr[[i]]
    fun <- NULL
    if (is.list(d) && is.function(d$q)) {
      fun <- d
    } else if (is.list(d) && length(d) >= 1L && is.list(d[[1]]) &&
               is.function(d[[1]]$q)) {
      fun <- d[[1]]
    }
    med <- if (is.null(fun)) {
      NA_real_
    } else {
      tryCatch(suppressWarnings(fun$q(0.5)), error = function(e) NA_real_)
    }
    if (length(med) != 1L || !is.finite(med)) {
      nm <- if (is.null(fun) || is.null(fun$name) || !nzchar(fun$name)) {
        "<unnamed>"
      } else {
        fun$name
      }
      stop(sprintf(
        "[%s] Basic variable %d ('%s') is not completely defined. Please define Mean and Sd (or Cov) or DistributionParameters.",
        method, i, nm
      ), call. = FALSE)
    }
  }
  invisible(TRUE)
}