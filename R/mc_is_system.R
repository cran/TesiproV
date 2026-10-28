# =====================================================================
# MC_IS_system.R - Importance Sampling for a system of LSFs
# =====================================================================

#' Monte-Carlo Importance Sampling for a system of limit-state functions
#'
#' @inheritParams MC_IS_single
#' @param sys_type  Either "parallel" (default) or "serial"
#' @param beta_l  Optional threshold for reliability index screening.
#'   Limit-state functions with \eqn{\beta > beta_l} may be excluded
#'   from mixture construction in system analysis.
#'
#' @details
#' The failure probability is estimated with the standard importance sampling
#' estimator \eqn{\hat{p}_f = n^{-1} \sum_i I_i w_i} with
#' \eqn{w = \varphi(u) / h(u)}, where \eqn{h} is a mixture of unit-variance
#' normal densities centred at the proposal points. \eqn{I_i = 1} if
#' \eqn{g_k(x_i) < 0} for at least one component (serial) or for all
#' components (parallel), otherwise \eqn{I_i = 0}. The coefficient of
#' variation is derived from the sample variance of \eqn{I w}.
#' \itemize{
#'   \item Serial systems: the component design points (FORM) are the
#'     centres. If they are close and nearly parallel, a single shift at their
#'     mean is used.
#'   \item Parallel systems: the design point of the joint failure domain
#'     \eqn{\max_i g_i < 0} is used as the main centre (weight 0.7). The
#'     component design points are kept as a defensive part of the mixture.
#' }
#' The returned \code{ESS} is the effective number of failure samples
#' \eqn{(\sum I w)^2 / \sum (I w)^2}. The argument \code{stability_mode}
#' no longer influences the result because all sums are accumulated in the
#' log domain.
#'
#' @keywords internal
#'
MC_IS_system <- function(lsf,
                         lDistr,
                         cov_user = 0.05,
                         n_batch,
                         n_max,
                         use_threads,
                         sys_type,
                         dataRecord = TRUE,
                         beta_l,
                         densityType,
                         dps = NULL,
                         debug.level = 0,
                         streams = NULL,
                         libPaths_local = .libPaths(),
                         seed = NULL,
                         backend = "future",
                         adaptive_alpha = FALSE,
                         alpha_update_rate = 0.1,
                         adaptive_batch = FALSE,
                         batch_control = list(),
                         min_adapt_samples = 10000,
                         alpha_min = 0.02,
                         stability_mode = c("robust", "fast")) {
  if (debug.level >= 1) {
    cat("DEBUG LEVEL INSIDE MC_IS_system:", debug.level, "\n")
  }

  # Basic input checks
  stability_mode <- match.arg(stability_mode)

  if (missing(debug.level) || is.null(debug.level)) {
    debug.level <- 0
  }

  if (!is.list(lsf) || length(lsf) == 0L) {
    stop("MC_IS_system: 'lsf' must be a non-empty list of functions.")
  }

  if (!is.list(lDistr) || length(lDistr) != length(lsf)) {
    stop("MC_IS_system: 'lDistr' must be a list with one entry per LSF.")
  }

  if (!sys_type %in% c("serial", "parallel")) {
    stop("MC_IS_system: 'sys_type' must be 'serial' or 'parallel'.")
  }

  # basic variables must be completely defined (one list per LSF)
  for (k in seq_along(lDistr)) .check_lDistr(lDistr[[k]], "MC_IS")
  
  debug.TAG <- "MC_IS_System"
  if (debug.level >= 1) message("[MC_IS_system] starting system importance sampling...")

  tic <- proc.time()


  # defaults for adaptive_batch, user values override single entries
  batch_control <- utils::modifyList(
    list(n_min = 2000, n_max = 100000, K_future = 5),
    batch_control
  )

  n_lsfs <- length(lsf)

  n_sim <- 0L

  ESS_global <- 0
  pf_hat <- NA_real_
  var_est <- NA_real_
  cov <- Inf

  # --------------------------------------------------
  # Normalize LSF input: accept functions or SYS_LSF objects
  # --------------------------------------------------
  lsf_fun <- lapply(lsf, function(obj) {
    if (is.function(obj)) {
      return(obj)
    } else if (!is.null(obj$func) && is.function(obj$func)) {
      return(obj$func)
    } else {
      stop("Invalid LSF format: must be function or SYS_LSF object.")
    }
  })


  # set seeds
  # if (!is.null(seed)) {
  #   RNGkind("L'Ecuyer-CMRG")
  #   set.seed(as.integer(seed))
  # }
  rng_mgr <- create_rng_manager(seed)


  ## -------------------------------------------------
  ## Cluster setup for the Windows parallel backend
  ##
  ## On Unix-like systems, run_parallel() uses mclapply().
  ## Therefore, no persistent cluster is required there.
  ## -------------------------------------------------

  cl <- NULL

  if (
    backend == "parallel" &&
      use_threads > 1L &&
      .Platform$OS.type == "windows"
  ) {
    cl <- make_parallel_cluster(
      use_threads = use_threads,
      seed = seed
    )

    on.exit(parallel::stopCluster(cl), add = TRUE)
  }


  ## ------------------------------------------------------------------
  ## Normalize lDistr -> ensure each marginal is a funlist (with $d/$p/$q/$r)
  ## ------------------------------------------------------------------
  normalize_marginal <- function(m) {
    if (is.list(m) && !is.null(m$p) && is.function(m$p)) {
      return(m)
    } # already funlist
    if (is.list(m) && length(m) >= 1 && is.list(m[[1]]) && !is.null(m[[1]]$p) && is.function(m[[1]]$p)) {
      return(m[[1]])
    }
    stop("MC_IS_system: lDistr contains unexpected marginal format", call. = FALSE)
  }

  lDistr <- lapply(lDistr, function(per_lsf) {
    lapply(per_lsf, normalize_marginal)
  })

  ## ------------------------------------------------------------------
  ## User supplied design points: one vector per LSF (physical space).
  ## A numeric vector with one value per LSF is accepted if every LSF
  ## has exactly one variable.
  ## ------------------------------------------------------------------
  if (!is.null(dps)) {
    n_var_lsf <- vapply(lDistr, length, integer(1))
    if (is.numeric(dps) && length(dps) == n_lsfs && all(n_var_lsf == 1L)) {
      dps <- as.list(dps)
    }
    dps_ok <- is.list(dps) && length(dps) == n_lsfs &&
      all(vapply(seq_len(n_lsfs), function(i) {
        is.numeric(dps[[i]]) && length(dps[[i]]) == n_var_lsf[i] &&
          all(is.finite(dps[[i]]))
      }, logical(1)))
    if (!dps_ok) {
      stop("MC_IS_system: 'dps' must be a list with one design point vector per limit-state function (one value per variable of that function).",
           call. = FALSE)
    }
  }

  ## ------------------------------------------------------------------
  ## 1) Build unique variable list and mappings
  ## ------------------------------------------------------------------

  var.names <- character(0)
  vars <- list() # store funlist per unique variable
  var.map <- integer(0)
  var.per.lsf <- integer(n_lsfs)

  for (i in seq_len(n_lsfs)) {
    nd <- length(lDistr[[i]])
    var.per.lsf[i] <- nd

    for (j in seq_len(nd)) {
      funlist <- lDistr[[i]][[j]]
      if (!is.list(funlist) || is.null(funlist$p) || !is.function(funlist$p)) {
        stop(sprintf("MC_IS_system: lDistr[[%d]][[%d]] is not a valid marginal (funlist)", i, j), call. = FALSE)
      }
      name_j <- funlist$name
      if (is.null(name_j) || name_j == "") name_j <- paste0("var_", length(var.names) + 1L)
      idx <- match(name_j, var.names, nomatch = 0L)
      if (idx == 0L) {
        var.names <- c(var.names, name_j)
        vars[[length(var.names)]] <- funlist
        var.map <- c(var.map, length(var.names))
      } else {
        var.map <- c(var.map, idx)
      }
    }
  }

  n_vars_unique <- length(vars)

  # var.in.lsf
  var.in.lsf <- vector("list", n_lsfs)
  k <- 1L
  for (i in seq_len(n_lsfs)) {
    var.in.lsf[[i]] <- var.map[k:(k + var.per.lsf[i] - 1L)]
    k <- k + var.per.lsf[i]
  }

  # lsf.in.var
  lsf.in.var <- vector("list", n_vars_unique)
  for (p in seq_len(n_vars_unique)) lsf.in.var[[p]] <- which(vapply(var.in.lsf, function(v) p %in% v, logical(1)))

  ## ------------------------------------------------------------------
  ## 2) FORM analysis -> design points per LSF
  ## ------------------------------------------------------------------
  res_form <- vector("list", n_lsfs)
  res_form.beta <- rep(NA_real_, n_lsfs)
  res_form.dp <- vector("list", n_lsfs)
  dp_u <- vector("list", n_lsfs)

  for (i in seq_len(n_lsfs)) {
      if (!is.null(dps)) {
      # --- user supplied design point (one value per variable) ---
      x_star <- as.numeric(dps[[i]])

      dp_u[[i]] <- vapply(seq_along(x_star), function(k) {
        Fxi <- tryCatch(lDistr[[i]][[k]]$p(x_star[k]),
          error = function(e) NA_real_
        )
        qnorm(min(max(as.numeric(Fxi), 1e-12), 1 - 1e-12))
      }, numeric(1))

      res_form[[i]] <- list(
        beta = NA_real_,
        pf = NA_real_,
        x_points = x_star
      )

      res_form.beta[i] <- NA_real_
      res_form.dp[[i]] <- x_star
    } else {
      # --- Extract funlists ---
      lDistr_funlist_i <- lapply(lDistr[[i]], function(dd) {
        if (is.list(dd) && !is.null(dd$p) && is.function(dd$p)) {
          return(dd)
        }
        if (is.list(dd) && length(dd) >= 1 &&
          is.list(dd[[1]]) && !is.null(dd[[1]]$p)) {
          return(dd[[1]])
        }
        stop("Invalid marginal format for FORM")
      })

      # --- Build wrapper ---
      lDistr_for_FORM_i <- lapply(lDistr_funlist_i, function(fun) {
        params <- c(NA_real_, NA_real_)
        if (!is.null(fun$mean)) params[1] <- as.numeric(fun$mean)
        if (!is.null(fun$Sd)) params[2] <- as.numeric(fun$Sd)
        list(fun, params)
      })

      # --- Call FORM ---
      res_form[[i]] <- TesiproV::FORM(
        lsf_fun[[i]],
        lDistr_for_FORM_i,
        n_optim = 20,
        loctol  = 1e-4
      )

      res_form.beta[i] <- res_form[[i]]$beta
      res_form.dp[[i]] <- as.numeric(unlist(res_form[[i]]$x_points))

      # --- U-Transformation ---
      x_star_vec <- res_form.dp[[i]]

      dp_u[[i]] <- numeric(length(lDistr_funlist_i))

      for (kidx in seq_along(lDistr_funlist_i)) {
        xsk <- as.numeric(x_star_vec[kidx])

        Fxi <- tryCatch(
          lDistr_funlist_i[[kidx]]$p(xsk),
          error = function(e) NA_real_
        )

        Fxi <- min(max(as.numeric(Fxi), 1e-12), 1 - 1e-12)

        dp_u[[i]][kidx] <- qnorm(Fxi)
      }
    }
  }


  # ------------------------------------------------------------
  # Automatic vectorization detection for each LSF
  # ------------------------------------------------------------

  lsf_vectorized_list <- logical(n_lsfs)

  for (k in seq_len(n_lsfs)) {
    # small test matrix (2 samples, correct variable subset)
    idxs <- var.in.lsf[[k]]

    x_star <- res_form.dp[[k]]

    test_mat <- matrix(
      rep(x_star, 2),
      nrow = 2,
      byrow = TRUE
    )

    test_result <- try(lsf_fun[[k]](test_mat), silent = TRUE)

    if (!inherits(test_result, "try-error") &&
      is.numeric(test_result) &&
      length(test_result) == 2 &&
      all(is.finite(test_result))) {
      lsf_vectorized_list[k] <- TRUE

      if (debug.level >= 1) {
        message(sprintf("[MC_IS_system] LSF %d detected as vectorized.", k))
      }
    }
  }

  # --------------------------------------------------
  # Global embedding of the design points in U-space
  # --------------------------------------------------

  dp_u_global <- vector("list", n_lsfs)

  for (k in seq_len(n_lsfs)) {
    u_star <- numeric(n_vars_unique)
    idxs <- var.in.lsf[[k]]
    u_star[idxs] <- dp_u[[k]]
    dp_u_global[[k]] <- u_star
  }

  if (debug.level >= 2) {
    dp_glob_str <- paste(
      sapply(seq_along(dp_u_global), function(i) {
        paste0(
          "LSF", i, ": (",
          paste(round(dp_u_global[[i]], 4), collapse = ", "),
          ")"
        )
      }),
      collapse = " | "
    )

    message(sprintf("[MC_IS_system] Global U-shifts: %s", dp_glob_str))
  }


  # --------------------------------------------------
  # Geometric analysis of failure regions (U-space)
  # --------------------------------------------------

  compute_geometry <- function(dp_u_list, res_form_list) {
    n_lsfs <- length(dp_u_list)

    if (n_lsfs < 2) {
      return(list(multimodal = FALSE, min_cos = NA_real_, max_dist = 0))
    }

    # gradient directions embedded in the global variable vector
    # (each LSF only uses its own subset of the variables)
    alpha_glob <- vector("list", n_lsfs)
    for (i in seq_len(n_lsfs)) {
      dy_i <- res_form_list[[i]]$dy
      if (is.null(dy_i) || length(dy_i) != length(var.in.lsf[[i]]) ||
        any(!is.finite(dy_i)) || sum(dy_i^2) == 0) {
        # no usable gradient information (e.g. user supplied design points)
        return(list(multimodal = TRUE, min_cos = NA_real_, max_dist = NA_real_))
      }
      a <- numeric(n_vars_unique)
      a[var.in.lsf[[i]]] <- dy_i
      alpha_glob[[i]] <- a / sqrt(sum(a^2))
    }

    min_cos <- 1
    max_dist <- 0

    for (i in 1:(n_lsfs - 1)) {
      for (j in (i + 1):n_lsfs) {
        # --- Distance between design points ---
        dist_ij <- sqrt(sum((dp_u_list[[i]] - dp_u_list[[j]])^2))

        max_dist <- max(max_dist, dist_ij)

        # --- Angle between gradients (alpha vectors) ---
        cos_ij <- sum(alpha_glob[[i]] * alpha_glob[[j]])

        min_cos <- min(min_cos, cos_ij)
      }
    }

    # single shift only if all design points are close and nearly parallel
    multimodal_flag <- !(min_cos > 0.85 && max_dist < 1.5)

    return(list(
      multimodal = multimodal_flag,
      min_cos = min_cos,
      max_dist = max_dist
    ))
  }

  ## decision multimodal distribution or single shift
  geom_info <- compute_geometry(dp_u_global, res_form)

  use_multimodal <- geom_info$multimodal

  if (debug.level >= 1) {
    message(sprintf(
      "[MC_IS_system] Geometry: min cos=%.3f, max dist=%.3f -> multimodal=%s",
      geom_info$min_cos,
      geom_info$max_dist,
      use_multimodal
    ))
  }


  # ------------------------------------------------------------------
  # Parallel system: system design point
  #
  # A parallel system fails if ALL component limit states are violated,
  # i.e. if max_i g_i < 0. The point of the joint failure domain closest to
  # the origin is used as an additional proposal centre. The mixture of the
  # component design points alone may miss the joint failure domain.
  # ------------------------------------------------------------------
  dp_u_sys <- NULL

  if (sys_type == "parallel") {
    g_sys <- function(x) {
      vals <- numeric(n_lsfs)
      for (i in seq_len(n_lsfs)) {
        vals[i] <- lsf_fun[[i]](x[var.in.lsf[[i]]])
      }
      max(vals)
    }

    lDistr_global <- lapply(vars, function(fun) {
      list(fun, c(fun$mean, fun$Sd))
    })

    res_form_sys <- tryCatch(
      suppressWarnings(TesiproV::FORM(
        g_sys,
        lDistr_global,
        n_optim = 20,
        loctol = 1e-4
      )),
      error = function(e) NULL
    )

    x_star_sys <- if (is.null(res_form_sys)) NULL else as.numeric(unlist(res_form_sys$x_points))

    if (!is.null(x_star_sys) && length(x_star_sys) == n_vars_unique &&
      all(is.finite(x_star_sys))) {
      u_sys <- numeric(n_vars_unique)
      for (j in seq_len(n_vars_unique)) {
        Fxi <- tryCatch(vars[[j]]$p(x_star_sys[j]), error = function(e) NA_real_)
        Fxi <- min(max(as.numeric(Fxi), 1e-12), 1 - 1e-12)
        u_sys[j] <- qnorm(Fxi)
      }
      if (all(is.finite(u_sys))) dp_u_sys <- u_sys
    }

    if (is.null(dp_u_sys) && debug.level >= 1) {
      message("[MC_IS_system] System design point not available - using component design points only.")
    }
  }

  # proposal centres (modes) and mixture weights
  mode_shifts <- dp_u_global
  if (sys_type == "parallel") {
    use_multimodal <- TRUE
    if (!is.null(dp_u_sys)) {
      # defensive mixture: system design point plus the component design points
      mode_shifts <- c(list(dp_u_sys), dp_u_global)
    }
  }
  n_modes <- length(mode_shifts)

  # ------------------------------------------------------------------
  # Initial mixture weights (serial only)
  # ------------------------------------------------------------------
  if (sys_type == "serial") {
    pf_single <- pnorm(-res_form.beta)
    # alpha_i <- pf_single / sum(pf_single)
    alpha_i <- rep(1 / n_lsfs, n_lsfs)
  }

  ## ------------------------------------------------------------------
  ## 3) Build dp matrix (global)
  ## ------------------------------------------------------------------
  dp <- matrix(NA_real_, nrow = n_vars_unique, ncol = n_lsfs)
  for (i in seq_len(n_lsfs)) {
    local_vars <- var.in.lsf[[i]]
    if (!is.null(res_form.dp[[i]])) {
      for (kl in seq_along(local_vars)) {
        global_id <- local_vars[kl]
        if (kl <= length(res_form.dp[[i]])) dp[global_id, i] <- res_form.dp[[i]][kl]
      }
    }
  }

  dp.means <- rowMeans(dp, na.rm = TRUE)

  ## ------------------------------------------------------------------
  ## 4) ai weights
  ## ------------------------------------------------------------------
  ai <- matrix(0, nrow = n_vars_unique, ncol = n_lsfs)
  for (p in seq_len(n_vars_unique)) {
    comps <- lsf.in.var[[p]]
    if (length(comps) == 0L) next
    betas <- res_form.beta[comps]
    betas[!is.finite(betas) | betas <= 0] <- Inf
    inv_b <- 1 / betas
    s <- sum(inv_b, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) ai[p, comps] <- 1 / length(comps) else ai[p, comps] <- inv_b / s
  }

  ## ------------------------------------------------------------------
  ## 5) Build hv mixtures per global variable
  ## ------------------------------------------------------------------
  hv_list <- vector("list", n_vars_unique)
  for (p in seq_len(n_vars_unique)) {
    comps <- lsf.in.var[[p]]
    if (length(comps) == 0L) {
      hv_list[[p]] <- list(rfun = function(n) rep(NA_real_, n), dfun = function(x) rep(0, length(x)))
      next
    }
    orig_funlist <- vars[[p]] # funlist (d,p,q,r,mean,Sd,...)
    sd_p <- orig_funlist$Sd
    comp_dfun <- list()
    comp_rfun <- list()
    comp_probs <- numeric(length(comps))
    for (ci in seq_along(comps)) {
      j <- comps[ci]
      mu_j <- dp[p, j]
      sd_loc <- ifelse(is.finite(sd_p) && sd_p > 0, sd_p, max(1e-8, abs(mu_j) * 0.1, 1e-3))
      dist_type <- tolower(if (is.null(orig_funlist$DistributionType) || orig_funlist$DistributionType == "") "norm" else orig_funlist$DistributionType)
      if (dist_type %in% c("lnorm", "lt", "weibull")) {
        if (!is.finite(mu_j) || mu_j <= 0) {
          mn <- mu_j
          sn <- sd_loc
          comp_rfun[[ci]] <- function(n, m = mn, s = sn) stats::rnorm(n, mean = m, sd = s)
          comp_dfun[[ci]] <- function(x, m = mn, s = sn) stats::dnorm(x, mean = m, sd = s)
        } else {
          sd_rel <- sd_loc / mu_j
          sn <- sqrt(log(1 + sd_rel^2))
          mn <- log(mu_j / sqrt(1 + sd_rel^2))
          comp_rfun[[ci]] <- function(n, mn = mn, sn = sn) stats::rlnorm(n, meanlog = mn, sdlog = sn)
          comp_dfun[[ci]] <- function(x, mn = mn, sn = sn) stats::dlnorm(x, meanlog = mn, sdlog = sn)
        }
      } else {
        comp_rfun[[ci]] <- function(n, m = mu_j, s = sd_loc) stats::rnorm(n, mean = m, sd = s)
        comp_dfun[[ci]] <- function(x, m = mu_j, s = sd_loc) stats::dnorm(x, mean = m, sd = s)
      }
      comp_probs[ci] <- ai[p, j]
    }
    if (!is.finite(sum(comp_probs)) || sum(comp_probs) <= 0) comp_probs <- rep(1 / length(comp_probs), length(comp_probs)) else comp_probs <- comp_probs / sum(comp_probs)
    hv_dfun <- function(x) {
      res <- numeric(length(x))
      for (ci in seq_along(comp_dfun)) res <- res + comp_probs[ci] * comp_dfun[[ci]](x)
      res
    }
    hv_rfun <- function(n) {
      if (n <= 0) {
        return(numeric(0))
      }
      comps_drawn <- sample(seq_along(comp_probs), size = n, replace = TRUE, prob = comp_probs)
      res <- numeric(n)
      for (ci in unique(comps_drawn)) {
        idx <- which(comps_drawn == ci)
        res[idx] <- comp_rfun[[ci]](length(idx))
      }
      res
    }
    hv_list[[p]] <- list(rfun = hv_rfun, dfun = hv_dfun, comps = comps, probs = comp_probs)
  }

  ## ------------------------------------------------------------------
  ## 6) Monte Carlo Loop (single-sample per iter)
  ## ------------------------------------------------------------------
  I_sum <- 0L
  pf_acc <- 0
  var_acc <- 0
  n_sim <- 0L
  cov <- Inf
  log_pf <- -Inf
  log_var <- NA_real_
  log_cov <- Inf


  if (dataRecord) {
    max_iter_est <- ceiling(n_max)
    rec_pf <- numeric(max_iter_est)
    rec_var <- numeric(max_iter_est)
    rec_cov <- numeric(max_iter_est)
    rec_nsim <- integer(max_iter_est)
    rec_I_sum <- integer(max_iter_est)
    rec_time_user <- numeric(max_iter_est)
    rec_time_sys <- numeric(max_iter_est)
    rec_time_elapsed <- numeric(max_iter_est)
  }

  orig_d_list <- lapply(vars, function(v) v$d)


  # --------------------------------------------------------------
  # Monte Carlo Loop (worker-invariant, SNIS, global stabilization)
  # --------------------------------------------------------------

  # --------------------------------------------------------------
  # Monte Carlo Loop - Multimodal, worker-invariant, global log-stable
  # --------------------------------------------------------------

  if (!exists("dp_u")) {
    stop("dp_u not constructed before Monte Carlo loop")
  }

  n_sim <- 0L

  ESS_modes <- NULL


  total_batch <- n_batch
  iter <- 0L

  if (sys_type == "serial") {
    if (all(is.finite(res_form.beta))) {
      pf_single <- pnorm(-res_form.beta)
      alpha_i <- pf_single / sum(pf_single)
    } else {
      alpha_i <- rep(1 / n_lsfs, n_lsfs)
    }
  }

  if (sys_type == "serial") {
    mix_w <- alpha_i
  } else if (!is.null(dp_u_sys)) {
    mix_w <- c(0.7, rep(0.3 / n_lsfs, n_lsfs))
  } else {
    mix_w <- rep(1 / n_lsfs, n_lsfs)
  }
  mix_w <- mix_w / sum(mix_w)

  # running sums in the log domain
  state <- list(
    log_sum_Iw = -Inf,
    log_sum_Iw2 = -Inf,
    n_sim = 0L
  )

  t_sample <- NA_real_
  t_parallel <- NA_real_
  while (TRUE) {
    iter <- iter + 1L
    total_batch <- n_batch

    if (debug.level >= 1) {
      cat("total_batch:", total_batch, "\n")
    }

    # ==========================================================
    # Sampling Block
    # ==========================================================

    # if (sys_type == "serial") {
    #   # Serial: multimodal sampling
    #   mode_ids <- sample(seq_len(n_lsfs),
    #     size = total_batch,
    #     replace = TRUE,
    #     prob = alpha_i
    #   )

    #   u_mat <- matrix(0, nrow = total_batch, ncol = n_vars_unique)

    #   for (i in seq_len(total_batch)) {
    #     shift <- dp_u_global[[mode_ids[i]]]
    #     u_mat[i, ] <- rnorm(n_vars_unique,
    #       mean = shift,
    #       sd = 1
    #     )
    #   }
    # } else { # parallel

    #   # Parallel: single shift at system design point
    #   u_mat <- matrix(rnorm(total_batch * n_vars_unique),
    #     ncol = n_vars_unique
    #   )

    #   u_mat <- sweep(u_mat, 2, dp_u_sys, "+")

    #   # Dummy mode ids (needed for ESS calculation)
    #   mode_ids <- rep(1L, total_batch)
    # }

    ## -------------------------------------------------
    ## Deterministic sampling (master only)
    ## -------------------------------------------------
    if (debug.level >= 1) {
      t_sample_start <- proc.time()
    }

    u_vec <- rng_mgr$generate_norm_block(total_batch * (n_vars_unique + 1L))
    u_mat <- matrix(u_vec[seq_len(total_batch * n_vars_unique)], ncol = n_vars_unique)

    # Draw the mixture mode per sample (same RNG stream as the normals)
    v_mode <- stats::pnorm(u_vec[total_batch * n_vars_unique + seq_len(total_batch)])
    mode_ids <- findInterval(v_mode, c(0, cumsum(mix_w)[-n_modes]))

    shift_mat <- do.call(rbind, mode_shifts) # n_modes x n_vars_unique
    if (use_multimodal) {
      u_mat <- u_mat + shift_mat[mode_ids, , drop = FALSE]
    } else {
      mode_ids <- rep(1L, total_batch)
      u_mat <- sweep(u_mat, 2, colMeans(shift_mat), "+")
    }
    # if (sys_type == "serial") {
    #   ## deterministic mode selection
    #   mode_probs <- alpha_i
    #   mode_ids <- sample(
    #     seq_len(n_lsfs),
    #     size = total_batch,
    #     replace = TRUE,
    #     prob = mode_probs
    #   )

    #   for (i in seq_len(total_batch)) {
    #     shift <- dp_u_global[[mode_ids[i]]]
    #     u_mat[i, ] <- u_mat[i, ] + shift
    #   }
    # } else {
    #   u_mat <- sweep(u_mat, 2, dp_u_sys, "+")
    #   mode_ids <- rep(1L, total_batch)
    # }

    # ==========================================================
    # Parallel Evaluation Block
    # ==========================================================

    # chunk_ids <- rep(seq_len(use_threads),
    #   length.out = total_batch
    # )

    # chunk_indices <- split(seq_len(total_batch), chunk_ids)
    if (use_threads <= 1L) {
      chunk_indices <- list(seq_len(total_batch))
    } else {
      chunk_indices <- split(
        seq_len(total_batch),
        cut(seq_len(total_batch), breaks = use_threads)
      )
    }

    worker_fun <- function(idx) {
      u_block <- u_mat[idx, , drop = FALSE]
      n_local <- nrow(u_block)

      # ----------------------------
      # 1) Transformation to physical space and evaluation of the LSFs
      # ----------------------------
      fail_count <- integer(n_local)

      for (k in seq_len(n_lsfs)) {
        idxs <- var.in.lsf[[k]]
        n_loc_vars <- length(idxs)

        x_sub <- matrix(NA_real_, n_local, n_loc_vars)

        for (j in seq_len(n_loc_vars)) {
          v <- idxs[j]

          p <- pnorm(u_block[, v])
          p <- pmin(pmax(p, 1e-14), 1 - 1e-14)

          x_sub[, j] <- vars[[v]]$q(p)
        }

        fun_k <- lsf_fun[[k]]

        if (lsf_vectorized_list[k]) {
          vals <- fun_k(x_sub)
        } else {
          vals <- apply(x_sub, 1, fun_k)
        }

        fail_count <- fail_count + (vals < 0)
      }

      if (sys_type == "serial") {
        I_vals <- as.integer(fail_count > 0)
      } else {
        I_vals <- as.integer(fail_count == n_lsfs)
      }

      # ----------------------------
      # 2) log phi(u)
      # ----------------------------
      log_phi <- rowSums(dnorm(u_block, log = TRUE))

      # ----------------------------
      # 3) log h(u): same density as used for sampling
      # ----------------------------
      if (use_multimodal) {
        # Mixture of normal distributions at the design points, shares mix_w
        log_h_mat <- matrix(NA_real_, n_local, n_modes)

        for (k in seq_len(n_modes)) {
          shifted <- sweep(u_block, 2, mode_shifts[[k]], "-")

          log_h_mat[, k] <- log(mix_w[k]) +
            rowSums(dnorm(shifted, log = TRUE))
        }

        max_log <- apply(log_h_mat, 1, max)

        log_h <- max_log +
          log(rowSums(exp(log_h_mat - max_log)))
      } else {
        # Single normal distribution at the mean of the design points
        shift <- rowMeans(do.call(cbind, dp_u_global))
        shifted <- sweep(u_block, 2, shift, "-")

        log_h <- rowSums(dnorm(shifted, log = TRUE))
      }

      log_w <- log_phi - log_h

      list(
        log_w = log_w,
        I_vals = I_vals,
        mode_ids = mode_ids[idx]
      )
    }


    # --------------------------------------------------
    # Parallel evaluation
    # --------------------------------------------------

    # if (backend == "future") {
    #   res_list <- future.apply::future_lapply(
    #     chunk_indices,
    #     worker_fun,
    #     future.seed = FALSE
    #   )
    # } else {
    #   res_list <- lapply(chunk_indices, worker_fun)
    # }
    if (debug.level >= 1) {
      t_parallel_start <- proc.time()
    }

    res_list <- run_parallel(
      X = chunk_indices,
      FUN = worker_fun,
      backend = backend,
      use_threads = use_threads,
      cl = cl
    )

    if (debug.level >= 1) {
      t_parallel <- (proc.time() - t_parallel_start)[3]
    }

    if (debug.level >= 1) {
      cat(sprintf(
        "[MC_IS_system]: Sampling: %.3f s | Parallel: %.3f s\n",
        t_sample, t_parallel
      ))
    }

    # ==========================================================
    # Weight Accumulation Block (refactored)
    # ==========================================================

    all_log_w <- unlist(lapply(res_list, `[[`, "log_w"))
    all_I_vals <- unlist(lapply(res_list, `[[`, "I_vals"))
    # all_I_matrix <- do.call(rbind, lapply(res_list, `[[`, "I_matrix"))
    all_mode_ids <- unlist(lapply(res_list, `[[`, "mode_ids"))

    if (debug.level >= 1) {
      cat("len(all_log_w):", length(all_log_w), "\n")
    }

    state <- compute_weight_update(
      all_log_w,
      all_I_vals,
      stability_mode,
      sys_type,
      state
    )

    n_sim <- state$n_sim

    # ==========================================================
    # Estimator Update Block
    # ==========================================================

    if (sys_type == "serial" && use_multimodal && stability_mode == "fast") {
      max_log_w <- max(all_log_w)
      shifted_weights <- exp(all_log_w - max_log_w)

      ESS_modes <- numeric(n_lsfs)

      for (k in seq_len(n_lsfs)) {
        idx_k <- which(all_mode_ids == k)

        if (length(idx_k) > 0) {
          w_k <- shifted_weights[idx_k]
          ESS_modes[k] <- (sum(w_k)^2) / sum(w_k^2)
        } else {
          ESS_modes[k] <- 0
        }
      }
    }
    # --------------------------------------------------
    # Adaptive batch size update
    # --------------------------------------------------
    if (adaptive_batch && ESS_global > 0) {
      RSE_target <- cov_user
      K_future <- batch_control$K_future
      n_min_batch <- batch_control$n_min
      n_max_batch <- batch_control$n_max

      # Estimate current efficiency factor
      kappa_hat <- ESS_global / n_sim

      if (kappa_hat > 0) {
        ESS_target <- 1 / (RSE_target^2)
        N_target <- ESS_target / kappa_hat
        N_remaining <- max(N_target - n_sim, 0)

        n_batch_new <- ceiling(N_remaining / K_future)

        # Ensure divisibility by threads
        n_batch_new <- ceiling(n_batch_new / use_threads) * use_threads

        # Apply bounds
        n_batch <- max(n_min_batch, min(n_batch_new, n_max_batch))

        if (debug.level >= 2) {
          message(sprintf(
            "[MC_IS_system] Adaptive batch updated to %d (ESS=%.2f)",
            n_batch, ESS_global
          ))
        }
      }
    }


    # ==========================================================
    # Estimator Update Block
    # ==========================================================

    est <- .mc_is_log_estimates(
      log_sum_Iw  = state$log_sum_Iw,
      log_sum_Iw2 = state$log_sum_Iw2,
      n           = n_sim
    )

    log_pf <- est$log_pf
    pf_hat <- est$pf
    log_var <- est$log_var
    var_est <- est$var
    log_cov <- est$log_cov
    cov <- est$cov

    # System ESS: effective number of failure samples.
    if (is.finite(state$log_sum_Iw) &&
      is.finite(state$log_sum_Iw2)) {
      log_ess <- 2 * state$log_sum_Iw - state$log_sum_Iw2
      ESS_global <- min(exp(log_ess), n_sim)
    } else {
      ESS_global <- 0
    }

    RSE_est <- cov


    if (adaptive_alpha && sys_type == "serial" && use_multimodal && n_sim > min_adapt_samples) {
      stop("adaptive_alpha currently requires I_matrix; disable adaptive_alpha or re-enable I_matrix.")

      # pf_i_hat <- numeric(n_lsfs)
      #
      # if (stability_mode == "fast") {
      #   pf_sys_hat <- sum(shifted_weights[all_I_vals == 1])
      #
      #   if (!is.na(pf_sys_hat) && is.finite(pf_sys_hat) && pf_sys_hat > 0) {
      #     for (k in seq_len(n_lsfs)) {
      #       pf_i_hat[k] <-
      #         sum((all_I_matrix[, k] * all_I_vals) *
      #           shifted_weights) / pf_sys_hat
      #     }
      #   }
      # }
      #
      # if (stability_mode == "robust") {
      #   # pf_sys_hat in log-domain
      #   pf_sys_hat <- exp(log_total_sum_Iw - log_total_sum_w)
      #
      #   if (!is.na(pf_sys_hat) && is.finite(pf_sys_hat) && pf_sys_hat > 0) {
      #     for (k in seq_len(n_lsfs)) {
      #       idx_k <- which(all_I_matrix[, k] == 1 &
      #         all_I_vals == 1)
      #
      #       if (length(idx_k) > 0) {
      #         log_contrib <- all_log_w[idx_k]
      #
      #         max_log_c <- max(log_contrib)
      #
      #         pf_i_hat[k] <-
      #           exp(
      #             max_log_c +
      #               log(sum(exp(log_contrib - max_log_c))) -
      #               log_total_sum_w
      #           )
      #       } else {
      #         pf_i_hat[k] <- 0
      #       }
      #     }
      #   }
      # }
      #
      # # update alpha
      # if (sum(pf_i_hat) > 0) {
      #   alpha_target <- pf_i_hat / sum(pf_i_hat)
      #
      #   alpha_i <- (1 - alpha_update_rate) * alpha_i +
      #     alpha_update_rate * alpha_target
      #
      #   alpha_i <- pmax(alpha_i, alpha_min)
      #   alpha_i <- alpha_i / sum(alpha_i)
      # }
    }


    if (debug.level >= 2) {
      if (adaptive_alpha) {
        print(alpha_i)
      }
    }


    # ==========================================================
    # Convergence Check
    # ==========================================================
    if (!exists("ESS_global")) ESS_global <- 0

    ESS_min <- if (sys_type == "parallel") 50 else 30
    n_min <- 5 * n_batch

    if (debug.level >= 1) {
      cat(
        "Iteration:", iter,
        "n_sim:", n_sim,
        "cov:", cov,
        "log_cov:", log_cov,
        "\n"
      )

      cat(sprintf(
        "[MC_IS_system]: ESS: %.3f \n",
        ESS_global
      ))
    }

    # Compare the CoV target on the log scale. This avoids making
    # the stopping decision from an underflowed or overflowed `cov`.
    cov_target_met <-
      is.numeric(cov_user) &&
        length(cov_user) == 1L &&
        is.finite(cov_user) &&
        cov_user > 0 &&
        !is.na(log_cov) &&
        log_cov < log(cov_user)

    stop_condition <-
      cov_target_met &&
        ESS_global > ESS_min &&
        n_sim > n_min

    if (stop_condition) {
      break
    }

    if (n_sim >= n_max) {
      break
    }
  }

  beta_final <- if (is.finite(log_pf) && log_pf <= 0) {
    -stats::qnorm(log_pf, log.p = TRUE)
  } else {
    NA_real_
  }

  out <- list(
    method = "MCIS_System",
    sys_type = sys_type,
    beta = beta_final,
    pf = pf_hat,
    log_pf = log_pf,
    var = var_est,
    log_var = log_var,
    cov_mc = cov,
    RSE_est = RSE_est,
    ESS = ESS_global,
    ESS_modes = if (sys_type == "serial") ESS_modes else NULL,
    n_mc = n_sim
  )

  if (sys_type == "parallel" && ESS_global < ESS_min) {
    warning(sprintf(
      "Parallel system: effective number of failure samples (ESS = %.1f) is low. Relative standard error ~= %.2f%%.",
      ESS_global,
      100 * RSE_est
    ))
  }

  if (debug.level >= 1) message("[MC_IS_system] finished.")

  return(out)
}


#' Monte Carlo Importance Sampling using smooth system aggregation
#'
#' Performs reliability estimation for systems of limit-state functions
#' by constructing a **smooth global system limit-state function**
#' and applying the existing single-LSF importance sampling algorithm.
#'
#' The system limit-state function is built using \code{build_gsys()},
#' which aggregates component limit states via smooth min/max operators.
#'
#' This approach converts a system reliability problem into a standard
#' single limit-state reliability problem.
#'
#' @param sys_input List of \code{SYS_LSF} objects defining the system.
#'
#' @param sys_structure System aggregation type: \code{"serial"},
#' \code{"parallel"}, or a custom function.
#'
#' @param smooth_kappa Smoothing parameter controlling the
#' approximation of min/max aggregation.
#'
#' @inheritParams MC_IS
#'
#' @return Result object returned by \code{MC_IS_single()} with additional
#' metadata describing the system configuration.
#'
#' @keywords internal
#'
MC_IS_gsys <- function(sys_input,
                       sys_structure = "serial",
                       smooth_kappa = 20,
                       cov_user = 0.05,
                       n_batch = 5000,
                       n_max = 1e6,
                       use_threads = 4,
                       backend = "future",
                       dataRecord = TRUE,
                       densityType = "norm",
                       dps = NULL,
                       debug.level = 0,
                       seed = NULL,
                       adaptive_alpha = FALSE,
                       alpha_update_rate = 0.1,
                       adaptive_batch = FALSE,
                       batch_control = list(),
                       min_adapt_samples = 10000,
                       alpha_min = 0.02,
                       stability_mode = "robust") {
  if (debug.level >= 1) {
    message("[MC_IS_gsys] building smooth system limit-state function...")
  }

  # -------------------------------------------------
  # Build system representation
  # -------------------------------------------------

  sys_obj <- build_gsys(
    sys_input = sys_input,
    sys_structure = sys_structure,
    smooth_kappa = smooth_kappa,
    debug.level = debug.level
  )

  # debug
  # print(body(sys_obj$g_sys))

  g_sys <- sys_obj$g_sys
  distr_flat <- sys_obj$distr_flat

  if (debug.level >= 1) {
    message(sprintf(
      "[MC_IS_gsys] system has %d LSFs and %d variables",
      sys_obj$n_lsfs,
      sys_obj$n_vars
    ))
  }

  # debug
  # message("kappa=", smooth_kappa)

  # -------------------------------------------------
  # Run importance sampling using existing routine
  # -------------------------------------------------

  res <- MC_IS_single(
    lsf = g_sys,
    lDistr = distr_flat,
    cov_user = cov_user,
    n_batch = n_batch,
    n_max = n_max,
    use_threads = use_threads,
    backend = backend,
    dataRecord = dataRecord,
    densityType = densityType,
    dps = dps,
    debug.level = debug.level,
    seed = seed,
    adaptive_alpha = adaptive_alpha,
    alpha_update_rate = alpha_update_rate,
    adaptive_batch = adaptive_batch,
    batch_control = batch_control,
    min_adapt_samples = min_adapt_samples,
    alpha_min = alpha_min,
    stability_mode = stability_mode
  )

  # -------------------------------------------------
  # Attach system metadata
  # -------------------------------------------------

  res$method <- "MCIS_System_gsys"
  res$sys_structure <- sys_structure
  res$smooth_kappa <- smooth_kappa
  res$n_lsfs <- sys_obj$n_lsfs
  res$n_vars <- sys_obj$n_vars

  if (debug.level >= 1) {
    message("[MC_IS_gsys] finished.")
  }

  return(res)
}


#' @title Internal weight accumulation for Monte Carlo Importance Sampling (System Case)
#'
#' @description
#' Internal helper used by `MC_IS_system()` to accumulate the importance
#' sampling sums that are needed for the standard estimator
#' \deqn{\hat{p}_f = \frac{1}{n}\sum_i I_i w_i, \qquad
#'   \widehat{\mathrm{Var}}(\hat{p}_f) =
#'   \frac{1}{n}\left(\frac{1}{n}\sum_i I_i w_i^2 - \hat{p}_f^2\right).}
#'
#' Both sums are accumulated in the log domain (log-sum-exp), so that very
#' large or very small weights cannot overflow or underflow.
#'
#' @param all_log_w Numeric vector of log-weights for the current batch.
#' @param all_I_vals Numeric vector of failure indicators (0/1).
#' @param stability_mode,sys_type Not used any more. Kept for compatibility
#'   of the call signature.
#' @param state List containing the current accumulated sums:
#'   \describe{
#'     \item{log_sum_Iw}{Log of the cumulative sum of the failure weights.}
#'     \item{log_sum_Iw2}{Log of the cumulative sum of the squared failure weights.}
#'     \item{n_sim}{Total number of samples processed.}
#'   }
#'
#' @return Updated `state` list.
#'
#' @details
#' This function is internal and not intended to be called directly by users.
#'
#' @keywords internal
compute_weight_update <- function(all_log_w,
                                  all_I_vals,
                                  stability_mode = NULL,
                                  sys_type = NULL,
                                  state) {
  logsumexp <- function(x) {
    if (length(x) == 0L) {
      return(-Inf)
    }
    m <- max(x)
    if (!is.finite(m)) {
      return(m)
    }
    m + log(sum(exp(x - m)))
  }

  log_add <- function(a, b) {
    m <- max(a, b)
    if (!is.finite(m)) {
      return(m)
    }
    m + log(exp(a - m) + exp(b - m))
  }

  lw_fail <- all_log_w[which(all_I_vals == 1)]

  state$log_sum_Iw <- log_add(state$log_sum_Iw, logsumexp(lw_fail))
  state$log_sum_Iw2 <- log_add(state$log_sum_Iw2, logsumexp(2 * lw_fail))
  state$n_sim <- state$n_sim + length(all_log_w)

  state
}


#' Build a smooth system limit-state function (g_sys)
#'
#' Constructs a **global system limit-state function** from a list of
#' individual limit-state functions (`SYS_LSF` objects). The resulting
#' function aggregates all component limit states using a **smooth
#' approximation** of the minimum or maximum operator.
#'
#' The smooth aggregation improves numerical stability and differentiability,
#' which is beneficial for gradient-based reliability methods such as
#' FORM, SORM, and importance sampling algorithms.
#'
#' For a serial system the aggregation is
#'
#' \deqn{
#' g_{\kappa}(x) =
#' -\frac{1}{\kappa}\log
#' \left(
#' \sum_{i=1}^{m} \exp(-\kappa g_i(x))
#' \right)
#' }
#'
#' which approximates
#'
#' \deqn{
#' \min_i g_i(x).
#' }
#'
#' For a parallel system the aggregation is
#'
#' \deqn{
#' g_{\kappa}(x) =
#' \frac{1}{\kappa}\log
#' \left(
#' \sum_{i=1}^{m} \exp(\kappa g_i(x))
#' \right)
#' }
#'
#' which approximates
#'
#' \deqn{
#' \max_i g_i(x).
#' }
#'
#' The approximation becomes exact for
#' \eqn{\kappa \rightarrow \infty}.
#'
#' In addition to the system function, this helper also constructs a
#' **flattened list of marginal distributions** and a **variable index map**
#' that links the global variable vector to the variables of each
#' individual limit-state function.
#'
#' @param sys_input List of \code{SYS_LSF} objects defining the component
#'   limit-state functions of the system.
#'
#' @param sys_structure System aggregation structure.
#'   Supported options are
#'   \describe{
#'     \item{\code{"serial"}}{System fails if any component fails.}
#'     \item{\code{"parallel"}}{System fails only if all components fail.}
#'     \item{\code{function(g_vals)}}{Custom aggregation function applied
#'     to the vector of component limit-state values.}
#'   }
#'
#' @param smooth_kappa Positive numeric value controlling the smoothness
#'   of the aggregation. Larger values approach the exact \code{min}
#'   or \code{max} operator.
#'
#' @param debug.level Optional verbosity level for diagnostic output.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{g_sys}}{Callable system limit-state function
#'   \code{function(x)} operating on the global variable vector.}
#'
#'   \item{\code{distr_flat}}{Flattened list of marginal distribution
#'   objects corresponding to all variables in the system.}
#'
#'   \item{\code{var_map}}{List mapping each limit-state function to the
#'   indices of its variables in the global vector.}
#'
#'   \item{\code{n_vars}}{Total number of basic variables in the system.}
#'
#'   \item{\code{n_lsfs}}{Number of limit-state functions.}
#' }
#'
#' @details
#' This function is primarily used internally by reliability algorithms
#' that operate on a **single limit-state function**, such as
#' importance sampling or crude Monte Carlo simulation.
#'
#' By constructing a smooth global limit-state function, the system
#' reliability problem can be treated as a standard single-LSF problem,
#' enabling the reuse of algorithms originally developed for
#' component reliability.
#'
#' @keywords internal
#'
#' @seealso
#' \code{smooth_min}, \code{smooth_max},
#' \code{\link{MC_IS}}, \code{\link{SYS_LSF}}
#'
#' @author
#' (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -
#' Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#'
build_gsys <- function(sys_input,
                       sys_structure = "serial",
                       smooth_kappa = 20,
                       debug.level = 0) {
  n_lsfs <- length(sys_input)

  is_custom <- is.function(sys_structure)

  if (n_lsfs < 1) {
    stop("build_gsys: no limit state functions provided.")
  }

  # -------------------------------------------------
  # 1. Extract LSF functions
  # -------------------------------------------------

  lsfs <- lapply(sys_input, function(obj) obj$getLSF())

  # -------------------------------------------------
  # 2. Build distribution list + variable mapping
  # -------------------------------------------------

  distr_flat <- list()
  var_map <- list()

  idx <- 1

  for (i in seq_along(sys_input)) {
    vars_i <- sys_input[[i]]$vars
    n_vars_i <- length(vars_i)

    var_map[[i]] <- idx:(idx + n_vars_i - 1)

    subdistr <- lapply(vars_i, function(v) {
      d <- v$getlDistr()

      if (is.list(d) && !is.null(d[[1]]$p)) {
        return(d[[1]])
      }

      if (is.list(d) && !is.null(d$p)) {
        return(d)
      }

      stop("build_gsys: invalid distribution object.")
    })

    distr_flat <- c(distr_flat, subdistr)

    idx <- idx + n_vars_i
  }

  n_vars <- length(distr_flat)


  # eval_lsfs <- function(x) {
  #   if (is.vector(x)) {
  #     g_vals <- numeric(n_lsfs)

  #     for (i in seq_len(n_lsfs)) {
  #       idxs <- var_map[[i]]

  #       xi <- as.numeric(x[idxs])

  #       g_vals[i] <- lsfs[[i]](xi)
  #     }

  #     return(g_vals)
  #   }

  #   if (is.matrix(x)) {
  #     n_samples <- nrow(x)

  #     G <- matrix(NA_real_, nrow = n_samples, ncol = n_lsfs)

  #     for (i in seq_len(n_lsfs)) {
  #       idxs <- var_map[[i]]

  #       Xi <- x[, idxs, drop = FALSE]

  #       G[, i] <- lsfs[[i]](Xi)
  #     }

  #     return(G)
  #   }

  #   stop("Unsupported input type in eval_lsfs")
  # }
  eval_lsfs <- function(x) {
    if (is.vector(x)) {
      g_vals <- numeric(n_lsfs)

      for (i in seq_len(n_lsfs)) {
        xi <- as.numeric(x[var_map[[i]]])

        g_vals[i] <- as.numeric(lsfs[[i]](xi))
      }

      return(g_vals)
    }

    if (is.matrix(x)) {
      n_samples <- nrow(x)

      G <- matrix(NA_real_, n_samples, n_lsfs)

      for (i in seq_len(n_lsfs)) {
        Xi <- x[, var_map[[i]], drop = FALSE] # n x k matrix
        G[, i] <- apply(Xi, 1, function(row) as.numeric(lsfs[[i]](row)))
      }

      return(G)
    }

    stop("Unsupported input type in eval_lsfs")
  }

  # -------------------------------------------------
  # 3. Build system LSF
  # -------------------------------------------------

  g_sys <- function(x) {
    G <- eval_lsfs(x)

    if (is.null(dim(G))) {
      if (is_custom) {
        return(sys_structure(G))
      }

      if (sys_structure == "serial") {
        m <- min(G)

        return(
          m - log(sum(exp(-smooth_kappa * (G - m)))) / smooth_kappa
        )
      }

      if (sys_structure == "parallel") {
        m <- max(G)

        return(
          m + log(sum(exp(smooth_kappa * (G - m)))) / smooth_kappa
        )
      }

      stop("Unsupported system structure.")
    }

    if (is.matrix(G)) {
      if (is_custom) {
        return(apply(G, 1, sys_structure))
      }

      if (sys_structure == "serial") {
        # m <- G[, 1]

        # if (ncol(G) > 1) {
        #   for (j in 2:ncol(G)) {
        #     m <- pmin(m, G[, j])
        #   }
        # }

        m <- G[, 1]
        for (j in seq_len(ncol(G) - 1) + 1) {
          m <- pmin(m, G[, j])
        }

        E <- exp(pmax(-smooth_kappa * (G - m), -700))
        S <- rowSums(E)

        S <- pmax(S, .Machine$double.xmin)

        res <- m - log(S) / smooth_kappa

        bad <- !is.finite(res)
        if (any(bad)) res[bad] <- m[bad]

        return(res)
      }

      if (sys_structure == "parallel") {
        # m <- G[, 1]

        # if (ncol(G) > 1) {
        #   for (j in 2:ncol(G)) {
        #     m <- pmax(m, G[, j])
        #   }
        # }
        m <- G[, 1]
        for (j in seq_len(ncol(G) - 1) + 1) {
          m <- pmax(m, G[, j])
        }

        E <- exp(pmin(smooth_kappa * (G - m), 700))
        S <- rowSums(E)

        S <- pmax(S, .Machine$double.xmin)

        res <- m + log(S) / smooth_kappa

        bad <- !is.finite(res)
        if (any(bad)) res[bad] <- m[bad]

        return(res)
      }

      stop("Unsupported system structure.")
    }

    stop("Unsupported system structure.")
  }

  # -------------------------------------------------
  # 4. Return structure
  # -------------------------------------------------

  list(
    g_sys = g_sys,
    distr_flat = distr_flat,
    var_map = var_map,
    n_vars = n_vars,
    n_lsfs = n_lsfs
  )
}
