#' @name FORM
#' @title First-Order Reliability Method (FORM)
#' @description
#' First-Order Reliability Method (FORM) for the approximation
#' of failure probabilities in structural reliability analysis.
#'
#' The FORM estimates the probability of failure by transforming
#' the basic random variables into standard normal space
#' and approximating the limit-state function by a first-order
#' (linear) Taylor expansion at the design point.
#'
#' The reliability index is defined as the minimum distance
#' from the origin to the limit-state surface in standard normal space.
#' The probability of failure is then approximated using
#' the standard normal cumulative distribution function.
#'
#' Failure is defined as \eqn{g(x) < 0}. The reliability index is positive if
#' the origin of the standard normal space lies in the safe domain and
#' negative otherwise. The probability of failure is \eqn{p_f = \Phi(-\beta)}.
#'
#' Solver selection: With an explicitly given \code{optim_type} exactly this
#' algorithm is run and nothing else. If it does not converge, a warning is
#' issued and the last iterate is returned with \code{converged = FALSE}.
#' With \code{optim_type = "auto"} (default) a robust chain is used: SLSQP
#' (\pkg{nloptr}), and iHLRF as a fallback. A warning is only issued if both fail.
#'
#' @param lsf Objective function representing the limit-state, e.g. \code{function(R,E){R-E}}.
#'   Supplied automatically by a SYS_ object - do not provide manually.
#' @param lDistr List of distribution objects created by TesiproV. Supplied automatically by
#'   a SYS_ object - do not provide manually.
#' @param n_optim Maximum number of iterations. For \code{"slsqp"} it is the
#'   maximum number of function evaluations of the optimizer. In \code{"auto"}
#'   SLSQP may use \code{max(200, 10 * n_optim)} evaluations.
#' @param loctol Relative tolerance. Converged means
#'   \eqn{|g(u)| \le loctol \cdot g_{scale}} and the distance of \eqn{u} from the
#'   line through the origin along the gradient is
#'   \eqn{\le loctol \cdot \max(1, \|u\|)}. SLSQP results are accepted with
#'   at most \code{max(loctol, 1e-5)}.
#' @param optim_type Solver:
#'   \itemize{
#'     \item \code{"auto"}     - (default) SLSQP, then iHLRF if needed
#'     \item \code{"ihlrf"}    - improved HLRF with merit-function line-search
#'                                (Zhang & Der Kiureghian 1997), no fallback
#'     \item \code{"rackfies"} - Rackwitz-Fiessler iteration (classic HLRF), no fallback
#'     \item \code{"slsqp"}    - SQP with BFGS (Kraft 1988, \pkg{nloptr}), no fallback
#'     \item \code{"auglag"}   - Augmented Lagrangian (nloptr), random starts
#'   }
#' @param lsf_list  Reserved for a future system-FORM (not active).
#' @param debug.level Verbosity level: 0 = silent, 1 = basic info, 2 = detailed output.
#'
#' @return A list (content depends on optim_type):
#'   * \code{beta}       - Hasofer-Lind reliability index
#'   * \code{pf}         - probability of failure \eqn{\Phi(-\beta)}
#'   * \code{x_points}   - design point in physical space
#'   * \code{u_points}   - design point in standard normal space
#'   * \code{dy}         - gradient in standard normal space at the design point
#'   * \code{alpha}      - direction cosines
#'   * \code{n_sim}      - number of (i)HLRF iterations performed
#'   * \code{converged}  - TRUE if the convergence criteria are met
#'   * \code{solver}     - algorithm that delivered the result ("HLRF", "iHLRF", "SLSQP")
#'   * \code{fallback}   - TRUE if \code{optim_type = "auto"} had to switch from SLSQP to iHLRF,
#'     otherwise FALSE
#'   * \code{conv_info}  - named vector of final convergence metrics
#'
#' @references
#' Hasofer, A. M., & Lind, N. C. (1974). An exact and invariant first-order reliability format.
#'   \emph{Journal of the Engineering Mechanics Division, ASCE}, 100(1), 111-121.
#'
#' Rackwitz, R., & Fieszler, B. (1978). Structural reliability under combined random load sequences.
#'   \emph{Computers & Structures}, 9(5), 489-494.
#'
#' Zhang, Y., & Der Kiureghian, A. (1997). Finite element reliability methods for inelastic
#'   structures. Report UCB/SEMM-97/05, University of California, Berkeley.
#'
#' Kraft, D. (1988). A software package for sequential quadratic programming.
#'   Technical Report DFVLR-FB 88-28, Institut fuer Dynamik der Flugsysteme,
#'   Oberpfaffenhofen.
#'
#' @import nloptr
#' @author (C) 2021-2026 K. Nille-Hauf, J.P. Schulze-Ardey, T. Feiri, M. Ricker, T. Lux
#'   -- Hochschule Biberach (until 2022),
#'   TU Dortmund University - Chair of Structural Concrete (since 2023)
#' @export
#'
FORM <- function(lsf      = NULL,
                 lDistr,
                 n_optim  = 100L,
                 loctol   = 1e-5,
                 optim_type = "auto",
                 lsf_list   = NULL,
                 debug.level = 0) {

  debug.TAG <- "FORM_OP"
  debug.print(debug.level, debug.TAG, c(optim_type),
              msg = "FORM Algorithm started with Solvingtype:")
  tic <- proc.time()

  # Check if lDistr ist properly defined
  .check_lDistr(lDistr, "FORM")

  # --- Default guards ---------------------------------------------------
  if (is.null(loctol) || is.na(loctol)) loctol <- 1e-4
  if (is.null(n_optim) || is.na(n_optim)) n_optim <- 100L
  n_optim <- as.integer(n_optim)
  optim_l <- tolower(optim_type)

  n_vars <- length(lDistr)
  U_C    <- 7                      # exact Rosenblatt map for |u_i| <= 7 (p ~ 1e-12)
  U_MAX  <- 100                    # iterates are limited to |u_i| <= 100

  # =====================================================================
  # Shared helper functions
  # =====================================================================

  ## Safe CDF / PDF wrappers --------------------------------------------
  safe_p <- function(fun, x) {
    val <- fun(x)
    min(max(val, .Machine$double.eps), 1 - .Machine$double.eps)
  }
  safe_d <- function(fun, x) {
    val <- fun(x)
    if (is.na(val) || val == 0) val <- .Machine$double.eps
    val
  }

  ## Central finite difference (optimal step ~ eps^(1/3)) ----------------
  finite_diff <- function(fun_obj, x) {
    dx <- numeric(length(x))
    for (i in seq_along(x)) {
      h  <- .Machine$double.eps^(1 / 3) * max(abs(x[i]), 1)
      xp <- x; xm <- x
      xp[i] <- x[i] + h
      xm[i] <- x[i] - h
      fp <- fun_obj(xp)
      fm <- fun_obj(xm)
      dx[i] <- if (!is.finite(fp) || !is.finite(fm)) NA_real_ else (fp - fm) / (xp[i] - xm[i])
    }
    dx
  }

  ## Rosenblatt transformation: physical x -> standard y ---------------
  x_to_y <- function(x) {
    y    <- numeric(n_vars)
    sd_s <- numeric(n_vars)
    m_s  <- numeric(n_vars)
    for (i in seq_len(n_vars)) {
      ld      <- lDistr[[i]][[1]]
      p_val   <- safe_p(ld$p, x[i])
      d_val   <- safe_d(ld$d, x[i])
      sd_s[i] <- stats::dnorm(stats::qnorm(p_val)) / d_val
      m_s[i]  <- x[i] - sd_s[i] * stats::qnorm(p_val)
      y[i]    <- (x[i] - m_s[i]) / sd_s[i]
    }
    list(y = y, sd_s = sd_s, m_s = m_s)
  }

  ## Standard u -> physical x and slope dx/du --------------------------
  ## For |u_i| <= U_C the exact map x = F^-1(Phi(u)) is used (slope dnorm(u)/f(x)).
  ## pnorm(u) is 1 in double precision for u > 8.3, so beyond U_C the map is
  ## continued linearly with the slope at U_C. This is exact for normal
  ## variables and the usual local normal approximation otherwise
  ## (relevant only for pf < 1e-12).
  map_u <- function(u) {
    x     <- numeric(n_vars)
    slope <- numeric(n_vars)
    for (i in seq_len(n_vars)) {
      ld <- lDistr[[i]][[1]]
      uc <- min(max(u[i], -U_C), U_C)
      xc <- ld$q(min(max(stats::pnorm(uc), 1e-300), 1 - 2^-53))
      sl <- stats::dnorm(uc) / safe_d(ld$d, xc)
      x[i]     <- xc + (u[i] - uc) * sl
      slope[i] <- sl
    }
    list(x = x, slope = slope)
  }
  y_to_x <- function(y) map_u(y)$x

  ## Limit state and gradient in standard normal space -------------------
  ## (d g / d u_i = d g / d x_i * sd_s_i is exact for the Rosenblatt map)
  G_u <- function(u) {
    g <- tryCatch(lsf(y_to_x(u)), error = function(e) NA_real_)
    if (length(g) != 1L || !is.finite(g)) NA_real_ else as.numeric(g)
  }
  grad_u <- function(u) {
    m <- map_u(u)
    finite_diff(lsf, m$x) * m$slope
  }

  ## Start point: mean (or median for heavy-tailed distributions) ------
  start_x <- function() {
    x <- numeric(n_vars)
    for (i in seq_len(n_vars)) {
      ld   <- lDistr[[i]][[1]]
      x[i] <- if (!is.null(ld$mean) && is.finite(ld$mean)) ld$mean else ld$q(0.5)
    }
    x
  }

  ## Is u a usable start (finite g, finite and non-vanishing gradient)? --
  usable_start <- function(u) {
    g  <- G_u(u)
    if (!is.finite(g)) return(FALSE)
    dg <- tryCatch(grad_u(u), error = function(e) rep(NA_real_, n_vars))
    all(is.finite(dg)) && sqrt(sum(dg^2)) > 1e-8 * max(1, abs(g))
  }

  ## Robust start point in u-space --------------------------------------
  ## 1) mean point; 2) points of growing radius on diagonals and axes
  find_start_u <- function() {
    xs <- start_x()
    u0 <- x_to_y(xs)$y
    u0[!is.finite(u0)] <- 0
    u0 <- pmin(pmax(u0, -U_MAX), U_MAX)
    if (usable_start(u0)) return(u0)

    if (debug.level >= 1) message("[FORM] Searching robust initial point")
    dirs <- list()
    for (sgn in c(1, -1)) dirs[[length(dirs) + 1L]] <- rep(sgn / sqrt(n_vars), n_vars)
    for (i in seq_len(n_vars)) for (sgn in c(1, -1)) {
      e <- numeric(n_vars); e[i] <- sgn; dirs[[length(dirs) + 1L]] <- e
    }
    for (r in c(0.5, 1, 2, 0.1, 4)) {
      for (d in dirs) {
        u_try <- u0 + r * d
        if (usable_start(u_try)) return(u_try)
      }
    }
    stop("[FORM] Could not determine valid initial point.")
  }

  ## Convergence check ----------------------------------------------------
  ##   crit_g   : |g(u)| / g_scale            < loctol
  ##   crit_kkt : ||u - (u.alpha) alpha|| / max(1, ||u||)  < loctol
  ##   crit_beta: |beta_new - beta_old| (information only)
  check_convergence <- function(u, G, dG, g_scale, beta_old = NA_real_, tol = loctol) {
    nd <- sqrt(sum(dG^2))
    if (!is.finite(G) || !all(is.finite(dG)) || nd < 1e-14) {
      return(list(converged = FALSE, crit1 = Inf, crit2 = Inf, crit3 = Inf,
                  beta = NA_real_, alpha = rep(NA_real_, n_vars)))
    }
    alpha  <- -dG / nd
    beta   <- sum(u * alpha)
    kkt    <- sqrt(sum((u - beta * alpha)^2)) / max(1, sqrt(sum(u^2)))
    crit_g <- abs(G) / g_scale
    list(converged = (crit_g < tol) && (kkt < tol),
         crit1 = if (is.na(beta_old)) Inf else abs(beta - beta_old),
         crit2 = kkt, crit3 = crit_g, beta = beta, alpha = alpha)
  }

  ## Largest step factor lam <= 1 with |u + lam d|_inf <= U_MAX ---------
  max_step <- function(u, d) {
    lam <- 1
    for (i in seq_along(u)) {
      if (d[i] > 0) lam <- min(lam, (U_MAX - u[i]) / d[i])
      if (d[i] < 0) lam <- min(lam, (-U_MAX - u[i]) / d[i])
    }
    max(lam, 0)
  }

  ## ---------------------------------------------------------------------
  ## Stage 1/2: HLRF (use_ls = FALSE) or iHLRF (use_ls = TRUE)
  ## Returns list(u, G, dG, conv, n_it, status)
  ## ---------------------------------------------------------------------
  hlrf_stage <- function(u, n_max, g_scale, use_ls, use_stall = FALSE) {
    G  <- G_u(u)
    dG <- grad_u(u)
    n_it <- 0L
    beta_old <- NA_real_
    status <- "max_iter"
    conv <- check_convergence(u, G, dG, g_scale, beta_old)
    best_crit <- max(conv$crit2, conv$crit3)   # stall detection
    stall <- 0L

    while (!conv$converged && n_it < n_max) {
      nd <- sqrt(sum(dG^2))
      if (!is.finite(G) || !all(is.finite(dG)) || nd < 1e-12) {
        status <- "gradient"
        break
      }
      alpha <- -dG / nd
      d     <- alpha * ((G - sum(u * dG)) / nd) - u     # HLRF direction
      lam_max <- max_step(u, d)

      if (!use_ls) {
        lam <- lam_max
      } else {
        # merit m(u) = 0.5 |u|^2 + c |g(u)|,  c > |u| / |grad g|
        c_pen <- max(2 * sqrt(sum(u^2)), 1) / nd
        merit <- function(v) {
          g <- G_u(v)
          if (!is.finite(g)) Inf else 0.5 * sum(v^2) + c_pen * abs(g)
        }
        m0    <- 0.5 * sum(u^2) + c_pen * abs(G)
        slope <- sum((u + c_pen * sign(G) * dG) * d)    # directional derivative (< 0)
        lam <- lam_max
        ok  <- FALSE
        for (ls in seq_len(25L)) {
          mt <- merit(u + lam * d)
          if (is.finite(mt) && mt <= m0 + 1e-1 * lam * slope) { ok <- TRUE; break }
          lam <- lam * 0.5
        }
        if (!ok) { status <- "line_search"; break }
      }

      u_new  <- u + lam * d
      G_new  <- G_u(u_new)
      dG_new <- tryCatch(grad_u(u_new), error = function(e) rep(NA_real_, n_vars))
      if (!is.finite(G_new) || !all(is.finite(dG_new))) { status <- "non_finite"; break }

      beta_old <- conv$beta
      u <- u_new; G <- G_new; dG <- dG_new
      n_it <- n_it + 1L
      conv <- check_convergence(u, G, dG, g_scale, beta_old)

      debug.print(debug.level, debug.TAG,
                  c("n_it", "beta", "step", "crit_kkt", "crit_g"),
                  c(n_it, conv$beta, lam, conv$crit2, conv$crit3))

      # no relevant progress for 10 iterations (oscillation / crawling)
      crit_now <- max(conv$crit2, conv$crit3)
      if (crit_now < 0.9 * best_crit) { best_crit <- crit_now; stall <- 0L } else stall <- stall + 1L
      if (use_stall && stall >= 10L) { status <- "stalled"; break }
    }
    if (conv$converged) status <- "converged"
    list(u = u, G = G, dG = dG, conv = conv, n_it = n_it, status = status)
  }

  ## ---------------------------------------------------------------------
  ## Stage 3: SLSQP   min 0.5|u|^2  s.t.  g(u) = 0
  ## ---------------------------------------------------------------------
  slsqp_stage <- function(u, g_scale, max_eval) {
    f_obj <- function(v) sum(v^2)
    f_gr  <- function(v) 2 * v
    n_ev  <- 0L
    h_eq  <- function(v) {
      n_ev <<- n_ev + 1L
      g <- G_u(v)
      if (!is.finite(g)) 1e6 else g / g_scale
    }
    h_jac <- function(v) {
      dg <- tryCatch(grad_u(v), error = function(e) rep(NA_real_, n_vars))
      dg[!is.finite(dg)] <- 0
      matrix(dg / g_scale, nrow = 1L)
    }
    res <- tryCatch(
      nloptr::slsqp(x0 = u, fn = f_obj, gr = f_gr,
                    lower = rep(-U_MAX, n_vars), upper = rep(U_MAX, n_vars),
                    heq = h_eq, heqjac = h_jac,
                    control = list(xtol_rel = 1e-12, ftol_rel = 1e-14,
                                   maxeval = as.integer(max_eval))),
      error = function(e) NULL)
    if (is.null(res)) return(NULL)
    u_new <- as.numeric(res$par)
    G  <- G_u(u_new)
    dG <- tryCatch(grad_u(u_new), error = function(e) rep(NA_real_, n_vars))
    list(u = u_new, G = G, dG = dG,
         # SLSQP reaches its optimum only to ~1e-6: do not demand more than 1e-5
         conv = check_convergence(u_new, G, dG, g_scale, tol = max(loctol, 1e-5)),
         n_it = n_ev)   # number of g evaluations of the optimizer
  }

  # =====================================================================
  # (1)  AUGMENTED LAGRANGIAN
  #      min ||u||  s.t.  g(u) = 0, random starts, local solver COBYLA.
  #      A run counts only if its end point lies on the limit-state
  #      surface. The return code of nloptr (usually 5, maxeval) does not
  #      tell this, and some runs stop near the origin.
  # =====================================================================
  if (optim_l == "auglag") {

    beta_vec <- vector("numeric", n_optim)
    conv_vec <- logical(n_optim)
    h_vec    <- numeric(n_optim)
    par_u    <- matrix(ncol = n_vars, nrow = n_optim)
    par_x    <- matrix(ncol = n_vars, nrow = n_optim)
    dy_mat   <- matrix(ncol = n_vars, nrow = n_optim)

    h_std <- function(y) lsf(y_to_x(y))
    unorm <- function(u) sqrt(sum(u^2))

    # tolerance for the constraint g(u) = 0, relative to g at the origin
    g0    <- G_u(rep(0, n_vars))
    g_tol <- max(10 * loctol, 1e-6) * max(1, if (is.finite(g0)) abs(g0) else 1)

    for (i in seq_len(n_optim)) {
      q_rand <- stats::qnorm(stats::runif(n_vars))
      info.print(debug.TAG, debug.level, c("optim_run_i", "/"), c(i, n_optim))
      res <- nloptr::auglag(x0 = q_rand, fn = unorm, heq = h_std, localtol = loctol)
      h_end       <- h_std(res$par)
      h_vec[i]    <- if (is.finite(h_end)) abs(h_end) else Inf
      conv_vec[i] <- h_vec[i] <= g_tol
      beta_vec[i] <- res$value
      par_u[i, ]  <- res$par
      par_x[i, ]  <- y_to_x(res$par)
      dy_mat[i, ] <- numDeriv::grad(h_std, res$par)
    }

    if (any(conv_vec)) {
      cand <- which(conv_vec)
      idx  <- cand[which.min(beta_vec[cand])]
    } else {
      warning("[FORM/auglag] no run reached the limit-state surface, the result is the run with the smallest constraint violation.",
              call. = FALSE)
      idx <- which.min(h_vec)
    }

    # beta is negative if the origin lies in the failure domain (g(0) < 0)
    beta <- if (is.finite(g0) && g0 < 0) -beta_vec[idx] else beta_vec[idx]

    duration <- proc.time() - tic
    return(list(
      method     = "FORM",
      beta       = beta,
      pf         = stats::pnorm(-beta),
      u_points   = par_u[idx, ],
      x_points   = par_x[idx, ],
      dy         = dy_mat[idx, ],
      converged  = any(conv_vec),
      optim_type = optim_type,
      runtime    = duration
    ))
  }

  # =====================================================================
  # (2)-(4)  rackfies / ihlrf / slsqp (single solver) and auto (chain)
  # =====================================================================
  if (optim_l %in% c("auto", "rackfies", "ihlrf", "slsqp")) {

    u_start  <- find_start_u()
    G_start  <- G_u(u_start)
    dG_start <- grad_u(u_start)
    g_scale  <- max(abs(G_start), sqrt(sum(dG_start^2)), 1e-12)
    n_sim    <- 0L
    fallback <- FALSE
    ev_fb    <- max(200L, 10L * n_optim)

    if (optim_l == "slsqp") {
      # --- explicit SLSQP: only SLSQP ------------------------------------
      best <- slsqp_stage(u_start, g_scale, n_optim)
      if (is.null(best)) {
        G <- G_u(u_start); dG <- dG_start
        best <- list(u = u_start, G = G, dG = dG, n_it = 0L,
                     conv = check_convergence(u_start, G, dG, g_scale))
      }
      n_sim  <- best$n_it
      solver <- "SLSQP"

    } else if (optim_l %in% c("rackfies", "ihlrf")) {
      # --- explicit (i)HLRF: only this algorithm, no stall abort ---------
      best   <- hlrf_stage(u_start, n_optim, g_scale, use_ls = (optim_l == "ihlrf"))
      n_sim  <- best$n_it
      solver <- if (optim_l == "ihlrf") "iHLRF" else "HLRF"

    } else {
      # --- auto: SLSQP first, then iHLRF (from start, then from SLSQP end point)
      sl <- slsqp_stage(u_start, g_scale, ev_fb)
      if (!is.null(sl) && isTRUE(sl$conv$converged)) {
        best   <- sl
        solver <- "SLSQP"
      } else {
        fallback <- TRUE
        if (debug.level >= 1) message("[FORM/auto] SLSQP did not converge - switching to iHLRF")
        ih <- hlrf_stage(u_start, n_optim, g_scale, use_ls = TRUE, use_stall = TRUE)
        best   <- ih
        solver <- "iHLRF"
        if (!ih$conv$converged && !is.null(sl) && all(is.finite(sl$u))) {
          ih2 <- hlrf_stage(sl$u, n_optim, g_scale, use_ls = TRUE, use_stall = TRUE)
          if (isTRUE(ih2$conv$converged)) {
            best <- ih2
          }
        }
        if (!best$conv$converged) {
          # no solver converged: return the point with the smallest g residual
          cands <- list(iHLRF = ih)
          if (!is.null(sl)) cands$SLSQP <- sl
          crit  <- vapply(cands, function(z) {
            if (is.finite(z$conv$crit3)) z$conv$crit3 else Inf
          }, numeric(1))
          best   <- cands[[which.min(crit)]]
          solver <- names(cands)[which.min(crit)]
        }
      }
      n_sim <- best$n_it
    }

    conv <- best$conv
    u    <- best$u
    if (!conv$converged) {
      if (optim_l == "auto") {
        warning(sprintf(
          "[FORM/auto] no convergence with SLSQP and iHLRF (crit_kkt = %.2e, crit_g = %.2e, tol = %.1e)",
          conv$crit2, conv$crit3, loctol), call. = FALSE)
      } else {
        warning(sprintf(
          "[FORM/%s] no convergence with this algorithm for this problem (crit_kkt = %.2e, crit_g = %.2e, tol = %.1e). Please use a different optim_type (e.g. 'auto' or 'slsqp').",
          optim_l, conv$crit2, conv$crit3, loctol), call. = FALSE)
      }
    }
    if (!is.finite(conv$beta)) {
      warning(sprintf("[FORM/%s] beta not finite.", optim_l), call. = FALSE)
    }

    x_fin <- y_to_x(u)
    beta  <- conv$beta
    duration <- proc.time() - tic

    return(list(
      method     = "FORM",
      beta       = beta,
      pf         = stats::pnorm(-beta),
      x_points   = x_fin,
      u_points   = u,
      dy         = best$dG,
      alpha      = conv$alpha,
      diff       = conv$crit2,
      n_sim      = n_sim,
      converged  = isTRUE(conv$converged),
      solver     = solver,
      fallback   = fallback,
      conv_info  = unlist(conv[c("converged", "crit1", "crit2", "crit3")]),
      optim_type = optim_type,
      runtime    = duration[1:5]
    ))
  }

  # =====================================================================
  # (4)  SYSTEM-FORM  (Rackwitz G.4.7 / G.4.8)
  #
  #  Supports multiple limit-state functions  g_1, ..., g_l.
  #  Active-set strategy: add/remove at most ONE constraint per step.
  #  Iteration formula (Rackwitz G.4.8):
  #
  #    lambda^{k+1} = N^{-1} R^{-1} (A^T u  +  N^{-1} g)
  #    u^{k+1}      = A      R^{-1} (A^T u  +  N^{-1} g)
  #
  #  where  A = matrix of normalised gradients (direction cosines),
  #         N = diagonal matrix of gradient norms,
  #         R = A^T A  (correlation matrix of active constraints).
  #
  #  For a series system  pf = P(union g_j <= 0)
  #  For a parallel system pf = P(intersect g_j <= 0)
  #  (simple bounds; full multinormal bounds require FORM-SORM system)
  # =====================================================================

  ### not yet tested!

  #
  #   else if (identical(tolower(optim_type), "system")) {
  #
  #     if (is.null(lsf_list) || length(lsf_list) == 0) {
  #       stop("[FORM/system] Please supply lsf_list = list(g1, g2, ...) for system FORM.")
  #     }
  #     n_lsf <- length(lsf_list)
  #
  #     # -- shared start point
  #     x        <- start_x()
  #     beta_vec <- rep(NA_real_, n_lsf)
  #     beta_old <- rep(Inf,      n_lsf)
  #     alpha_mat <- matrix(0, nrow = n_lsf, ncol = n_vars)
  #
  #     # Transform once at start
  #     tr   <- x_to_y(x)
  #     y    <- tr$y
  #     sd_s <- tr$sd_s
  #     m_s  <- tr$m_s
  #
  #     # Active set  J* :  indices of active constraints
  #     # Initialisation: start with constraint with largest positive g value
  #     g_vals  <- sapply(lsf_list, function(gi) gi(x))
  #     J_star  <- which(g_vals == max(g_vals[g_vals > 0]))
  #     if (length(J_star) == 0) J_star <- which.max(g_vals)
  #
  #     conv_info <- NULL
  #
  #     repeat {
  #       # -----------------------------------------------------------------
  #       # Evaluate all active constraints and their gradients
  #       # -----------------------------------------------------------------
  #       g_active  <- numeric(length(J_star))   # g_j(u^k)
  #       A_cols    <- matrix(0, nrow = n_vars, ncol = length(J_star))  # normalised gradients
  #       N_diag    <- numeric(length(J_star))   # gradient norms
  #
  #       for (jj in seq_along(J_star)) {
  #         j     <- J_star[jj]
  #         g_active[jj] <- lsf_list[[j]](x)
  #
  #         dy_x  <- finite_diff(lsf_list[[j]], x)
  #         dy_y  <- dy_x * sd_s                       # gradient in u-space
  #
  #         norm_j       <- sqrt(sum(dy_y^2))
  #         N_diag[jj]   <- if (norm_j > 1e-14) norm_j else 1e-14
  #         A_cols[, jj] <- -dy_y / N_diag[jj]         # alpha_j = -nabla g / ||nabla g||
  #       }
  #
  #       # -----------------------------------------------------------------
  #       # Rackwitz G.4.8:
  #       #   R     = A^T A
  #       #   rhs   = A^T u  +  N^{-1} g
  #       #   lambda = N^{-1} R^{-1} rhs
  #       #   u_new  = A      R^{-1} rhs
  #       # -----------------------------------------------------------------
  #       R_mat  <- t(A_cols) %*% A_cols                   # (l x l) correlation matrix
  #       rhs    <- t(A_cols) %*% y + g_active / N_diag    # (l x 1)
  #
  #       # Check conditioning of R
  #       if (rcond(R_mat) < 1e-12) {
  #         warning("[FORM/system] Correlation matrix R is near-singular ",
  #                 "(nearly parallel constraints). Attempting Moore-Penrose inverse.")
  #         R_inv <- MASS::ginv(R_mat)
  #       } else {
  #         R_inv <- solve(R_mat)
  #       }
  #
  #       lambda_new <- (1 / N_diag) * (R_inv %*% rhs)    # (l x 1)
  #       y_new      <- A_cols %*% (R_inv %*% rhs)        # (n x 1)
  #
  #       # Back-transform y_new -> x_new
  #       x_neu <- y_to_x(as.vector(y_new))
  #       delta_x <- sqrt(sum((x_neu - x)^2))
  #
  #       # -----------------------------------------------------------------
  #       # Active-set update (Rackwitz: add/remove at most ONE per step)
  #       # -----------------------------------------------------------------
  #       g_all    <- sapply(lsf_list, function(gi) gi(x_neu))
  #
  #       # (a) Add constraint with largest g > 0 not yet in J*
  #       not_active <- setdiff(seq_len(n_lsf), J_star)
  #       g_not      <- g_all[not_active]
  #       if (any(g_not > 0)) {
  #         add_idx <- not_active[which.max(g_not)]
  #         J_star  <- sort(c(J_star, add_idx))
  #         if (debug.level >= 1)
  #           message(sprintf("[FORM/system] Added constraint %d to active set", add_idx))
  #       }
  #
  #       # (b) Remove constraint with most negative lambda (not feasible as active)
  #       lambda_vec <- as.vector(lambda_new)
  #       if (any(lambda_vec < 0)) {
  #         rm_jj  <- which.min(lambda_vec)
  #         rm_idx <- J_star[rm_jj]
  #         J_star <- J_star[-rm_jj]
  #         if (debug.level >= 1)
  #           message(sprintf("[FORM/system] Removed constraint %d from active set", rm_idx))
  #       }
  #
  #       # Safety: J* must not be empty
  #       if (length(J_star) == 0) {
  #         J_star <- which.max(g_all)
  #       }
  #
  #       # -----------------------------------------------------------------
  #       # Per-constraint beta values  (distance to each hyperplane)
  #       # -----------------------------------------------------------------
  #       for (j in seq_len(n_lsf)) {
  #         dy_x_j  <- finite_diff(lsf_list[[j]], x_neu)
  #         tr_neu  <- x_to_y(x_neu)
  #         dy_y_j  <- dy_x_j * tr_neu$sd_s
  #         gn      <- sqrt(sum(dy_y_j^2))
  #         if (gn > 1e-14) {
  #           beta_old[j] <- beta_vec[j]
  #           g_j          <- lsf_list[[j]](x_neu)
  #           beta_vec[j]  <- (g_j - sum(tr_neu$y * dy_y_j)) / gn
  #           alpha_mat[j, ] <- -dy_y_j / gn
  #         }
  #       }
  #
  #       # -----------------------------------------------------------------
  #       # 3-criteria convergence (check over ALL active constraints)
  #       # -----------------------------------------------------------------
  #       tr_neu   <- x_to_y(x_neu)
  #       g_active_new <- sapply(J_star, function(j) lsf_list[[j]](x_neu))
  #
  #       crit1 <- max(abs(beta_vec[J_star] - beta_old[J_star]))
  #       crit3 <- max(abs(g_active_new))
  #
  #       # KT: u + A lambda  (using current active set)
  #       A_cur <- A_cols
  #       kt    <- max(abs(as.vector(tr_neu$y) + as.vector(A_cur %*% lambda_new)))
  #       crit2 <- kt
  #
  #       eps1 <- loctol
  #       eps2 <- loctol * 0.1
  #       eps3 <- loctol * 0.01
  #       converged <- (crit1 < eps1) && (crit2 < eps2) && (crit3 < eps3)
  #       conv_info <- list(converged = converged,
  #                         crit1 = crit1, crit2 = crit2, crit3 = crit3)
  #
  #       x   <- x_neu
  #       y   <- tr_neu$y
  #       sd_s <- tr_neu$sd_s
  #       m_s  <- tr_neu$m_s
  #       n_sim <- n_sim + 1L
  #
  #       debug.print(debug.level, debug.TAG,
  #                   c("n_sim", "J_star", "crit1", "crit2", "crit3"),
  #                   c(n_sim, paste(J_star, collapse = ","),
  #                     crit1, crit2, crit3))
  #
  #       if (converged || n_sim >= n_optim) break
  #     }
  #
  #     duration <- proc.time() - tic
  #
  #     # --- System failure probability (series / parallel) ----------------
  #     # Series:   Pf = 1 - Phi_n(beta; R)  -> upper/lower bounds
  #     # Simple first-order bounds (Ditlevsen):
  #     pf_individual <- stats::pnorm(-beta_vec)
  #     pf_series_ub  <- min(1, sum(pf_individual))           # union upper bound
  #     pf_series_lb  <- max(pf_individual)                   # union lower bound
  #     pf_parallel   <- prod(stats::pnorm(-beta_vec))        # intersection (independent)
  #
  #     return(list(
  #       method          = "FORM",
  #       beta            = beta_vec,             # one per limit-state
  #       beta_min        = min(beta_vec),        # governing (series) index
  #       pf_individual   = pf_individual,
  #       pf_series_lb    = pf_series_lb,
  #       pf_series_ub    = pf_series_ub,
  #       pf_parallel     = pf_parallel,
  #       x_points        = x,
  #       alpha           = alpha_mat,
  #       J_star          = J_star,
  #       n_sim           = n_sim,
  #       conv_info       = unlist(conv_info),
  #       optim_type      = optim_type,
  #       runtime         = duration[1:5]
  #     ))
  #   }

  # =====================================================================
  # Unknown optim_type
  # =====================================================================
  stop(paste0(
    "[FORM] Unknown optim_type: '", optim_type,
    "'. Valid options: 'auto', 'ihlrf', 'rackfies', 'slsqp', 'auglag'."
  ))
}
