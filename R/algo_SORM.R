#' @name SORM
#' @title Second-Order Reliability Method (SORM)
#' @description
#' Second-Order Reliability Method (SORM) for the approximation
#' of failure probabilities in structural reliability analysis.
#'
#' The SORM extends the First-Order Reliability Method (FORM)
#' by incorporating second-order information of the limit-state function
#' at the design point. By accounting for local curvature effects,
#' SORM generally provides more accurate failure probability estimates
#' for nonlinear limit-state functions.
#'
#' Failure is defined as \eqn{g(x) < 0} (see \code{\link{FORM}}). The
#' curvature correction is applied to the FORM result.
#'
#' The method evaluates the principal curvatures of the limit-state surface
#' in standard normal space and applies an asymptotic correction
#' to the FORM probability of failure.
#'
#' The curvatures are determined numerically from second derivatives of the
#' limit-state function in standard normal space. For distributions that are
#' not smooth, e.g. empirical distributions (\code{"emp"}), finite differences
#' are dominated by sampling noise. If the second derivatives depend strongly on
#' the step size, a local quadratic regression around the design point is used
#' instead. The curvatures, and thus the SORM result, are then less precise than
#' for smooth distributions and should be checked against FORM or a simulation
#' method.
#'
#' @param lsf objective function with limit state function in form of
#'   \code{function(x) \{ x[1] + x[2] + ... \}}.
#' @param lDistr list of distributions (distribution objects of TesiproV)
#' @param debug.level If 0 no additional info if 2 high output during calculation
#'
#' @return SORM returns an object containing the following elements:
#' \itemize{
#'   \item \code{beta}: Hasofer-Lind reliability index.
#'   \item \code{pf}: Estimated probability of failure.
#'   \item \code{u_points}: Design point(s) in the standard normal space (\eqn{u}-space).
#'   \item \code{dy}: Gradient(s) of the limit-state function at the design point.
#' }
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#' @references
#' Breitung, K. (1989).
#' Asymptotic approximations for probability integrals.
#' \emph{Probabilistic Engineering Mechanics}, 4(4), 187-190.
#'
#' Cai, G. Q., & Elishakoff, I. (1994).
#' Refined second-order reliability analysis.
#' \emph{Structural Safety}, 14(4), 267-276.
#'
#' Hohenbichler, M., Gollwitzer, S., Kruse, W., & Rackwitz, R. (1987).
#' New light on first- and second-order reliability methods.
#' \emph{Structural Safety}, 4, 267-284.
#'
#' Tvedt, L. (1990).
#' Distribution of quadratic forms in normal space -
#' Applications to structural reliability.
#' \emph{Journal of Engineering Mechanics}, 116(6), 1183-1197.
#'
#' Marelli, S., & Sudret, B. (2014).
#' UQLab: A framework for uncertainty quantification in Matlab.
#' In \emph{Proceedings of the 2nd International Conference on
#' Vulnerability, Risk Analysis and Management (ICVRAM2014)},
#' Liverpool, United Kingdom, 2554-2563.
#'
#' Lacaze, S., & Missoum, S. (2015).
#' CODES: A toolbox for computational design (Version 1.0).
#' Retrieved from \url{https://sites.arizona.edu/codes/}
#'
#' Wu, X. Z. (2017).
#' Implementing statistical fitting and reliability analysis
#' for geotechnical engineering problems in R.
#' \emph{Georisk: Assessment and Management of Risk for Engineered Systems and Geohazards},
#' 11(2), 173-188.
#' @export

SORM <- function(lsf,
                 lDistr,
                 debug.level = 0) {
  debug.TAG <- "SORM_RABU"
  debug.print(debug.level, debug.TAG, c(TRUE), msg = "SORM started...")
  tic <- proc.time()

  # ceck if lDistr is properly defined
  .check_lDistr(lDistr, "SORM")
  
  # u-space -> x-space, then evaluation of the limit-state function
  Tlsf <- function(u) {
    x <- numeric(length(u))
    for (i in seq_along(u)) {
      p <- min(max(stats::pnorm(u[i]), 1e-15), 1 - 1e-15)
      x[i] <- lDistr[[i]][[1]]$q(p)
    }
    lsf(x)
  }

  # Orthonormal rotation matrix, last row = alphas
  GRAM_RABU <- function(alphas) {
    dimension <- length(alphas)
    Base <- diag(dimension)
    Base[, dimension] <- alphas
    NewBase <- Base
    Rot <- matrix(rep(0, dimension * dimension), nrow = dimension)
    Rot[dimension, dimension] <- 1
    for (k in dimension:1) {
      Z <- matrix(rep(0, dimension), nrow = dimension)
      j <- k + 1
      while (j <= dimension) {
        Z <- Z + pracma::dot(Base[, k], NewBase[, j]) * NewBase[, j]
        j <- j + 1
      }
      NewBase[, k] <- Base[, k] - Z
      Rot[k, k] <- pracma::Norm(NewBase[, k], 2)
      NewBase[, k] <- NewBase[, k] / Rot[k, k]
      if (k < dimension) {
        Rot[k, (k + 1):dimension] <- t(Base[, k]) %*% NewBase[, (k + 1):dimension]
      }
    }
    return(t(NewBase))
  }

  # --- FORM: an accurate design point is required for the curvature --------
  res <- FORM(lsf, lDistr, loctol = 1e-6, n_optim = 50, optim_type = "iHLRF")

  n <- length(res$u_points)
  if (n < 2L) stop("SORM requires at least two random variables.", call. = FALSE)

  # --- Curvatures in u-space at the design point ----------------------------
  gradient <- res$dy # u-space gradient at the design point
  alphas <- gradient / pracma::Norm(gradient)
  R <- GRAM_RABU(alphas)
  H <- pracma::hessian(Tlsf, res$u_points)
  # Finite differences are unreliable for non-smooth transformations (e.g.
  # empirical distributions). If the estimate depends strongly on the step
  # size, a local quadratic regression is used instead.
  H_wide <- pracma::hessian(Tlsf, res$u_points, h = 0.1)
  if (anyNA(H) || anyNA(H_wide) ||
      max(abs(H - H_wide)) > 0.05 * max(abs(H_wide), 1e-8)) {
    H <- .hessian_regression(Tlsf, res$u_points, r = 1)
  }
  A <- R %*% H %*% t(R) / pracma::Norm(gradient)
  A_red <- A[seq_len(n - 1L), seq_len(n - 1L), drop = FALSE]

  if (anyNA(A_red) || any(!is.finite(A_red))) {
    warning("No curvature found. Maybe bad design point from FORM.", call. = FALSE)
    curvatures <- rep(0, n - 1L)
  } else {
    curvatures <- eigen((A_red + t(A_red)) / 2,
                        symmetric = TRUE, only.values = TRUE)$values
  }
  if (any(1 + res$beta * curvatures <= 0)) {
    warning("Breitung/Hohenbichler not applicable: 1 + beta * kappa <= 0.", call. = FALSE)
  }
  debug.print(debug.level, debug.TAG, c("curvatures"), curvatures)

  # --- Breitung ------------------------------------------------------------
  SORMpfBreitung <- res$pf / sqrt(prod(1 + res$beta * curvatures))
  SORMbetaBreitung <- -stats::qnorm(SORMpfBreitung)

  # --- Hohenbichler --------------------------------------------------------
  CDFminusBeta <- stats::pnorm(-res$beta, 0, 1)
  RatioCDFPDF <- stats::dnorm(-res$beta, 0, 1) / CDFminusBeta
  HohenbichlerVec <- 1 / sqrt(1 + RatioCDFPDF * curvatures)
  SORMpfHB <- CDFminusBeta * prod(HohenbichlerVec)
  SORMbetaHB <- (-1) * stats::qnorm(SORMpfHB)

  # --- Tvedt ---------------------------------------------------------------
  A11 <- 1 / sqrt(prod(1 + res$beta * curvatures))
  A12 <- res$beta * res$pf - stats::dnorm(-res$beta)
  A1 <- res$pf * A11
  A2 <- A12 * (A11 - 1 / sqrt(prod(1 + (res$beta + 1) * curvatures)))
  A3 <- (res$beta + 1) * A12 *
    (A11 - Re(1 / sqrt(prod(1 + (res$beta + 1i) * curvatures))))
  SORMpfTvedt <- A1 + A2 + A3
  SORMbetaTvedt <- -stats::qnorm(SORMpfTvedt)

  duration <- proc.time() - tic

  output <- list(
    "method" = "SORM",
    "beta" = SORMbetaBreitung,
    "pf" = SORMpfBreitung,
    "FORM_pf" = res$pf,
    "FORMbeta" = res$beta,
    "curvatures" = curvatures,
    "design.point_X" = res$x_points,
    "SORMpfTvedt" = SORMpfTvedt,
    "SORMbetaTvedt" = SORMbetaTvedt,
    "SORMpfHB" = SORMpfHB,
    "SORMbetaHB" = SORMbetaHB,
    "runtime" = duration[1:5]
  )

  debug.print(debug.level, debug.TAG, c(duration), msg = "SORM finished in [s]: ")
  return(output)
}

# Hessian from a quadratic least-squares fit in a ball of radius r around x0.
# The symmetric design makes the fit robust against noise and kinks.
.hessian_regression <- function(f, x0, r = 1, m = 20) {
  n <- length(x0)
  idx <- which(upper.tri(diag(n), diag = TRUE), arr.ind = TRUE)
  N <- m * (n + 1) * (n + 2) / 2
  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = globalenv()))
  set.seed(1)
  D <- matrix(stats::rnorm(N * n), N, n)
  D <- D / sqrt(rowSums(D^2)) * r * stats::runif(N)^(1 / n)
  D <- rbind(D, -D)
  y <- apply(D, 1, function(d) f(x0 + d))
  X <- cbind(1, D, vapply(seq_len(nrow(idx)),
                          function(k) D[, idx[k, 1]] * D[, idx[k, 2]],
                          numeric(nrow(D))))
  cf <- tryCatch(qr.solve(X, y), error = function(e) rep(NA_real_, ncol(X)))
  cq <- cf[-seq_len(n + 1)]
  H <- matrix(0, n, n)
  for (k in seq_len(nrow(idx))) {
    i <- idx[k, 1]; j <- idx[k, 2]
    if (i == j) H[i, i] <- 2 * cq[k] else H[i, j] <- H[j, i] <- cq[k]
  }
  H
}
