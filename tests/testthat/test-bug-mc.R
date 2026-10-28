# Verifikation der Befunde zu MC_CRUDE, MC_IS_single und MC_IS_system
# Ausfuehren: devtools::test(filter = "bug-mc")
#
# Relative Abweichungen werden mit expect_lt(abs(pf / pf_ref - 1), tol) geprueft,
# weil expect_equal(tolerance = ) bei kleinen Sollwerten (< tolerance) absolut
# vergleicht und dann wirkungslos ist.

# ---------------------------------------------------------------- MC_CRUDE
test_that("CRUDE-0 (Kontrolle): pf stimmt bei vorhandenen Versagensereignissen", {
  R <- PROB_BASEVAR$new(Name = "R", DistributionType = "norm", Mean = 10, Sd = 1)
  S <- PROB_BASEVAR$new(Name = "S", DistributionType = "norm", Mean = 5, Sd = 2)
  lDistr <- list(R$getlDistr(), S$getlDistr())

  res <- suppressWarnings(suppressMessages(MC_CRUDE(
    function(x) x[1] - x[2], lDistr,
    cov_user = 0.02, n_batch = 10000, n_max = 5e5, use_threads = 1,
    backend = "parallel", dataRecord = FALSE, debug.level = 0, seed = 1
  )))

  expect_lt(abs(res$pf / pnorm(-5 / sqrt(5)) - 1), 0.08)
})

test_that("CRUDE-1: ohne Versagensereignis ist pf nicht 1 und beta nicht -Inf", {
  R <- PROB_BASEVAR$new(Name = "R", DistributionType = "norm", Mean = 10, Sd = 1)
  S <- PROB_BASEVAR$new(Name = "S", DistributionType = "norm", Mean = 5, Sd = 2)
  lDistr <- list(R$getlDistr(), S$getlDistr())

  res <- suppressWarnings(suppressMessages(MC_CRUDE(
    function(x) 100 + x[1] - x[2], lDistr,
    cov_user = 0.05, n_batch = 100, n_max = 1000, use_threads = 1,
    backend = "parallel", dataRecord = FALSE, debug.level = 0, seed = 1
  )))

  expect_true(is.na(res$pf) || res$pf < 0.5)
  expect_false(identical(res$beta, -Inf))
})

# ------------------------------------------------------------ MC_IS_single
test_that("ISSINGLE-0 (Kontrolle): pf stimmt", {
  A <- PROB_BASEVAR$new(Name = "A", DistributionType = "norm", Mean = 0, Sd = 1)
  lDistr <- list(A$getlDistr()[[1]])

  res <- suppressWarnings(suppressMessages(MC_IS_single(
    lsf = function(x) 3 - x[1], lDistr = lDistr,
    cov_user = 0.02, n_batch = 10000, n_max = 1e6, use_threads = 1,
    backend = "parallel", dataRecord = FALSE, densityType = "norm", dps = NULL,
    debug.level = 0, streams = NULL, libPaths_local = .libPaths(),
    seed = 1, stability_mode = "robust"
  )))

  expect_lt(abs(res$pf / pnorm(-3) - 1), 0.08)
})

test_that("ISSINGLE-1: ESS wird berechnet (nicht konstant 0)", {
  A <- PROB_BASEVAR$new(Name = "A", DistributionType = "norm", Mean = 0, Sd = 1)
  lDistr <- list(A$getlDistr()[[1]])

  res <- suppressWarnings(suppressMessages(MC_IS_single(
    lsf = function(x) 3 - x[1], lDistr = lDistr,
    cov_user = 0.02, n_batch = 10000, n_max = 1e6, use_threads = 1,
    backend = "parallel", dataRecord = FALSE, densityType = "norm", dps = NULL,
    debug.level = 0, streams = NULL, libPaths_local = .libPaths(),
    seed = 1, stability_mode = "robust"
  )))

  expect_gt(res$ESS, 0)
})

test_that("ISSINGLE-2: vorgegebener Designpunkt (dps) ist nutzbar", {
  A <- PROB_BASEVAR$new(Name = "A", DistributionType = "norm", Mean = 0, Sd = 1)
  lDistr <- list(A$getlDistr()[[1]])

  expect_no_error(res <- suppressWarnings(suppressMessages(MC_IS_single(
    lsf = function(x) 3 - x[1], lDistr = lDistr,
    cov_user = 0.02, n_batch = 10000, n_max = 1e6, use_threads = 1,
    backend = "parallel", dataRecord = FALSE, densityType = "norm", dps = 3,
    debug.level = 0, streams = NULL, libPaths_local = .libPaths(),
    seed = 1, stability_mode = "robust"
  ))))

  expect_lt(abs(res$pf / pnorm(-3) - 1), 0.08)
})

# ------------------------------------------------------------ MC_IS_system
test_that("ISSYS-1: serielles System, pf stimmt und ESS_modes sind positiv", {
  A <- PROB_BASEVAR$new(Name = "A", DistributionType = "norm", Mean = 0, Sd = 1)
  B <- PROB_BASEVAR$new(Name = "B", DistributionType = "norm", Mean = 0, Sd = 1)
  lsfs <- list(function(x) 3 - x[1], function(x) 3.2 - x[1])
  lDistr <- list(list(A$getlDistr()[[1]]), list(B$getlDistr()[[1]]))
  p_exact <- 1 - (1 - pnorm(-3)) * (1 - pnorm(-3.2))

  res <- suppressWarnings(suppressMessages(MC_IS_system(
    lsf = lsfs, lDistr = lDistr, cov_user = 0.015, n_batch = 10000, n_max = 2e6,
    use_threads = 1, sys_type = "serial", dataRecord = FALSE, beta_l = 8,
    densityType = "norm", debug.level = 0, seed = 1, backend = "parallel",
    stability_mode = "fast"
  )))

  expect_lt(abs(res$pf / p_exact - 1), 0.06)
  expect_true(all(res$ESS_modes > 0))
})

test_that("ISSYS-2: paralleles System, pf stimmt (Gewichte passen zur Proposal)", {
  A <- PROB_BASEVAR$new(Name = "A", DistributionType = "norm", Mean = 0, Sd = 1)
  B <- PROB_BASEVAR$new(Name = "B", DistributionType = "norm", Mean = 0, Sd = 1)
  lsfs <- list(function(x) 1 - x[1], function(x) 1.2 - x[1])
  lDistr <- list(list(A$getlDistr()[[1]]), list(B$getlDistr()[[1]]))

  res <- suppressWarnings(suppressMessages(MC_IS_system(
    lsf = lsfs, lDistr = lDistr, cov_user = 0.015, n_batch = 10000, n_max = 2e6,
    use_threads = 1, sys_type = "parallel", dataRecord = FALSE, beta_l = 8,
    densityType = "norm", debug.level = 0, seed = 1, backend = "parallel",
    stability_mode = "robust"
  )))

  expect_lt(abs(res$pf / (pnorm(-1) * pnorm(-1.2)) - 1), 0.06)
})

test_that("ISSYS-3: gemeldeter CoV passt zur empirischen Streuung", {
  skip_on_cran()
  A <- PROB_BASEVAR$new(Name = "A", DistributionType = "norm", Mean = 0, Sd = 1)
  B <- PROB_BASEVAR$new(Name = "B", DistributionType = "norm", Mean = 0, Sd = 1)
  lsfs <- list(function(x) 3 - x[1], function(x) 3.2 - x[1])
  lDistr <- list(list(A$getlDistr()[[1]]), list(B$getlDistr()[[1]]))

  runs <- vapply(1:15, function(i) {
    r <- suppressWarnings(suppressMessages(MC_IS_system(
      lsf = lsfs, lDistr = lDistr, cov_user = 0.1, n_batch = 10000, n_max = 2e5,
      use_threads = 1, sys_type = "serial", dataRecord = FALSE, beta_l = 8,
      densityType = "norm", debug.level = 0, seed = i, backend = "parallel",
      stability_mode = "fast"
    )))
    c(r$pf, r$cov_mc)
  }, numeric(2))

  empirical_cv <- sd(runs[1, ]) / mean(runs[1, ])
  reported_cv <- mean(runs[2, ])
  expect_lt(empirical_cv, 1.5 * reported_cv)
})

test_that("ISSYS-4: parallel system with distant design points, shared variables", {
  Z1 <- PROB_BASEVAR$new(Name = "Z1", DistributionType = "norm", Mean = 0, Sd = 1)
  Z2 <- PROB_BASEVAR$new(Name = "Z2", DistributionType = "norm", Mean = 0, Sd = 1)
  lsfs <- list(function(x) 3 - x[1], function(x) 3 - x[2])
  lDistr <- list(
    list(Z1$getlDistr()[[1]], Z2$getlDistr()[[1]]),
    list(Z1$getlDistr()[[1]], Z2$getlDistr()[[1]])
  )
  p_exact <- pnorm(-3)^2

  res <- suppressWarnings(suppressMessages(MC_IS_system(
    lsf = lsfs, lDistr = lDistr, cov_user = 0.03, n_batch = 10000, n_max = 2e5,
    use_threads = 1, sys_type = "parallel", dataRecord = FALSE, beta_l = 8,
    densityType = "norm", debug.level = 0, seed = 1, backend = "parallel",
    stability_mode = "robust"
  )))

  expect_gt(res$pf, 0)
  expect_lt(abs(res$pf / p_exact - 1), 0.1)
})

test_that("ISSYS-5: serial system with correlated failure domains, shared variables", {
  Z1 <- PROB_BASEVAR$new(Name = "Z1", DistributionType = "norm", Mean = 0, Sd = 1)
  Z2 <- PROB_BASEVAR$new(Name = "Z2", DistributionType = "norm", Mean = 0, Sd = 1)
  c12 <- 0.9
  s12 <- sqrt(1 - c12^2)
  lsfs <- list(function(x) 3 - x[1], function(x) 3 - (c12 * x[1] + s12 * x[2]))
  lDistr <- list(
    list(Z1$getlDistr()[[1]], Z2$getlDistr()[[1]]),
    list(Z1$getlDistr()[[1]], Z2$getlDistr()[[1]])
  )
  # P(Z1 > 3 or W > 3) = 2 * P(Z > 3) - P(Z1 > 3 and W > 3), W = c12 * Z1 + s12 * Z2
  p_both <- integrate(
    function(z) dnorm(z) * pnorm((c12 * z - 3) / s12),
    lower = 3, upper = Inf
  )$value
  p_exact <- 2 * pnorm(-3) - p_both

  res <- suppressWarnings(suppressMessages(MC_IS_system(
    lsf = lsfs, lDistr = lDistr, cov_user = 0.02, n_batch = 10000, n_max = 1e6,
    use_threads = 1, sys_type = "serial", dataRecord = FALSE, beta_l = 8,
    densityType = "norm", debug.level = 0, seed = 1, backend = "parallel",
    stability_mode = "robust"
  )))

  expect_lt(abs(res$pf / p_exact - 1), 0.06)
})

test_that("ISSYS-6: estimator is accurate at a moderate sample budget", {
  Z1 <- PROB_BASEVAR$new(Name = "Z1", DistributionType = "norm", Mean = 0, Sd = 1)
  Z2 <- PROB_BASEVAR$new(Name = "Z2", DistributionType = "norm", Mean = 0, Sd = 1)
  c12 <- 0.9
  s12 <- sqrt(1 - c12^2)
  lsfs <- list(function(x) 3 - x[1], function(x) 3 - (c12 * x[1] + s12 * x[2]))
  lDistr <- list(
    list(Z1$getlDistr()[[1]], Z2$getlDistr()[[1]]),
    list(Z1$getlDistr()[[1]], Z2$getlDistr()[[1]])
  )
  p_both <- integrate(
    function(z) dnorm(z) * pnorm((c12 * z - 3) / s12),
    lower = 3, upper = Inf
  )$value
  p_exact <- 2 * pnorm(-3) - p_both

  pf <- vapply(1:4, function(i) {
    r <- suppressWarnings(suppressMessages(MC_IS_system(
      lsf = lsfs, lDistr = lDistr, cov_user = 1e-6, n_batch = 10000, n_max = 6e4,
      use_threads = 1, sys_type = "serial", dataRecord = FALSE, beta_l = 8,
      densityType = "norm", debug.level = 0, seed = i, backend = "parallel",
      stability_mode = "robust"
    )))
    r$pf
  }, numeric(1))

  expect_lt(max(abs(pf / p_exact - 1)), 0.03)
})

test_that("ISSYS-7: nearly identical design points use a single shift", {
  Z1 <- PROB_BASEVAR$new(Name = "Z1", DistributionType = "norm", Mean = 0, Sd = 1)
  lsfs <- list(function(x) 3 - x[1], function(x) 3.1 - x[1])
  lDistr <- list(list(Z1$getlDistr()[[1]]), list(Z1$getlDistr()[[1]]))
  p_exact <- pnorm(-3)

  msgs <- character(0)
  out <- capture.output(
    res <- withCallingHandlers(
      suppressWarnings(MC_IS_system(
        lsf = lsfs, lDistr = lDistr, cov_user = 0.03, n_batch = 10000, n_max = 5e5,
        use_threads = 1, sys_type = "serial", dataRecord = FALSE, beta_l = 8,
        densityType = "norm", debug.level = 1, seed = 1, backend = "parallel",
        stability_mode = "robust"
      )),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
  )

  expect_true(any(grepl("multimodal=FALSE", msgs)))
  expect_lt(abs(res$pf / p_exact - 1), 0.1)
})
