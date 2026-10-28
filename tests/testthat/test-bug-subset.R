# Verifikation der Vermutungen zu MC_SubSam (Subset Simulation)
# Ausfuehren: devtools::test(filter = "bug-subset")
# Hinweis: SUS-4 bis SUS-6 rechnen je 12 Laeufe und brauchen einige Minuten.

test_that("SUS-1: eine einzelne Basisvariable wird unterstuetzt", {
  A <- PROB_BASEVAR$new(Name = "A", DistributionType = "norm", Mean = 0, Sd = 1)
  lDistr <- list(A$getlDistr())

  expect_no_error(res <- suppressMessages(
    MC_SubSam(function(x) 3 - x[1], lDistr, Nsubset = 2000, seed = 1)
  ))
  expect_lt(abs(res$pf / pnorm(-3) - 1), 0.5)
})


test_that("SUS-2: Seed liefert unabhaengig von der vorher gesetzten RNG-Art dasselbe Ergebnis", {
  A <- PROB_BASEVAR$new(Name = "A", DistributionType = "norm", Mean = 0, Sd = 1)
  B <- PROB_BASEVAR$new(Name = "B", DistributionType = "norm", Mean = 0, Sd = 1)
  lDistr <- list(A$getlDistr(), B$getlDistr())
  lsf <- function(x) 3.5 - (x[1] + x[2]) / sqrt(2)

  RNGkind("L'Ecuyer-CMRG")
  res_1 <- suppressMessages(MC_SubSam(lsf, lDistr, Nsubset = 2000, seed = 1))
  RNGkind("Mersenne-Twister", "Inversion")
  res_2 <- suppressMessages(MC_SubSam(lsf, lDistr, Nsubset = 2000, seed = 1))

  expect_equal(res_1$pf, res_2$pf)
})


test_that("SUS-3: Nsubset, das kein Vielfaches von 1/p0 ist, laeuft durch", {
  A <- PROB_BASEVAR$new(Name = "A", DistributionType = "norm", Mean = 0, Sd = 1)
  B <- PROB_BASEVAR$new(Name = "B", DistributionType = "norm", Mean = 0, Sd = 1)
  lDistr <- list(A$getlDistr(), B$getlDistr())
  lsf <- function(x) 3.5 - (x[1] + x[2]) / sqrt(2)

  expect_no_error(res <- suppressMessages(MC_SubSam(lsf, lDistr, Nsubset = 1005, seed = 1)))
  expect_lt(abs(res$pf / pnorm(-3.5) - 1), 1)
})


test_that("SUS-4: glatter Grenzzustand, Mittel ueber 12 Seeds ist erwartungstreu", {
  skip_on_cran()
  A <- PROB_BASEVAR$new(Name = "A", DistributionType = "norm", Mean = 0, Sd = 1)
  B <- PROB_BASEVAR$new(Name = "B", DistributionType = "norm", Mean = 0, Sd = 1)
  lDistr <- list(A$getlDistr(), B$getlDistr())
  lsf <- function(x) 3.5 - (x[1] + x[2]) / sqrt(2)

  pf <- sapply(1:12, function(s) {
    suppressMessages(MC_SubSam(lsf, lDistr, Nsubset = 4000, p0 = 0.1, seed = s))$pf
  })

  expect_lt(abs(mean(pf) / pnorm(-3.5) - 1), 0.25)
})


test_that("SUS-5: p0 = 0.29 (1/p0 nicht ganzzahlig) ist erwartungstreu", {
  skip_on_cran()
  A <- PROB_BASEVAR$new(Name = "A", DistributionType = "norm", Mean = 0, Sd = 1)
  B <- PROB_BASEVAR$new(Name = "B", DistributionType = "norm", Mean = 0, Sd = 1)
  lDistr <- list(A$getlDistr(), B$getlDistr())
  lsf <- function(x) 3.5 - (x[1] + x[2]) / sqrt(2)

  pf <- sapply(1:12, function(s) {
    suppressMessages(MC_SubSam(lsf, lDistr, Nsubset = 4000, p0 = 0.29, seed = s))$pf
  })

  expect_lt(abs(mean(pf) / pnorm(-3.5) - 1), 0.25)
})


test_that("SUS-6: diskreter Grenzzustand mit vielen Gleichstaenden ist erwartungstreu", {
  skip_on_cran()
  A <- PROB_BASEVAR$new(Name = "A", DistributionType = "norm", Mean = 0, Sd = 1)
  B <- PROB_BASEVAR$new(Name = "B", DistributionType = "norm", Mean = 0, Sd = 1)
  lDistr <- list(A$getlDistr(), B$getlDistr())
  # g < 0  <=>  x1 > 3; the LSF values lie on a 0.5 grid
  lsf <- function(x) floor(2 * (3 - x[1])) / 2 + 0 * x[2]

  pf <- sapply(1:12, function(s) {
    suppressMessages(MC_SubSam(lsf, lDistr, Nsubset = 4000, seed = s))$pf
  })

  expect_lt(abs(mean(pf) / (1 - pnorm(3.0)) - 1), 0.25)
})


test_that("MC_SubSam agrees with independent estimators for a discrete boundary", {
  x <- PROB_BASEVAR$new(
    Name = "X",
    DistributionType = "norm",
    Mean = 0,
    Sd = 1
  )

  # MC_SubSam and MC_CRUDE use the wrapped getlDistr() format.
  lDistr_wrapped <- list(x$getlDistr())

  # MC_IS_single accepts the marginal funlist directly.
  marginal <- x$getlDistr()[[1]]

  # g(x) = 0 on a non-zero-width interval.
  # Under the strict failure convention g(x) < 0, failure is X > 3.
  lsf <- function(x) {
    if (is.matrix(x)) {
      floor(2 * (3 - x[, 1L])) / 2
    } else {
      floor(2 * (3 - x[1L])) / 2
    }
  }

  pf_exact <- 1 - stats::pnorm(3)

  # Explicitly protect the boundary semantics this test is intended to cover.
  expect_identical(lsf(2.75), 0)
  expect_identical(lsf(3), 0)
  expect_lt(lsf(3.01), 0)
  expect_equal(pf_exact, 0.0013498980316301, tolerance = 1e-13)

  res_subsam <- suppressMessages(
    MC_SubSam(
      lsf = lsf,
      lDistr = lDistr_wrapped,
      Nsubset = 20000,
      p0 = 0.1,
      MaxSubsets = 10,
      seed = 1701
    )
  )

  res_crude <- suppressMessages(
    MC_CRUDE(
      lsf = lsf,
      lDistr = lDistr_wrapped,
      cov_user = 0,
      n_batch = 10000,
      n_max = 400000,
      use_threads = 1,
      backend = "parallel",
      dataRecord = FALSE,
      seed = 1702
    )
  )

  res_is <- suppressMessages(
    MC_IS(
      lsf = lsf,
      lDistr = list(marginal),
      backend = "parallel",
      use_threads = 1,
      dps = 3,
      seed = 1703,
      cov_user = 0,
      n_batch = 20000,
      n_max = 100000,
      dataRecord = FALSE,
      debug.level = 0
    )
  )

  expect_true(is.finite(res_subsam$pf))
  expect_true(is.finite(res_crude$pf))
  expect_true(is.finite(res_is$pf))

  # Statistical regression bounds against the analytical probability.
  # Keep these fixed-seed checks broad enough for sampling variability,
  # but tight enough to catch a boundary or estimator regression.
  expect_lt(abs(res_subsam$pf / pf_exact - 1), 0.25)
  expect_lt(abs(res_crude$pf / pf_exact - 1), 0.15)
  expect_lt(abs(res_is$pf / pf_exact - 1), 0.10)

  # Independent estimates should also be mutually compatible.
  expect_lt(abs(res_subsam$pf / res_crude$pf - 1), 0.30)
  expect_lt(abs(res_subsam$pf / res_is$pf - 1), 0.30)
})
