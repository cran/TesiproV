# Verifikation der Befunde zu getlDistr() und prepare()
# Ausfuehren: devtools::test(filter = "bug-distributions")

test_that("DIST-1: gamma nutzt scale (nicht rate)", {
  G <- PROB_BASEVAR$new(Name = "G", DistributionType = "gamma", Mean = 10, Sd = 2)
  ld <- G$getlDistr()[[1]]  # shape = 25, scale = 0.4
  expect_equal(ld$p(10), pgamma(10, shape = 25, scale = 0.4))
  expect_equal(ld$q(0.5), qgamma(0.5, shape = 25, scale = 0.4))
})

test_that("DIST-2: exp liefert Dichte (nicht log-Dichte) und r() laeuft", {
  E <- PROB_BASEVAR$new(Name = "E", DistributionType = "exp", Mean = 2)
  ld <- E$getlDistr()[[1]]
  expect_equal(ld$d(1), dexp(1, rate = 0.5))
  expect_no_error(x <- ld$r(5))
  expect_length(x, 5)
})

test_that("DIST-3: weibull reproduziert Mean und Sd", {
  W <- PROB_BASEVAR$new(Name = "W", DistributionType = "weibull", Mean = 10, Sd = 2)
  k <- W$DistributionParameters[1]
  lam <- W$DistributionParameters[2]
  m <- lam * gamma(1 + 1 / k)
  s <- lam * sqrt(gamma(1 + 2 / k) - gamma(1 + 1 / k)^2)
  expect_equal(m, 10, tolerance = 1e-3)
  expect_equal(s, 2, tolerance = 1e-3)
})

test_that("DIST-4: gumbel ist mit Default-Package nutzbar", {
  Q <- PROB_BASEVAR$new(Name = "Q", DistributionType = "gumbel", Mean = 10, Sd = 2)
  ld <- Q$getlDistr()[[1]]
  expect_no_error(ld$d(10))
})

test_that("DIST-5 (Kontrolle): lnorm hat den vorgegebenen Mittelwert", {
  L <- PROB_BASEVAR$new(Name = "L", DistributionType = "lnorm", Mean = 10, Sd = 2)
  ld <- L$getlDistr()[[1]]
  expect_equal(integrate(function(p) ld$q(p), 0, 1)$value, 10, tolerance = 1e-3)
})
