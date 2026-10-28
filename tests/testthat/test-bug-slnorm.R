# Eigene Implementierung der verschobenen Lognormalverteilung (ersetzt brms)
# Ausfuehren: devtools::test(filter = "bug-slnorm")

test_that("SLN-1: slnorm reproduziert Mean und Sd", {
  v <- PROB_BASEVAR$new(Name = "X", DistributionType = "slnorm", Mean = 10, Sd = 2, x0 = 2)
  ld <- v$getlDistr()[[1]]
  m <- integrate(function(p) ld$q(p), 0, 1)$value
  s <- sqrt(integrate(function(p) (ld$q(p) - m)^2, 0, 1)$value)

  expect_equal(m, 10, tolerance = 1e-3)
  expect_equal(s, 2, tolerance = 1e-3)
})

test_that("SLN-2: d, p, q und r sind konsistent und liegen oberhalb der Verschiebung", {
  v <- PROB_BASEVAR$new(Name = "X", DistributionType = "slnorm", Mean = 10, Sd = 2, x0 = 2)
  ld <- v$getlDistr()[[1]]
  mn <- v$DistributionParameters[1]
  sn <- v$DistributionParameters[2]

  expect_equal(ld$d(9), dlnorm(9 - 2, mn, sn))
  expect_equal(ld$p(9), plnorm(9 - 2, mn, sn))
  expect_equal(ld$q(0.3), 2 + qlnorm(0.3, mn, sn))
  expect_equal(ld$p(ld$q(c(0.01, 0.5, 0.99))), c(0.01, 0.5, 0.99))
  expect_equal(ld$p(2), 0)

  set.seed(1)
  x <- ld$r(1e5)
  expect_true(all(x > 2))
  expect_equal(mean(x), 10, tolerance = 0.02)
})
