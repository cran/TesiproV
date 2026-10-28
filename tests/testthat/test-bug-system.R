# Verifikation der Befunde zu SYS_PROB (Simple Bounds, MVFOSM-Maschine)
# Ausfuehren: devtools::test(filter = "bug-system")

test_that("SYS-1: Simple Bounds parallel, obere Beta-Schranke aus Produkt der pf", {
  sys <- SYS_PROB$new(sys_type = "parallel")
  sys$beta_single <- matrix(c(3, 3.5), nrow = 1, dimnames = list("FORM", c("L1", "L2")))
  sys$calculateSystemProbability("simpleBounds")
  res <- sys$res_sys[[1]]

  expect_equal(res$min_beta, 3.5)
  expect_equal(res$max_beta, -qnorm(pnorm(-3) * pnorm(-3.5)))
})

test_that("SYS-2: Simple Bounds mit mehreren Maschinen bricht nicht ab", {
  sys <- SYS_PROB$new(sys_type = "parallel")
  sys$beta_single <- matrix(c(3, 3.1, 3.5, 3.4), nrow = 2,
                            dimnames = list(c("FORM", "SORM"), c("L1", "L2")))

  expect_no_error(sys$calculateSystemProbability("simpleBounds"))
})

test_that("SYS-3 (Kontrolle): Simple Bounds seriell", {
  sys <- SYS_PROB$new(sys_type = "serial")
  sys$beta_single <- matrix(c(3, 3.5), nrow = 1, dimnames = list("FORM", c("L1", "L2")))
  sys$calculateSystemProbability("simpleBounds")
  res <- sys$res_sys[[1]]

  expect_equal(res$max_beta, 3)
  expect_equal(res$min_beta, -qnorm(1 - pnorm(3) * pnorm(3.5)))
})

test_that("SYS-4: MVFOSM-Maschine ohne Option isExpression laeuft in runMachines", {
  R <- PROB_BASEVAR$new(Name = "R", DistributionType = "norm", Mean = 10, Sd = 1)
  S <- PROB_BASEVAR$new(Name = "S", DistributionType = "norm", Mean = 5, Sd = 2)
  lsf <- SYS_LSF$new(name = "g", func = function(R, S) R - S, vars = list(R, S))
  m <- PROB_MACHINE$new(name = "MVFOSM", fCall = "MVFOSM")
  sys <- SYS_PROB$new(sys_input = list(lsf), sys_type = "serial", probMachines = list(m))

  expect_no_error(sys$runMachines())
  expect_equal(unname(sys$beta_single[1, 1]), 5 / sqrt(5), tolerance = 1e-3)
})
