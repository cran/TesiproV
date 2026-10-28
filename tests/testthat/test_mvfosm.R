test_that("MVFOSM uses Mean and Sd, not the native parameters (lnorm)", {
  R <- PROB_BASEVAR(Id = 1, Name = "R", DistributionType = "lnorm", Mean = 30, Sd = 3)
  S <- PROB_BASEVAR(Id = 2, Name = "S", DistributionType = "norm", Mean = 20, Sd = 2)
  lsf <- function(x) x[1] - x[2]
  lDistr <- list(R$getlDistr(), S$getlDistr())

  res <- MVFOSM(lsf, lDistr)

  expect_equal(res$beta, 10 / sqrt(13), tolerance = 1e-6)
  expect_equal(res$pf, pnorm(-10 / sqrt(13)), tolerance = 1e-6)
})

test_that("MVFOSM gives correct results in SYS_PARAM (lnorm)", {
  R <- PARAM_BASEVAR(Id = 1, Name = "R", DistributionType = "lnorm", Sd = 3,
                     ParamType = "Mean", ParamValues = c(30, 35))
  S <- PROB_BASEVAR(Id = 2, Name = "S", DistributionType = "norm", Mean = 20, Sd = 2)
  lsf <- SYS_LSF(name = "RS", vars = list(R, S), func = function(x) x[1] - x[2])
  ps <- SYS_PARAM(
    sys_input = list(lsf),
    probMachines = list(PROB_MACHINE(name = "MVFOSM", fCall = "MVFOSM")),
    sys_type = "serial"
  )

  ps$runMachines()

  expect_equal(ps$beta_params[[1]][1, 1], 10 / sqrt(13), tolerance = 1e-6)
  expect_equal(ps$beta_params[[2]][1, 1], 15 / sqrt(13), tolerance = 1e-6)
})

test_that("MVFOSM stops for variables without finite Mean and Sd (lt)", {
  fc <- PROB_BASEVAR(Id = 1, Name = "fc", DistributionType = "lt",
                     DistributionParameters = c(3.75, 0.105, 3, 10))
  S <- PROB_BASEVAR(Id = 2, Name = "S", DistributionType = "norm", Mean = 20, Sd = 2)
  lsf <- function(x) x[1] - x[2]
  lDistr <- list(fc$getlDistr(), S$getlDistr())

  expect_error(MVFOSM(lsf, lDistr), "finite Mean and Sd")
})

test_that("MVFOSM reproduces reference value with lnorm variables (ffpack example)", {
  X1 <- PROB_BASEVAR(Id = 1, Name = "X1", DistributionType = "lnorm", Mean = 1, Sd = 3)
  X2 <- PROB_BASEVAR(Id = 2, Name = "X2", DistributionType = "lnorm", Mean = 1, Sd = 4)
  lsf <- function(x) 3 * x[1] - 2 * x[2]
  lDistr <- list(X1$getlDistr(), X2$getlDistr())

  res <- MVFOSM(lsf, lDistr)

  expect_equal(res$beta, 0.0830455, tolerance = 1e-6)
})

test_that("Nowak and Collins (2000), Example 5.3", {
  P <- PROB_BASEVAR(Id = 1, Name = "P", Mean = 10, Sd = 2)
  L <- PROB_BASEVAR(Id = 2, Name = "L", Mean = 8, Sd = sqrt(10*10^(-3)))
  Z <- PROB_BASEVAR(Id = 3, Name = "Z", Mean = 100*10^(-6), Sd = sqrt(400*10^(-12)))
  Fy <- PROB_BASEVAR(Id = 4, Name = "Fy", Mean = 600*10^(3), Sd = sqrt(10*10^(9)))

  lsf1 <- SYS_LSF(
    vars = list(P, L, Z,Fy),
    name = "Nowak & Collins Exp5.3"
  )
  lsf1$func <- function(P, L, Z,Fy) Z * Fy - P*L/4
  lsf1$check()

  lsf2 <- SYS_LSF(
    vars = list(P, L, Z,Fy),
    name = "Nowak & Collins Exp5.3"
  )
  lsf2$func <- function(P, L, Z,Fy) Fy - P*L/4/Z
  lsf2$check()

  ps <- SYS_PROB(
    sys_input = list(lsf1,lsf2),
    probMachines = list(PROB_MACHINE(name = "MVFOSM", fCall = "MVFOSM"))
  )

  ps$runMachines()



  expect_equal(ps$beta_single[1], 2.48, tolerance = 1e-2)
  expect_equal(ps$beta_single[2], 3.48, tolerance = 1e-2)
})


