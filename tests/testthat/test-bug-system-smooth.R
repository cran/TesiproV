# System analysis by importance sampling with the smoothed system function
# (system_method = "smooth"). The system function is the log-sum-exp
# approximation of min (serial) or max (parallel) of the single functions.
#
# Reference: for the independent functions g_i = R_i - S_i with normal
# variables, the failure probabilities p_i = Phi(-beta_i) are exact and
#   serial:   p = 1 - prod(1 - p_i)
#   parallel: p = prod(p_i)

test_that("SMOOTH-1: smooth_min and smooth_max approximate min and max", {
  x <- c(1, 2, 3)

  expect_equal(smooth_min(x, kappa = 1000), 1, tolerance = 1e-3)
  expect_equal(smooth_max(x, kappa = 1000), 3, tolerance = 1e-3)

  # log-sum-exp bounds: smooth_min <= min and smooth_max >= max
  expect_lte(smooth_min(x, kappa = 20), min(x))
  expect_gte(smooth_max(x, kappa = 20), max(x))

  # sharper approximation for larger kappa
  expect_lt(abs(smooth_min(x, kappa = 100) - 1), abs(smooth_min(x, kappa = 5) - 1))
  expect_lt(abs(smooth_max(x, kappa = 100) - 3), abs(smooth_max(x, kappa = 5) - 3))
})

test_that("SMOOTH-2: serial system with independent functions matches the exact solution", {
  skip_on_cran()

  R1 <- PROB_BASEVAR$new(Id = 1, Name = "R1", DistributionType = "norm", Mean = 2, Sd = 0.5)
  S1 <- PROB_BASEVAR$new(Id = 2, Name = "S1", DistributionType = "norm", Mean = 1, Sd = 0.5)
  R2 <- PROB_BASEVAR$new(Id = 3, Name = "R2", DistributionType = "norm", Mean = 2.3, Sd = 0.5)
  S2 <- PROB_BASEVAR$new(Id = 4, Name = "S2", DistributionType = "norm", Mean = 1, Sd = 0.5)

  g1 <- SYS_LSF$new(name = "g1", vars = list(R1, S1))
  g1$func <- function(x) x[1] - x[2]
  g2 <- SYS_LSF$new(name = "g2", vars = list(R2, S2))
  g2$func <- function(x) x[1] - x[2]

  p_exact <- 1 - (1 - pnorm(-1 / sqrt(0.5))) * (1 - pnorm(-1.3 / sqrt(0.5)))
  beta_exact <- -qnorm(p_exact)

  ps <- SYS_PROB$new(sys_input = list(g1, g2), probMachines = list(), sys_type = "serial")
  expect_equal(ps$sys_structure, "serial")

  suppressMessages(ps$calculateSystemProbability(
    calcType = "MCIS",
    params = list(
      system_method = "smooth", cov_user = 0.03, n_batch = 5000,
      n_max = 3e5, seed = 1, backend = "future"
    )
  ))
  res <- ps$res_sys[[length(ps$res_sys)]]

  expect_lt(abs(res$beta - beta_exact), 0.05)
})

test_that("SMOOTH-3: parallel system with independent functions matches the exact solution", {
  skip_on_cran()

  R1 <- PROB_BASEVAR$new(Id = 1, Name = "R1", DistributionType = "norm", Mean = 2, Sd = 0.5)
  S1 <- PROB_BASEVAR$new(Id = 2, Name = "S1", DistributionType = "norm", Mean = 1, Sd = 0.5)
  R2 <- PROB_BASEVAR$new(Id = 3, Name = "R2", DistributionType = "norm", Mean = 2.3, Sd = 0.5)
  S2 <- PROB_BASEVAR$new(Id = 4, Name = "S2", DistributionType = "norm", Mean = 1, Sd = 0.5)

  g1 <- SYS_LSF$new(name = "g1", vars = list(R1, S1))
  g1$func <- function(x) x[1] - x[2]
  g2 <- SYS_LSF$new(name = "g2", vars = list(R2, S2))
  g2$func <- function(x) x[1] - x[2]

  p_exact <- pnorm(-1 / sqrt(0.5)) * pnorm(-1.3 / sqrt(0.5))
  beta_exact <- -qnorm(p_exact)

  ps <- SYS_PROB$new(sys_input = list(g1, g2), probMachines = list(), sys_type = "parallel")
  ps$smooth_kappa <- 50

  suppressMessages(ps$calculateSystemProbability(
    calcType = "MCIS",
    params = list(
      system_method = "smooth", cov_user = 0.03, n_batch = 5000,
      n_max = 3e5, seed = 1, backend = "future"
    )
  ))
  res <- ps$res_sys[[length(ps$res_sys)]]

  expect_lt(abs(res$beta - beta_exact), 0.05)
})

test_that("SMOOTH-4: a variable shared by two functions enters the system only once", {
  skip_on_cran()

  # g_i = R_i - S with the common load S. Given S = s the resistances are
  # independent, so the exact solution follows by numerical integration.
  R1 <- PROB_BASEVAR$new(Id = 1, Name = "R1", DistributionType = "norm", Mean = 2, Sd = 0.5)
  R2 <- PROB_BASEVAR$new(Id = 2, Name = "R2", DistributionType = "norm", Mean = 2.3, Sd = 0.5)
  S <- PROB_BASEVAR$new(Id = 3, Name = "S", DistributionType = "norm", Mean = 1, Sd = 0.5)

  g1 <- SYS_LSF$new(name = "g1", vars = list(R1, S))
  g1$func <- function(x) x[1] - x[2]
  g2 <- SYS_LSF$new(name = "g2", vars = list(R2, S))
  g2$func <- function(x) x[1] - x[2]

  p_exact <- integrate(
    function(s) dnorm(s, 1, 0.5) * pnorm((s - 2) / 0.5) * pnorm((s - 2.3) / 0.5),
    lower = -10, upper = 10, rel.tol = 1e-10
  )$value
  beta_exact <- -qnorm(p_exact)

  ps <- SYS_PROB$new(sys_input = list(g1, g2), probMachines = list(), sys_type = "parallel")
  ps$smooth_kappa <- 50

  suppressMessages(ps$calculateSystemProbability(
    calcType = "MCIS",
    params = list(
      system_method = "smooth", cov_user = 0.03, n_batch = 5000,
      n_max = 3e5, seed = 1, backend = "future"
    )
  ))
  res <- ps$res_sys[[length(ps$res_sys)]]

  expect_lt(abs(res$beta - beta_exact), 0.1)
})

test_that("SMOOTH-5: an unknown system_method is rejected", {
  R1 <- PROB_BASEVAR$new(Id = 1, Name = "R1", DistributionType = "norm", Mean = 2, Sd = 0.5)
  S1 <- PROB_BASEVAR$new(Id = 2, Name = "S1", DistributionType = "norm", Mean = 1, Sd = 0.5)
  g1 <- SYS_LSF$new(name = "g1", vars = list(R1, S1))
  g1$func <- function(x) x[1] - x[2]
  g2 <- SYS_LSF$new(name = "g2", vars = list(R1, S1))
  g2$func <- function(x) x[1] - 0.9 * x[2]

  ps <- SYS_PROB$new(sys_input = list(g1, g2), probMachines = list(), sys_type = "serial")

  expect_error(
    ps$calculateSystemProbability(
      calcType = "MCIS",
      params = list(system_method = "foo", use_threads = 1)
    ),
    "system_method must be"
  )
})
