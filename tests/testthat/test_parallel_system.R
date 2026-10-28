# Parallel system with a common load: reference solution from conditional independence
#
# A parallel system of n components fails only if ALL components fail. All
# components are loaded by the same load S, the resistances R_i are independent.
# Given S = s, all components fail with probability F_R(s)^n, hence
#   P_f = integral of f_S(s) * F_R(s)^n ds
# which is evaluated by numerical integration and serves as an independent
# reference for the parallel system.

test_that("PARSYS-1: parallel system with a common load matches the exact solution (MCIS)", {
  n <- 3
  mu_R <- 3
  sd_R <- 0.5
  mu_S <- 1.5
  sd_S <- 0.5

  mR <- vector("list", n)
  for (i in seq_len(n)) {
    mR[[i]] <- PROB_BASEVAR$new(
      Id = i, Name = paste0("R", i),
      DistributionType = "norm", Mean = mu_R, Sd = sd_R
    )
  }
  mS <- PROB_BASEVAR$new(
    Id = n + 1, Name = "S",
    DistributionType = "norm", Mean = mu_S, Sd = sd_S
  )

  lsfs <- vector("list", n)
  for (i in seq_len(n)) {
    lsfs[[i]] <- SYS_LSF$new(name = paste0("g", i), vars = list(mR[[i]], mS))
    lsfs[[i]]$func <- function(x) x[1] - x[2]
    lsfs[[i]]$check()
  }

  p_exact <- integrate(
    function(s) dnorm(s, mu_S, sd_S) * pnorm((s - mu_R) / sd_R)^n,
    lower = mu_S - 10 * sd_S, upper = mu_S + 10 * sd_S,
    rel.tol = 1e-10, subdivisions = 1000
  )$value
  beta_exact <- -qnorm(p_exact)

  ps <- SYS_PROB$new(sys_input = lsfs, probMachines = list(), sys_type = "parallel")
  suppressWarnings(suppressMessages(with_test_future_plan(future::sequential, {
    ps$calculateSystemProbability(
      calcType = "MCIS",
      params = list(
        cov_user = 0.03, n_batch = 10000, n_max = 5e5, use_threads = 1,
        seed = 1, backend = "future"
      )
    )
  })))
  res <- ps$res_sys[[length(ps$res_sys)]]

  expect_gt(res$pf, 0)
  expect_lt(abs(res$pf / p_exact - 1), 0.05)
  expect_lt(abs(res$beta - beta_exact), 0.01)
})

test_that("PARSYS-2: simple bounds enclose the exact reliability index of the parallel system", {
  n <- 3
  mu_R <- 3
  sd_R <- 0.5
  mu_S <- 1.5
  sd_S <- 0.5

  mR <- vector("list", n)
  for (i in seq_len(n)) {
    mR[[i]] <- PROB_BASEVAR$new(
      Id = i, Name = paste0("R", i),
      DistributionType = "norm", Mean = mu_R, Sd = sd_R
    )
  }
  mS <- PROB_BASEVAR$new(
    Id = n + 1, Name = "S",
    DistributionType = "norm", Mean = mu_S, Sd = sd_S
  )

  lsfs <- vector("list", n)
  for (i in seq_len(n)) {
    lsfs[[i]] <- SYS_LSF$new(name = paste0("g", i), vars = list(mR[[i]], mS))
    lsfs[[i]]$func <- function(x) x[1] - x[2]
    lsfs[[i]]$check()
  }

  p_exact <- integrate(
    function(s) dnorm(s, mu_S, sd_S) * pnorm((s - mu_R) / sd_R)^n,
    lower = mu_S - 10 * sd_S, upper = mu_S + 10 * sd_S,
    rel.tol = 1e-10, subdivisions = 1000
  )$value
  beta_exact <- -qnorm(p_exact)

  m_form <- PROB_MACHINE$new(name = "FORM", fCall = "FORM")
  ps <- SYS_PROB$new(sys_input = lsfs, probMachines = list(m_form), sys_type = "parallel")
  ps$runMachines()
  suppressWarnings(suppressMessages(ps$calculateSystemProbability("simpleBounds")))
  res <- ps$res_sys[[length(ps$res_sys)]]

  # lower bound: largest single beta, upper bound: independent components
  expect_lte(res$min_beta, beta_exact)
  expect_gte(res$max_beta, beta_exact)
})
