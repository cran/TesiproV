# Weakest-link chain: reference solution from conditional independence
#
# A chain of n links fails if the weakest link fails. All links are loaded by
# the same load S. Given S = s, the link resistances R_i are independent, so
#   P_f = integral of f_S(s) * (1 - (1 - F_R(s))^n) ds
# which is evaluated by numerical integration and serves as an independent
# reference for the serial system.

test_that("CHAIN-1: serial system of 100 links with a common load matches the exact solution", {
  skip_on_cran()

  n <- 100
  mu_R <- 1.5
  sd_R <- 0.2
  mu_S <- 0.5
  sd_S <- 0.1

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

  ps <- SYS_PROB$new(sys_input = lsfs, probMachines = list(), sys_type = "serial")
  suppressMessages(ps$calculateSystemProbability(
    calcType = "MCIS",
    params = list(
      cov_user = 0.03, n_batch = 5000, n_max = 2e5,
      seed = 1, backend = "future"
    )
  ))
  res <- ps$res_sys[[length(ps$res_sys)]]

  p_lit <- 3.845 * 10^(-4) # Spaethe S.175
  beta_lit <- 3.364 # Spaethe S.175

  expect_lt(abs(res$pf / p_lit - 1), 0.001)
  expect_lt(abs(res$beta - beta_lit), 0.001)
})
