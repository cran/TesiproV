# Long running example, start manually, e.g.
#   testthat::test_file("tests/testthat/long_test_empdist.R")
# The file name does not start with "test", so devtools::test() and
# R CMD check do not run it.
#
# Example Ex. 5.11 from Nowak & Collins (2000): reliability analysis with
# student-t, empirical (log-normal sample) and Gumbel variables.

# -----------------------------------------------------------------------
# Common setup, used by all tests of this file
# -----------------------------------------------------------------------
if (.Platform$OS.type == "windows") {
  backend <- "future"
  # use_threads <- 1L
} else {
  backend <- "parallel"
  use_threads <- 20L
}

# Fy: empirical distribution from a log-normal sample (Mean = 40, Cov = 0.1)
set.seed(1)
sig_fy <- sqrt(log(1 + 0.1^2))
mue_fy <- log(40) - sig_fy^2 / 2
data_fy <- rlnorm(100000, mue_fy, sig_fy)

X1 <- PROB_BASEVAR(
  Id = 1, Name = "Z",
  DistributionType = "st",
  DistributionParameters = c(100, 4, 10000, 10000)
)

X2 <- PROB_BASEVAR(
  Id = 2, Name = "Fy",
  DistributionType = "emp",
  Package = "EnvStats",
  DistributionParameters = data_fy
)

X3 <- PROB_BASEVAR(
  Id = 3, Name = "M",
  DistributionType = "gumbel",
  Package = "evd",
  Mean = 2000, Cov = 0.1
)

lsf1 <- SYS_LSF(
  vars = list(X1, X2, X3),
  name = "Nowak & Collins Exp5.11_1 (beta=4.03)"
)
lsf1$func <- function(M, Fy, Z) Z * Fy - M
lsf1$check()

# Reference values
beta_form <- 4.02211455568222
beta_mcis <- 3.9585857628075
beta_mcc <- 3.97628577254636
beta_sus <- 4.01123480350189


# -----------------------------------------------------------------------
# FORM
# SORM is not tested here: the empirical distribution of Fy is not smooth,
# so the curvature estimate of SORM is not reliable (NaN in a test run).
# -----------------------------------------------------------------------
test_that("Ex. 5.11 Nowak & Collins (2000): FORM", {
  form_rf <- PROB_MACHINE(
    name = "FORM Rack.-Fiess.", fCall = "FORM",
    options = list(n_optim = 20, loctol = 1e-4, optim_type = "rackfies")
  )

  ps <- SYS_PROB(sys_input = list(lsf1), probMachines = list(form_rf), debug.level = 0)

  message("FORM started")
  t0 <- proc.time()[["elapsed"]]
  ps$runMachines()
  message(sprintf(
    "FORM finished after %.1f s, beta = %.4f",
    proc.time()[["elapsed"]] - t0, ps$beta_single[1]
  ))

  expect_false(is.na(ps$beta_single[1]))

  # The reference value was determined with a different sample of Fy
  expect_equal(ps$beta_single[1], beta_form, tolerance = 5e-3)
})


# -----------------------------------------------------------------------
# MC_IS
# -----------------------------------------------------------------------
test_that("Ex. 5.11 Nowak & Collins (2000): MC-IS", {
  mcis <- PROB_MACHINE(
    name = "MC IS", fCall = "MC_IS",
    options = list(
      cov_user = 0.05, n_max = 1e7, n_batch = 10000,
      use_threads = 1L, seed = 1234, backend = backend
    )
  )

  ps <- SYS_PROB(sys_input = list(lsf1), probMachines = list(mcis), debug.level = 0)

  message("MC-IS started")
  t0 <- proc.time()[["elapsed"]]
  ps$runMachines()
  message(sprintf(
    "MC-IS finished after %.1f s, beta = %.4f",
    proc.time()[["elapsed"]] - t0, ps$beta_single[1]
  ))

  expect_false(is.na(ps$beta_single[1]))
  expect_lt(abs(ps$beta_single[1] - beta_mcis) / beta_mcis, 0.02)
})


# -----------------------------------------------------------------------
# MC_CRUDE (Pf about 3e-5: cov_user = 0.1 needs about 3e6 samples)
# -----------------------------------------------------------------------
test_that("Ex. 5.11 Nowak & Collins (2000): MC-CRUDE", {
  mcc <- PROB_MACHINE(
    name = "MCC", fCall = "MC_CRUDE",
    options = list(
      cov_user = 0.1, n_max = 4e6,
      use_threads = if (backend == "future") NULL else use_threads, seed = 161268, backend = backend
    )
  )

  ps <- SYS_PROB(sys_input = list(lsf1), probMachines = list(mcc), debug.level = 0)

  message("MC-CRUDE started")
  t0 <- proc.time()[["elapsed"]]
  ps$runMachines()
  message(sprintf(
    "MC-CRUDE finished after %.1f s, beta = %.4f",
    proc.time()[["elapsed"]] - t0, ps$beta_single[1]
  ))

  expect_false(is.na(ps$beta_single[1]))
  expect_lt(abs(ps$beta_single[1] - beta_mcc) / beta_mcc, 0.03)
})


# -----------------------------------------------------------------------
# MC_SubSam (subset simulation)
# With Nsubset = 1e5 (default) this test needed more than one hour, each call of
# the empirical distribution of Fy is slow. Nsubset = 1e4 reduces the effort.
# -----------------------------------------------------------------------
test_that("Ex. 5.11 Nowak & Collins (2000): MC-SuS", {
  mcsus <- PROB_MACHINE(
    name = "MC SuS", fCall = "MC_SubSam",
    options = list(seed = 1234, Nsubset = 1e4)
  )

  ps <- SYS_PROB(sys_input = list(lsf1), probMachines = list(mcsus), debug.level = 0)

  message("MC-SuS started")
  t0 <- proc.time()[["elapsed"]]
  ps$runMachines()
  message(sprintf(
    "MC-SuS finished after %.1f s, beta = %.4f",
    proc.time()[["elapsed"]] - t0, ps$beta_single[1]
  ))

  expect_false(is.na(ps$beta_single[1]))
  expect_lt(abs(ps$beta_single[1] - beta_sus) / beta_sus, 0.05)
})

# SORM with a non-smooth (empirical) distribution:
# finite-difference Hessians are dominated by noise, a local quadratic
# regression is used as fallback. The result must be finite and close to the
# SORM result of the underlying smooth distribution.

test_that("SORM liefert mit emp-Verteilung einen endlichen Wert nahe am glatten Fall", {
  skip_if_not_installed("EnvStats")
  skip_if_not_installed("evd")

  sig_fy <- sqrt(log(1 + 0.1^2))
  mue_fy <- log(40) - sig_fy^2 / 2

  X1 <- PROB_BASEVAR(
    Id = 1, Name = "Z", DistributionType = "st",
    DistributionParameters = c(100, 4, 10000, 10000)
  )
  X3 <- PROB_BASEVAR(
    Id = 3, Name = "M", DistributionType = "gumbel",
    Package = "evd", Mean = 2000, Cov = 0.1
  )

  set.seed(1)
  X2_emp <- PROB_BASEVAR(
    Id = 2, Name = "Fy", DistributionType = "emp",
    Package = "EnvStats",
    DistributionParameters = rlnorm(20000, mue_fy, sig_fy)
  )
  X2_ln <- PROB_BASEVAR(
    Id = 2, Name = "Fy", DistributionType = "lnorm",
    DistributionParameters = c(mue_fy, sig_fy)
  )

  lsf_emp <- SYS_LSF(vars = list(X1, X2_emp, X3), name = "emp")
  lsf_emp$func <- function(M, Fy, Z) Z * Fy - M
  lsf_ln <- SYS_LSF(vars = list(X1, X2_ln, X3), name = "ln")
  lsf_ln$func <- function(M, Fy, Z) Z * Fy - M

  ps <- SYS_PROB(
    sys_input = list(lsf_emp, lsf_ln),
    probMachines = list(PROB_MACHINE(name = "SORM", fCall = "SORM")),
    debug.level = 0
  )
  suppressWarnings(ps$runMachines())

  beta_emp <- ps$beta_single["SORM", "emp"]
  beta_ln <- ps$beta_single["SORM", "ln"]

  expect_true(is.finite(beta_emp))
  expect_lt(abs(beta_emp - beta_ln), 0.05)
})
