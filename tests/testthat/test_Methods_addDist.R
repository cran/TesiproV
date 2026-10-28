# -----------------------------------------------------------------
# Example Ex.5.11 from Nowak & Collins (2000): reliability analysis
# with student-t and log-student-t
# -----------------------------------------------------------------
test_that("Ex. 5.11 Nowak & Collins (2000)", {
  skip_on_cran()
  testthat::local_edition(3)

  # The backend is irrelevant for this example. With the sequential plan
  # from setup.R, backend = "future" runs identically on all systems.
  backend <- "future"

  # Base variables ------------------------------------------------
  X1 <- PROB_BASEVAR(
    Id = 1, Name = "Z",
    DistributionType = "st",
    Package = "TesiproV",
    DistributionParameters = c(100, 4, 1000000000, 1000000000)
  ) # kN/m^2

  sig_fy <- sqrt(log(1 + 0.1^2))
  mue_fy <- log(40) - sig_fy^2 / 2

  X2 <- PROB_BASEVAR(
    Id = 2, Name = "Fy",
    DistributionType = "lt",
    Package = "TesiproV",
    DistributionParameters = c(mue_fy, sig_fy, 1000000000, 1000000000)
  ) # m

  X3 <- PROB_BASEVAR(
    Id = 3, Name = "M",
    DistributionType = "gumbel",
    Package = "evd",
    Mean = 2000, Cov = 0.1
  ) # m

  # Limit State Function -------------------------------------------
  lsf1 <- SYS_LSF(
    vars = list(X1, X2, X3),
    name = "Nowak & Collins Exp5.11_1 (beta=4.03)"
  )

  lsf1$func <- function(M, Fy, Z) Z * Fy - M

  lsf1$check()

  # Machines -------------------------------------------------------
  form_rf <- PROB_MACHINE(
    name = "FORM Rack.-Fiess.",
    fCall = "FORM",
    options = list(
      n_optim = 20,
      loctol = 1e-4,
      optim_type = "rackfies"
    )
  )

  sorm <- PROB_MACHINE(name = "SORM", fCall = "SORM")

  mcis <- PROB_MACHINE(
    name = "MC IS",
    fCall = "MC_IS",
    options = list(
      cov_user = 0.02,
      n_max = 1e7,
      n_batch = 100000,
      seed = 1234,
      backend = backend
    )
  )

  mcc <- PROB_MACHINE(
    name = "MCC",
    fCall = "MC_CRUDE",
    options = list(
      cov_user = 0.05,
      n_max = 1e8,
      seed = 161268,
      backend = backend
    )
  )

  mcsus <- PROB_MACHINE(
    name = "MC SuS",
    fCall = "MC_SubSam",
    options = list(seed = 1234)
  )

  # System ---------------------------------------------------------
  ps <- SYS_PROB(
    sys_input = list(lsf1),
    probMachines = list(form_rf, sorm, mcis, mcc, mcsus),
    debug.level = 0
  )

  ps$runMachines()

  # --------------------------------------------------------------
  # Basic plausibility checks
  # --------------------------------------------------------------
  expect_length(ps$beta_single, 5)
  expect_false(any(is.na(ps$beta_single)))

  beta_exact <- 4.009179

  # FORM: known linearisation error, but close to the exact solution
  expect_equal(ps$beta_single[1], 4.02211455568222, tolerance = 1e-8)
  expect_lt(abs(ps$beta_single[1] - beta_exact) / beta_exact, 0.005)

  # SORM: g = Z * Fy - M is nonlinear, SORM must be closer to the
  # exact solution than FORM
  expect_lt(abs(ps$beta_single[2] - beta_exact),
            abs(ps$beta_single[1] - beta_exact))
  expect_lt(abs(ps$beta_single[2] - beta_exact), 0.005)

  # MC-IS: maximum 2 % relative deviation
  rel_err_1 <- abs(ps$beta_single[3] - 4.03) / 4.03
  expect_lt(rel_err_1, 0.021)

  # MC-CRUDE: maximum 2 % relative deviation
  ref_mc_crude <- 4.02211455568236
  tol_mc_crude <- 0.02
  rel_err_2 <- abs(ps$beta_single[4] - ref_mc_crude) / ref_mc_crude
  expect_lt(rel_err_2, tol_mc_crude)

  # MC-SuS: literature value
  expect_equal(ps$beta_single[5], 4.0, tolerance = 5e-3)
})
