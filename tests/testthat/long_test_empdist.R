
# -----------------------------------------------------------------
# Example Ex.5.11 from Nowak & Collins (2000): reliability analysis with student-t and log-student-t
# -----------------------------------------------------------------
test_that("Ex. 5.11 Nowak & Collins (2000)", {
  testthat::local_edition(3)

  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }

  # Platform-dependent thread count
  if (.Platform$OS.type == "windows") {
    use_threads <- 1L
  } else {
    use_threads <- 20L
  }

  # Base variables ------------------------------------------------
  X1 <- PROB_BASEVAR(
    Id = 1, Name = "Z",
    DistributionType = "st",
    DistributionParameters = c(100,4,10000,10000)
  ) # kN/m²

  sig_fy <- sqrt(log(1+0.1^2))
  mue_fy <- log(40) - sig_fy^2/2
  data_fy <- rlnorm(1000000,mue_fy,sig_fy)
  X2 <- PROB_BASEVAR(
    Id = 2, Name = "Fy",
    DistributionType = "emp",
    Package = "EnvStats",
    DistributionParameters = data_fy)
  #Mean = 40, Cov = 0.1
  # X2 <- PROB_BASEVAR(
  #   Id = 2, Name = "Fy",
  #   DistributionType = "lt",
  #   DistributionParameters = c(mue_fy,sig_fy,10000,10000)
  # ) # m
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
      use_threads = 1L,
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
      use_threads = use_threads,
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

  ## FORM vs SORM – nearly identical results -----------------------
  rel_diff_form_sorm <- abs(ps$beta_single[1] - ps$beta_single[2]) /
    max(abs(ps$beta_single[1]), 1e-12)

  message("FORM and SORM should be identical for linear LSFs")
  expect_lt(rel_diff_form_sorm, 1e-8)

  ## Literature reference values -----------------------------------
  message("FORM beta must match literature value within 1e-8")
  expect_equal(ps$beta_single[1],
               4.02211455568222,
               tolerance = 1e-8
  )

  message("SORM beta must match literature value within 1e-8")
  expect_equal(ps$beta_single[2],
               4.02211455568236,
               tolerance = 1e-8
  )

  ## MC‑IS - maximum 2% relative deviation allowed
  rel_err_1 <- abs(ps$beta_single[3] - 3.9585857628075) / 3.9585857628075
  message("MC-IS must not deviate more than 2 % from the reference beta")
  expect_lt(rel_err_1, 0.02)

  ## MC‑CRUDE - maximum 2% relative deviation allowed
  ref_mc_crude <- 3.97628577254636 # 3.9160810558165
  tol_mc_crude <- 0.02 #  ≈ 2%

  rel_err_2 <- abs(ps$beta_single[4] - ref_mc_crude) / ref_mc_crude

  message(sprintf(
    "MC-CRUDE (%s) must lie within %.1f%% tolerance range",
    .Platform$OS.type, tol_mc_crude * 100
  ))

  expect_lt(rel_err_2, tol_mc_crude)
  # info=sprintf("MC-CRUDE (%s) must lie within %.1f%% tolerance range",
  #            .Platform$OS.type,tol_mc_crude*100))

  ## MC‑SuS – exact literature value -------------------------------
  message("MC-SuS beta must match literature value within 1e-8")
  expect_equal(ps$beta_single[5],
               4.01123480350189,
               tolerance = 1e-8
  )
})
