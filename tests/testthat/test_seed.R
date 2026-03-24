# tests/testthat/test-seed.R -----------------------------------------

#library(testthat)

# Load the development version of the package (only needed when the
# tests are executed in the source tree, not when the package is
# installed).
# if (requireNamespace("devtools", quietly = TRUE)) {
#   devtools::load_all("..")    # path to the package root
# }

# --------------------------------------------------------------
# Simple 3-RV limit-state function (the one you already use)
# --------------------------------------------------------------
make_demo_lsf <- function() {
  X1 <- PROB_BASEVAR(Id = 1, Name = "Z",
                     DistributionType = "norm", Mean = 100, Cov = 0.04)
  X2 <- PROB_BASEVAR(Id = 2, Name = "Fy",
                     DistributionType = "lnorm", Mean = 40, Cov = 0.10)
  X3 <- PROB_BASEVAR(Id = 3, Name = "M",
                     DistributionType = "gumbel", Package = "evd",
                     Mean = 2000, Cov = 0.10)

  lsf <- SYS_LSF(vars = list(X1, X2, X3),
                 name = "Demo LSF")
  lsf$func <- function(M, Fy, Z) Z * Fy - M
  lsf$check()
  lsf
}

# --------------------------------------------------------------
# Wrapper that starts one MC-IS run with the *parallel* backend.
# --------------------------------------------------------------
run_one_mc_is <- function(use_threads = 2,
                          seed = NULL,
                          n_batch = 2e3,
                          n_max   = 5e4,
                          backend = "parallel") {

  # 1) build the machine (no fixed seed if seed == NULL)
  machine_mc_is <- PROB_MACHINE(
    name    = "MC IS (test)",
    fCall   = "MC_IS",
    options = list(
      n_max       = n_max,
      cov_user    = 0.02,
      use_threads = use_threads,
      n_batch     = n_batch,
      backend     = backend,
      seed        = seed,          # <- the argument we are testing
      dataRecord  = FALSE,
      densityType = "norm"
    )
  )

  # 2) build the SYS_PROB object that holds everything
  ps <- SYS_PROB(
    sys_input    = list(make_demo_lsf()),
    probMachines = list(machine_mc_is),
    debug.level  = 0
  )

  # 3) run the machines and return the *single* beta value
  ps$runMachines()
  # `beta_single` is a matrix, we pick the first (and only) entry
  as.numeric(ps$beta_single)
}



# -----------------------------------------------------------------
# 1)  Random seed (seed = NULL) - beta values must *vary*
# -----------------------------------------------------------------
test_that("seed = NULL gives non‑identical beta values", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }
  # Run the simulation twice with the SAME arguments but with seed = NULL
  beta1 <- run_one_mc_is(seed = NULL, backend = backend)
  beta2 <- run_one_mc_is(seed = NULL, backend = backend)

  # The two results must be *different* (within a reasonable tolerance)
  # – because a new master RNG stream is created for every call.
  expect_false(isTRUE(all.equal(beta1, beta2, tolerance = .Machine$double.eps)))

  # As an additional sanity check we run the simulation a few more times
  # and verify that the empirical variance is > 0.
  betas <- replicate(5, run_one_mc_is(seed = NULL, backend = backend ))
  expect_gt(var(betas), 0)
})

# -----------------------------------------------------------------
# 2)  Fixed seed (seed = 1234) - beta values must be reproducible
# -----------------------------------------------------------------
test_that("fixed seed yields identical beta values", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }
  fixed_seed <- 1234L

  # Two independent calls with the *same* fixed seed
  beta_a <- run_one_mc_is(seed = fixed_seed, backend = backend)
  beta_b <- run_one_mc_is(seed = fixed_seed, backend = backend)

  # They have to be exactly the same (up to machine precision)
  expect_equal(beta_a, beta_b, tolerance = .Machine$double.eps)

  # And a third call should also match – this catches accidental
  # side‑effects (e.g. the seed being consumed somewhere else).
  beta_c <- run_one_mc_is(seed = fixed_seed, backend = backend)
  expect_equal(beta_a, beta_c, tolerance = .Machine$double.eps)
})
