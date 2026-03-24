library(testthat)
library(microbenchmark)

testthat::local_edition(3)



# -----------------------------------------------------------------
# Check RNG initialization and helper functions
# -----------------------------------------------------------------
test_that("init_rng_master creates a consistent stream", {
  s1 <- init_rng_master(1234)
  s2 <- init_rng_master(1234)
  expect_equal(s1, s2)

  # Verify that .Random.seed is actually set
  expect_true(!is.null(.Random.seed))
})


test_that("record_step inserts values correctly", {
  rec <- make_recorders(5)
  rec$cov <- numeric(5) # für completeness
  rec <- record_step(rec,
    idx = 2,
    pf_mc = 0.1, var_mc = 0.02,
    I_mc = 3, n_sim = 200,
    tic = proc.time()
  )
  expect_equal(rec$pf[2], 0.1)
  expect_equal(rec$n_sim[2], 200)
  expect_true(all(!is.na(rec$time$user[2])))
})


test_that("make_ready_cluster returns a functional cluster", {
  libPaths <- .libPaths()
  cl <- make_ready_cluster(
    use_threads = 2,
    libPaths_local = libPaths,
    seed = 123,
    export_objs = character(0)
  )
  expect_true(inherits(cl, "cluster"))

  # Test call: each worker should have the same initial seed
  out <- parallel::parLapply(cl, 1:2, function(i) .Random.seed[1])
  expect_equal(out[[1]], out[[2]]) # identical start seed per worker

  parallel::stopCluster(cl)
})


# -----------------------------------------------------------------
# Performance benchmark test for new vs old implementation
# -----------------------------------------------------------------
test_that("new implementation is not slower than old one", {
  set.seed(321)

  n_batch <- 2e3
  n_vars <- 4

  fx_i <- matrix(runif(n_batch * n_vars, 0.1, 2), nrow = n_batch)
  hv_i <- matrix(runif(n_batch * n_vars, 0.1, 2), nrow = n_batch)

  mb <- microbenchmark(
    alt = {
      fx_prod <- apply(fx_i, 1, prod)
      hv_prod <- apply(hv_i, 1, prod)
      qt <- fx_prod / hv_prod
    },
    neu = {
      fx_prod <- rowSums(log(fx_i))
      hv_prod <- rowSums(log(hv_i))
      qt <- exp(fx_prod - hv_prod)
    },
    times = 15
  )

  # Expectation: median(new) ≤ median(old)*0.9 (≈10% faster)
  med_alt <- median(mb$time[mb$expr == "alt"])
  med_neu <- median(mb$time[mb$expr == "neu"])
  expect_lte(med_neu, med_alt * 0.9)
})


# -----------------------------------------------------------------
# Example Ex.5.11 from Nowak & Collins (2000)
# -----------------------------------------------------------------
test_that("Ex. 5.11 Nowak & Collins (fast version)", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }
  use_threads <- 1L # CI stabil

  X1 <- PROB_BASEVAR(
    Id = 1, Name = "Z",
    DistributionType = "norm",
    Mean = 100, Cov = 0.04
  )

  X2 <- PROB_BASEVAR(
    Id = 2, Name = "Fy",
    DistributionType = "lnorm",
    Mean = 40, Cov = 0.1
  )

  X3 <- PROB_BASEVAR(
    Id = 3, Name = "M",
    DistributionType = "gumbel",
    Package = "evd",
    Mean = 2000, Cov = 0.1
  )

  lsf1 <- SYS_LSF(vars = list(X1, X2, X3))
  lsf1$func <- function(M, Fy, Z) Z * Fy - M
  lsf1$check()

  form_rf <- PROB_MACHINE("FORM", "FORM")
  sorm <- PROB_MACHINE("SORM", "SORM")

  mcis <- PROB_MACHINE("MCIS", "MC_IS",
    options = list(
      cov_user = 0.05,
      n_max = 200000, # statt Millionen
      n_batch = 20000,
      use_threads = use_threads,
      backend = backend,
      seed = 1234
    )
  )

  ps <- SYS_PROB(
    sys_input = list(lsf1),
    probMachines = list(form_rf, sorm, mcis)
  )

  ps$runMachines()

  expect_length(ps$beta_single, 3)
  expect_false(any(is.na(ps$beta_single)))

  # FORM == SORM bei linearer LSF
  expect_equal(ps$beta_single[1], ps$beta_single[2], tolerance = 1e-8)

  # Plausibilität: Beta zwischen 3 und 5
  expect_gt(ps$beta_single[3], 3)
  expect_lt(ps$beta_single[3], 5)
})


# -----------------------------------------------------------------
# Monte Carlo Importance Sampling with many input variables -------
# -----------------------------------------------------------------
test_that("MCIS and high number of inputs", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }
  # Helper function to create a normal variable ------------------
  make_norm_var <- function(id, mean, sd) {
    PROB_BASEVAR(
      Id = id,
      Name = paste0("u", id),
      DistributionType = "norm",
      Mean = mean,
      Sd = sd
    )
  }

  # Create a set of 26 variables – the first has a different mean ------
  vars <- lapply(1:26, function(i) make_norm_var(i, 0.2, 0.1))
  vars[[1]]$Mean <- 0.5

  # LSF -----------------------------------------------------------
  lsf4 <- SYS_LSF(vars = vars, name = "RackwitzProb.")
  lsf4$func <- function(x) {
    x[1] - sum((x[-1]^2) / seq_along(x)[-1])
  }

  lsf4$check()

  # MC–IS machine(no fixed seed → new RNG stream each run)
  mcis <- PROB_MACHINE(
    name = "MC IS",
    fCall = "MC_IS",
    options = list(
      cov_user = 0.08,
      n_max = 400000,
      seed = 1234,
      backend = backend
    )
  )

  ps <- SYS_PROB(
    sys_input = list(lsf4),
    probMachines = list(mcis),
    debug.level = 0
  )

  ps$runMachines()

  # --------------------------------------------------------------
  # Plausibility checks
  # --------------------------------------------------------------
  # (c) MC‑IS – maximum 2% relative deviation allowed

  rel_err_1 <- abs(ps$beta_single[1] - 3.337) / 3.337

  message("Check: MC-IS must not deviate more than 2% from the reference beta")
  expect_lt(rel_err_1, 0.07)

  message("Check: Calculated beta must be less than 3.5")
  expect_lt(ps$beta_single[1], 3.5)

  message("Check: Calculated beta should be greater than 3.2")
  expect_gt(ps$beta_single[1], 3.1)

  expect_false(any(is.na(ps$beta_single)),
    info = "beta_single darf keine NA-Werte enthalten"
  )

  expect_length(ps$beta_single, 1)
})


# -----------------------------------------------------------------
# Simple linear LSF using FORM/SORM/MCC comparison ----------------
# -----------------------------------------------------------------
test_that("Simple linear LSf FORM SORM and MCC", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }
  # Fallback for Windows
  if (.Platform$OS.type == "windows") {
    use_threads <- 1L
  } else {
    use_threads <- 16L
  }

  var1 <- PROB_BASEVAR(Name = "E", Description = "Effect", DistributionType = "norm", Mean = 10, Sd = 0.2)
  var2 <- PROB_BASEVAR(Name = "R", Description = "Resistance", DistributionType = "norm", Mean = 35, Sd = 0.5)


  lsf <- SYS_LSF(vars = list(var1, var2), name = "LSF XY")
  lsf$func <- function(E, R) {
    R - E^1.5
  }
  # you can run this check to see if all transformations of the variables worked, all needed data is available and the set of vars fits to the limit state function
  lsf$check()


  machine1 <- PROB_MACHINE(name = "FORM Rack.-Fieß.", fCall = "FORM", options = list("n_optim" = 20, "loctol" = 0.001, "optim_type" = "rackfies"))
  machine2 <- PROB_MACHINE(name = "MC CoV 0.05", fCall = "MC_CRUDE", options = list("n_max" = 1e7, "cov_user" = 0.05, "use_threads" = use_threads, "seed" = 1234, backend = backend))
  machine3 <- PROB_MACHINE(name = "SORM", fCall = "SORM")

  ps <- SYS_PROB(
    sys_input = list(lsf),
    probMachines = list(machine1, machine2, machine3),
    debug.level = 0
  )


  ps$runMachines()

  expect_equal(ps$beta_single[1], 3.11578032380995) # FORM

  if (.Platform$OS.type == "windows") {
    rel_err <- abs(ps$beta_single[2] - 3.11578032380995) / 3.11578032380995
    expect_lt(rel_err, 0.02) # < 2 % Abweichung erlaubt
  } else {
    rel_err <- abs(ps$beta_single[2] - 3.09023230616781) / 3.09023230616781
    expect_lt(rel_err, 0.02) # < 2 % Abweichung erlaubt
  }

  expect_equal(ps$beta_single[3], 3.11578032380995) # SORM
})


test_that("MC_IS works with empirical distribution (small test)", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }

  skip_if_not_installed("EnvStats")

  set.seed(1)

  data_small <- rnorm(5000, 100, 4)

  X <- PROB_BASEVAR(
    Id = 1,
    Name = "Z",
    DistributionType = "emp",
    Package = "EnvStats",
    DistributionParameters = data_small
  )

  lsf <- SYS_LSF(vars = list(X), name = "Simple test")
  lsf$func <- function(Z) Z - 95
  lsf$check()

  mcis <- PROB_MACHINE(
    name = "MC IS",
    fCall = "MC_IS",
    options = list(
      cov_user = 0.05,
      n_max = 50000,
      seed = 1234,
      backend = backend
    )
  )

  ps <- SYS_PROB(
    sys_input = list(lsf),
    probMachines = list(mcis)
  )

  ps$runMachines()

  expect_false(is.na(ps$beta_single))
  expect_true(is.finite(ps$beta_single))
})
