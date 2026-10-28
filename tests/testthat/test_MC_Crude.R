test_that("MC_CRUDE does not change the active future plan", {
  oplan <- future::plan(future::sequential)
  withr::defer(future::plan(oplan))

  x_var <- PROB_BASEVAR(Id = 1, Name = "X", DistributionType = "norm", Mean = 0, Sd = 1)

  MC_CRUDE(
    lsf = function(x) 2 - x[1],
    lDistr = list(x_var$getlDistr()),
    cov_user = 0.1,
    n_batch = 1000L,
    n_max = 1000L,
    backend = "future",
    seed = 1L,
    dataRecord = FALSE
  )

  expect_true(inherits(future::plan(), "sequential"))
})

test_that("MC_CRUDE does not exceed n_max", {
  x_var <- PROB_BASEVAR(Id = 1, Name = "X", DistributionType = "norm", Mean = 0, Sd = 1)

  res <- MC_CRUDE(
    lsf = function(x) 2 - x[1],
    lDistr = list(x_var$getlDistr()),
    cov_user = 1e-6,
    n_batch = 300L,
    n_max = 1000L,
    backend = "future",
    seed = 1L,
    dataRecord = FALSE
  )

  expect_equal(res$n_mc, 1000)
})

test_that("MC_CRUDE result does not depend on the number of workers", {
  x_var <- PROB_BASEVAR(Id = 1, Name = "X", DistributionType = "norm", Mean = 0, Sd = 1)

  res_1 <- MC_CRUDE(
    lsf = function(x) 2 - x[1],
    lDistr = list(x_var$getlDistr()),
    cov_user = 1e-6,
    n_batch = 1000L,
    n_max = 5000L,
    use_threads = 1L,
    backend = "parallel",
    seed = 42L,
    dataRecord = FALSE
  )

  res_2 <- MC_CRUDE(
    lsf = function(x) 2 - x[1],
    lDistr = list(x_var$getlDistr()),
    cov_user = 1e-6,
    n_batch = 1000L,
    n_max = 5000L,
    use_threads = 2L,
    backend = "parallel",
    seed = 42L,
    dataRecord = FALSE
  )

  expect_identical(res_1$pf, res_2$pf)
  expect_identical(res_1$n_mc, res_2$n_mc)
})

test_that("MC_CRUDE warns if use_threads is supplied with backend = 'future'", {
  oplan <- future::plan(future::sequential)
  withr::defer(future::plan(oplan))

  x_var <- PROB_BASEVAR(Id = 1, Name = "X", DistributionType = "norm", Mean = 0, Sd = 1)

  expect_warning(
    MC_CRUDE(
      lsf = function(x) 2 - x[1],
      lDistr = list(x_var$getlDistr()),
      cov_user = 0.1,
      n_batch = 1000L,
      n_max = 1000L,
      use_threads = 1L,
      backend = "future",
      seed = 1L,
      dataRecord = FALSE
    ),
    "use_threads"
  )
})
