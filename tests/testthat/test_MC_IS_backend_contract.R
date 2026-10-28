test_that("MC_IS warns if use_threads is supplied with backend = 'future'", {
  oplan <- future::plan(future::sequential)
  withr::defer(future::plan(oplan))

  x_var <- PROB_BASEVAR(
    Id = 1,
    Name = "X",
    DistributionType = "norm",
    Mean = 0,
    Sd = 1
  )

  expect_warning(
    MC_IS(
      lsf = function(x) 3 - x[1],
      lDistr = list(x_var$getlDistr()[[1]]),
      n_batch = 1000L,
      n_max = 1000L,
      use_threads = 1L,
      backend = "future",
      seed = 1L,
      dataRecord = FALSE
    ),
    "backend = 'future' ignores 'use_threads = 1' and uses the active Future plan"
  )
})


test_that("MC_IS defaults to parallel with one worker", {
  observed <- new.env(parent = emptyenv())

  # MC_IS_single is replaced here, so that only the backend selection
  # of the public MC_IS() function is tested.
  stub_mc_is_single <- function(...) {
    args <- list(...)
    observed$backend <- args$backend
    observed$use_threads <- args$use_threads
    list(method = "stub")
  }

  testthat::local_mocked_bindings(
    MC_IS_single = stub_mc_is_single,
    .package = "TesiproV"
  )

  result <- MC_IS(
    lsf = function(x) x[1],
    lDistr = list(list()),
    n_batch = 10L,
    n_max = 10L,
    seed = 1L
    # backend and use_threads intentionally not supplied
  )

  expect_identical(result$method, "stub")
  expect_identical(observed$backend, "parallel")
  expect_identical(observed$use_threads, 1L)
})


test_that("MC_IS with future respects the active Future plan", {
  observed <- new.env(parent = emptyenv())

  stub_mc_is_single <- function(...) {
    args <- list(...)
    observed$backend <- args$backend
    observed$use_threads <- args$use_threads
    list(method = "stub")
  }

  testthat::local_mocked_bindings(
    MC_IS_single = stub_mc_is_single,
    .package = "TesiproV"
  )

  with_test_future_plan(future::multisession,
    {
      plan_before <- future::plan()
      workers_before <- future::nbrOfWorkers()

      expect_equal(as.integer(workers_before), 2L)

      # An explicit use_threads must not override the worker count of the
      # Future plan and is reported as ignored.
      expect_warning(
        result <- MC_IS(
          lsf = function(x) x[1],
          lDistr = list(list()),
          backend = "future",
          use_threads = 9L,
          n_batch = 10L,
          n_max = 10L,
          seed = 1L
        ),
        "uses the active Future plan with 2 worker"
      )

      expect_identical(result$method, "stub")
      expect_identical(observed$backend, "future")
      expect_identical(observed$use_threads, 2L)

      # TesiproV must not change the plan set by the user.
      expect_identical(future::plan(), plan_before)
      expect_equal(as.integer(future::nbrOfWorkers()), workers_before)
    },
    workers = 2L
  )
})


test_that("parallel backend creates and uses a PSOCK cluster on Windows", {
  if (.Platform$OS.type != "windows") {
    skip("PSOCK-specific check; run on Windows.")
  }

  cl <- make_parallel_cluster(use_threads = 2L, seed = 1L)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  expect_true(inherits(cl, "SOCKcluster"))

  result <- run_parallel(
    X = list(1L, 2L),
    FUN = function(x) x + 1L,
    backend = "parallel",
    use_threads = 2L,
    cl = cl
  )

  expect_equal(unname(unlist(result)), c(2L, 3L))
})


test_that("MC_IS uses the Future plan, not use_threads, for chunking", {
  x_var <- PROB_BASEVAR(
    Id = 1,
    Name = "X",
    DistributionType = "norm",
    Mean = 0,
    Sd = 1
  )
  marginal <- x_var$getlDistr()[[1]]

  lsf <- function(x) {
    if (is.matrix(x)) {
      4 - x[, 1]
    } else {
      4 - x[1]
    }
  }

  run_case <- function(supply_use_threads = FALSE) {
    args <- list(
      lsf = lsf,
      lDistr = list(marginal),
      backend = "future",
      dps = 4,
      seed = 12345L,
      n_batch = 4000L,
      n_max = 4000L,
      cov_user = 0,
      dataRecord = FALSE,
      debug.level = 0
    )

    if (supply_use_threads) {
      args$use_threads <- 8L
    }

    do.call(MC_IS, args)
  }

  results <- with_test_future_plan(
    future::multisession,
    {
      expect_equal(as.integer(future::nbrOfWorkers()), 2L)

      without_use_threads <- run_case()

      expect_warning(
        with_use_threads <- run_case(supply_use_threads = TRUE),
        "active Future plan with 2 worker"
      )

      list(
        without_use_threads = without_use_threads,
        with_use_threads = with_use_threads
      )
    },
    workers = 2L
  )

  fields_to_compare <- c(
    "method",
    "beta",
    "pf",
    "log_pf",
    "var",
    "cov_mc",
    "ESS",
    "ESS_ratio",
    "n_mc",
    "design_points"
  )

  expect_identical(
    results$with_use_threads[fields_to_compare],
    results$without_use_threads[fields_to_compare]
  )
})


test_that("MC_IS_single is reproducible across backends with the same seed and worker count", {
  x_var <- PROB_BASEVAR(
    Id = 1,
    Name = "X",
    DistributionType = "norm",
    Mean = 0,
    Sd = 1
  )

  marginal <- x_var$getlDistr()[[1]]

  # Failure occurs when X >= 4, support both vector and matrix input.
  lsf <- function(x) {
    if (is.matrix(x)) {
      4 - x[, 1L]
    } else {
      4 - x[1L]
    }
  }

  common_args <- list(
    lsf = lsf,
    lDistr = list(marginal),
    dps = 4,
    seed = 314159L,
    cov_user = 0,
    n_batch = 8192L,
    n_max = 8192L,
    dataRecord = FALSE,
    debug.level = 0L
  )

  # The parallel backend takes its worker count from use_threads.
  res_parallel <- do.call(
    MC_IS,
    c(
      common_args,
      list(
        backend = "parallel",
        use_threads = 2L
      )
    )
  )

  # The Future backend takes its worker count from the active Future plan.
  res_future <- with_test_future_plan(
    future::multisession,
    {
      expect_equal(
        as.integer(future::nbrOfWorkers()),
        2L
      )

      do.call(
        MC_IS,
        c(
          common_args,
          list(backend = "future")
        )
      )
    },
    workers = 2L
  )

  expect_equal(res_parallel$n_mc, 8192L)
  expect_equal(res_future$n_mc, 8192L)

  # Compare the probability estimates and diagnostics using the same
  # fixed seed, sample budget, and worker count.
  fields_to_compare <- c(
    "pf",
    "log_pf",
    "beta",
    "var",
    "cov_mc",
    "ESS"
  )

  for (field in fields_to_compare) {
    expect_equal(
      res_parallel[[field]],
      res_future[[field]],
      tolerance = 1e-14,
      info = paste("Mismatch in result field:", field)
    )
  }
})