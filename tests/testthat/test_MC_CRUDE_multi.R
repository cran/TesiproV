library(testthat)

test_that("MC_CRUDE defaults to parallel with appropriate fallback cores and warns on future conflict", {
    observed <- new.env(parent = emptyenv())

    # Setup basic dummy random variable matching the true SYS_PROB environment structure
    x_var <- PROB_BASEVAR(Id = 1, Name = "X", DistributionType = "norm", Mean = 0, Sd = 1)
    marginal_wrapped <- x_var$getlDistr()

    # Determine active workers in the current test runner environment (usually 1 for sequential)
    current_workers <- as.integer(future::nbrOfWorkers())

    # Intentionally pass a conflicting thread count to force and verify the architectural warning
    expect_warning(
        res <- MC_CRUDE(
            lsf = function(x) 1 - x,
            lDistr = list(marginal_wrapped),
            backend = "future",
            use_threads = 80L,
            n_batch = 10L,
            n_max = 20L,
            cov_user = 0,
            dataRecord = FALSE
        ),
        "ignores 'use_threads = 80' and uses the active Future plan"
    )

    expect_identical(res$method, "MC_CRUDE")
    expect_equal(res$n_threads, current_workers)
})


test_that("MC_CRUDE with future backend respects the active global plan", {
    x_var <- PROB_BASEVAR(Id = 1, Name = "X", DistributionType = "norm", Mean = 0, Sd = 1)
    marginal_wrapped <- x_var$getlDistr()

    # Use the package utility to enforce a temporary 2-worker plan
    with_test_future_plan(future::multisession,
        {
            plan_before <- future::plan()
            workers_before <- future::nbrOfWorkers()

            expect_equal(as.integer(workers_before), 2L)

            # Wrap expect_warning around the entire assignment block to preserve the return object
            expect_warning(
                res <- MC_CRUDE(
                    lsf = function(x) 4 - x,
                    lDistr = list(marginal_wrapped),
                    backend = "future",
                    use_threads = 8L,
                    n_batch = 100L,
                    n_max = 200L,
                    cov_user = 0,
                    dataRecord = FALSE
                ),
                "ignores 'use_threads = 8' and uses the active Future plan with 2 workers"
            )

            expect_identical(res$method, "MC_CRUDE")
            expect_equal(res$n_threads, 2L)

            # Verify the global workspace environment has not been corrupted
            expect_identical(future::plan(), plan_before)
            expect_equal(as.integer(future::nbrOfWorkers()), workers_before)
        },
        workers = 2L
    )
})


test_that("MC_CRUDE uses the Future plan when use_threads is omitted", {
    x_var <- PROB_BASEVAR(
        Id = 1,
        Name = "X",
        DistributionType = "norm",
        Mean = 0,
        Sd = 1
    )
    marginal_wrapped <- x_var$getlDistr()

    with_test_future_plan(
        future::multisession,
        {
            expect_equal(as.integer(future::nbrOfWorkers()), 2L)

            expect_no_warning(
                res <- MC_CRUDE(
                    lsf = function(x) 1, # stets sicher: keine Versagensereignisse
                    lDistr = list(marginal_wrapped),
                    backend = "future",
                    n_batch = 10L,
                    n_max = 10L,
                    cov_user = 0.05,
                    dataRecord = FALSE,
                    seed = 123L
                    # use_threads absichtlich weggelassen
                )
            )

            expect_identical(res$n_threads, 2L)
        },
        workers = 2L
    )
})


test_that("MC_CRUDE execution is reproducible across backends under strict configurations", {
    x_var <- PROB_BASEVAR(Id = 1, Name = "X", DistributionType = "norm", Mean = 0, Sd = 1)
    marginal_wrapped <- x_var$getlDistr()

    lsf <- function(x) {
        if (is.matrix(x)) 4 - x[, 1L] else 4 - x[1L]
    }

    common_args <- list(
        lsf = lsf,
        lDistr = list(marginal_wrapped),
        seed = 42L,
        cov_user = 0,
        n_batch = 2000L,
        n_max = 2000L,
        dataRecord = FALSE,
        debug.level = 0L
    )

    # Standard calculation via native parallel backend (explicitly drives 2 threads)
    res_parallel <- do.call(MC_CRUDE, c(common_args, list(backend = "parallel", use_threads = 2L)))

    # Target calculation via active future plan
    # Execute the warning check strictly wrapped around the assignment to prevent NULL conversion
    with_test_future_plan(future::multisession,
        {
            expect_warning(
                res_future <<- do.call(MC_CRUDE, c(common_args, list(backend = "future", use_threads = 80L))),
                "ignores 'use_threads = 80' and uses the active Future plan with 2 workers"
            )
        },
        workers = 2L
    )

    # n_batch is the total number of samples per batch, independent of the workers
    expect_equal(res_parallel$n_mc, 2000L)
    expect_equal(res_future$n_mc, 2000L)

    # Compare critical fields to ensure identical cross-platform stochastic generation
    fields_to_compare <- c("pf", "beta", "var", "cov_mc")
    for (field in fields_to_compare) {
        expect_equal(
            res_parallel[[field]],
            res_future[[field]],
            tolerance = 1e-14,
            info = paste("Stochastic mismatch found in field:", field)
        )
    }
})


test_that("MC_CRUDE uses the Future plan without warning when use_threads is omitted", {
    x_var <- PROB_BASEVAR(
        Id = 1,
        Name = "X",
        DistributionType = "norm",
        Mean = 0,
        Sd = 1
    )
    marginal <- x_var$getlDistr()

    with_test_future_plan(
        future::multisession,
        {
            expect_no_warning(
                res <- MC_CRUDE(
                    lsf = function(x) 4 - x,
                    lDistr = list(marginal),
                    backend = "future",
                    n_batch = 10L,
                    n_max = 20L,
                    cov_user = 0,
                    dataRecord = FALSE
                    # use_threads absichtlich weggelassen
                )
            )

            expect_equal(res$n_threads, 2L)
        },
        workers = 2L
    )
})
