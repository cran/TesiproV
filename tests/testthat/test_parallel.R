library(testthat)

test_that("make_parallel_cluster works on current OS", {
    cl <- make_parallel_cluster(use_threads = 2)
    expect_true(inherits(cl, "cluster"))
    parallel::stopCluster(cl)
})


library(testthat)

test_that("run_parallel works with future (default plan)", {
    chunks <- split(1:10, cut(1:10, 2))

    out <- run_parallel(
        X = chunks,
        FUN = function(v) sum(v),
        backend = "future",
        use_threads = 2
    )

    expect_equal(
        unlist(out),
        unlist(lapply(chunks, sum))
    )
})


test_that("RNG manager reproducibility", {
    rng1 <- create_rng_manager(seed = 123)
    rng2 <- create_rng_manager(seed = 123)

    block1 <- rng1$generate_norm_block(1000)
    block2 <- rng2$generate_norm_block(1000)

    expect_equal(block1, block2)
})


test_that("RNG manager advances state correctly", {
    rng <- create_rng_manager(seed = 123)

    block1 <- rng$generate_norm_block(10)
    block2 <- rng$generate_norm_block(10)

    set.seed(123)
    full <- rnorm(20)

    expect_equal(block1, full[1:10])
    expect_equal(block2, full[11:20])
})
