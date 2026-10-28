test_that("u_to_x preserves extreme tails for normal marginals", {
    marginal <- list(
        DistributionType = "norm",
        mean = 0,
        Sd = 1,
        q = stats::qnorm
    )

    u <- c(-40, 0, 40)

    expect_equal(
        u_to_x(u, marginal),
        u
    )
})


test_that("u_to_x applies normal mean and standard deviation directly", {
    marginal <- list(
        DistributionType = "norm",
        mean = 10,
        Sd = 2,
        q = function(p) stats::qnorm(p, mean = 10, sd = 2)
    )

    u <- c(-40, 0, 40)

    expect_equal(
        u_to_x(u, marginal),
        c(-70, 10, 90)
    )
})


test_that("u_to_x uses the generic quantile transformation for non-normal marginals", {
    marginal <- list(
        DistributionType = "lnorm",
        mean = 1,
        Sd = 0.2,
        q = function(p) {
            stats::qlnorm(
                p,
                meanlog = 0,
                sdlog = 0.2
            )
        }
    )

    u <- c(-1, 0, 1)

    expected <- marginal$q(stats::pnorm(u))

    expect_equal(
        u_to_x(u, marginal),
        expected
    )
})
