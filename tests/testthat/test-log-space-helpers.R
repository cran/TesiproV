test_that("log_sum_exp avoids underflow for very negative values", {
    log_w <- c(-900, -901, -905)

    expect_equal(sum(exp(log_w)), 0)

    result <- log_sum_exp(log_w)

    expect_true(is.finite(result))
    expect_equal(
        result,
        -899.6818245707525,
        tolerance = 1e-12
    )
})


test_that("log_add_exp combines log-scale values stably", {
    a <- -900
    b <- -901

    result <- log_add_exp(a, b)

    expect_true(is.finite(result))
    expect_equal(
        result,
        -899.6867383124818,
        tolerance = 1e-12
    )
})


test_that("log_sum_exp returns minus infinity for empty input", {
    expect_identical(log_sum_exp(numeric()), -Inf)
})


test_that("log_add_exp treats minus infinity as log zero", {
    expect_identical(log_add_exp(-Inf, -123), -123)
    expect_identical(log_add_exp(-123, -Inf), -123)
})
