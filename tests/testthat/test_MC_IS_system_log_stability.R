# The mock replaces only the accumulation step. This lets the tests
# exercise the system estimator's post-processing with known log sums,
# including values whose exponentials underflow to zero.
.run_mc_is_system_with_mocked_log_sums <- function(
  log_pf,
  log_second_moment,
  n_batch = 20L
) {
    mock_compute_weight_update <- function(all_log_w,
                                           all_I_vals,
                                           stability_mode = NULL,
                                           sys_type = NULL,
                                           state) {
        n_new <- state$n_sim + length(all_log_w)

        state$log_sum_Iw <- if (is.finite(log_pf)) {
            log(n_new) + log_pf
        } else {
            -Inf
        }

        state$log_sum_Iw2 <- if (is.finite(log_second_moment)) {
            log(n_new) + log_second_moment
        } else {
            -Inf
        }

        state$n_sim <- n_new
        state
    }

    testthat::local_mocked_bindings(
        compute_weight_update = mock_compute_weight_update,
        .package = "TesiproV"
    )

    v <- PROB_BASEVAR(
        Id = 1,
        Name = "X",
        DistributionType = "norm",
        Mean = 0,
        Sd = 1
    )
    marginal <- v$getlDistr()[[1]]

    lsf <- function(x) {
        if (is.matrix(x)) {
            1 - x[, 1L]
        } else {
            1 - x[1L]
        }
    }

    MC_IS(
        lsf = list(lsf),
        lDistr = list(list(marginal)),
        sys_type = "serial",
        backend = "parallel",
        use_threads = 1L,
        dps = 1,
        seed = 123L,
        n_batch = n_batch,
        n_max = n_batch,
        cov_user = 0,
        dataRecord = FALSE,
        debug.level = 0
    )
}


test_that("system importance-sampling log sums remain stable at extreme scales", {
    state <- list(
        log_sum_Iw = -Inf,
        log_sum_Iw2 = -Inf,
        n_sim = 0L
    )

    update <- TesiproV:::compute_weight_update(
        all_log_w = c(-1000, -1001),
        all_I_vals = c(1L, 1L),
        state = state
    )

    expect_equal(
        update$log_sum_Iw,
        -1000 + log1p(exp(-1))
    )
    expect_equal(
        update$log_sum_Iw2,
        -2000 + log1p(exp(-2))
    )
    expect_identical(update$n_sim, 2L)
})


test_that("system estimator agrees with ordinary-scale results when representable", {
    n <- 20L
    p <- 0.2

    res <- .run_mc_is_system_with_mocked_log_sums(
        log_pf = log(p),
        log_second_moment = log(p),
        n_batch = n
    )

    expected_var <- (p - p^2) / n
    expected_cov <- sqrt(expected_var) / p

    expected_log_var <- log(expected_var)

    expect_equal(res$log_var, expected_log_var, tolerance = 1e-14)

    expect_equal(res$pf, p, tolerance = 1e-14)
    expect_equal(res$log_pf, log(p), tolerance = 1e-14)
    expect_equal(res$beta, -stats::qnorm(p), tolerance = 1e-14)
    expect_equal(res$var, expected_var, tolerance = 1e-14)
    expect_equal(res$cov_mc, expected_cov, tolerance = 1e-14)
    expect_equal(res$ESS, n * p, tolerance = 1e-14)
})


test_that("system estimator preserves finite log probability after pf underflows", {
    n <- 20L
    log_pf <- -1000
    log_second_moment <- -1990

    res <- .run_mc_is_system_with_mocked_log_sums(
        log_pf = log_pf,
        log_second_moment = log_second_moment,
        n_batch = n
    )

    log_var <- log_second_moment +
        log1p(-exp(2 * log_pf - log_second_moment)) -
        log(n)

    expected_cov <- exp(0.5 * log_var - log_pf)
    expected_ess <- exp(log(n) + 2 * log_pf - log_second_moment)

    expect_identical(res$pf, 0)
    expect_equal(res$log_pf, log_pf, tolerance = 1e-14)
    expect_true(is.finite(res$beta))
    expect_equal(
        res$beta,
        -stats::qnorm(log_pf, log.p = TRUE),
        tolerance = 1e-12
    )
    expect_true(is.finite(res$cov_mc))
    expect_equal(res$cov_mc, expected_cov, tolerance = 1e-12)
    expect_true(is.finite(res$ESS))
    expect_gt(res$ESS, 0)
    expect_lte(res$ESS, res$n_mc)
    expect_equal(res$ESS, expected_ess, tolerance = 1e-12)

    expect_equal(res$log_var, log_var, tolerance = 1e-12)
    expect_identical(res$var, 0)
})


test_that("system estimator reports no observed failure mass consistently", {
    res <- .run_mc_is_system_with_mocked_log_sums(
        log_pf = -Inf,
        log_second_moment = -Inf,
        n_batch = 20L
    )

    expect_identical(res$pf, 0)
    expect_identical(res$log_pf, -Inf)
    expect_true(is.na(res$beta))
    expect_identical(res$ESS, 0)

    expect_identical(res$log_var, NA_real_)
    expect_identical(res$var, NA_real_)
})
