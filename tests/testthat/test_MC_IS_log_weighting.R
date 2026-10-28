# tests/testthat/test_MC_IS_log_weighting.R

test_that("MC_IS_single: Log-Wichtung stimmt mit direkter Berechnung überein", {
    x_var <- PROB_BASEVAR(
        Id = 1,
        Name = "X",
        DistributionType = "norm",
        Mean = 0,
        Sd = 1
    )

    marginal <- x_var$getlDistr()[[1]]

    # g(x) <= 0 genau dann, wenn x >= 4.
    lsf <- function(x) 4 - x[1]

    seed <- 12345L
    n <- 4000L
    dps <- 4

    res <- with_test_future_plan(future::sequential, {
        MC_IS(
            lsf = lsf,
            lDistr = list(marginal),
            dps = dps,
            cov_user = 0, # nur bei Erreichen von n_max stoppen
            n_batch = n,
            n_max = n,
            backend = "future",
            seed = seed,
            dataRecord = FALSE,
            stability_mode = "robust"
        )
    })

    # Dieselben standardnormalen Zufallszahlen wie im Single-LSF-Lauf erzeugen.
    rng <- create_rng_manager(seed)
    z <- rng$generate_norm_block(n)

    # Bei dps = 4 ist die Vorschlagsverteilung N(4, 4^2).
    u <- dps + 4 * z

    log_w <- stats::dnorm(u, log = TRUE) -
        stats::dnorm(u, mean = dps, sd = 4, log = TRUE)

    failure <- as.numeric(vapply(
        u,
        function(ui) lsf(c(ui)) <= 0,
        logical(1)
    ))
    weighted_failure <- failure * exp(log_w)

    pf_ref <- mean(weighted_failure)
    second_moment_ref <- mean(weighted_failure^2)

    var_ref <- max(second_moment_ref - pf_ref^2, 0) / n
    cov_ref <- sqrt(var_ref) / pf_ref
    beta_ref <- -stats::qnorm(pf_ref)

    expect_equal(res$n_mc, n)
    expect_gt(pf_ref, 0)

    expect_equal(res$pf, pf_ref, tolerance = 1e-10)
    expect_equal(res$log_pf, log(pf_ref), tolerance = 1e-10)
    expect_equal(res$var, var_ref, tolerance = 1e-10)
    expect_equal(res$cov_mc, cov_ref, tolerance = 1e-10)
    expect_equal(res$beta, beta_ref, tolerance = 1e-10)
})


test_that("MC_IS_single: n_max wird auch bei Teil-Batches eingehalten", {
    x_var <- PROB_BASEVAR(
        Id = 1,
        Name = "X",
        DistributionType = "norm",
        Mean = 0,
        Sd = 1
    )

    marginal <- x_var$getlDistr()[[1]]

    # Kein Versagen: Der CoV-Abbruch greift nicht, daher muss n_max stoppen.
    lsf <- function(x) 1

    res <- with_test_future_plan(future::sequential, {
        MC_IS(
            lsf = lsf,
            lDistr = list(marginal),
            dps = 3,
            cov_user = 0.05,
            n_batch = 1000L,
            n_max = 2500L,
            backend = "future",
            seed = 123L,
            dataRecord = FALSE
        )
    })

    expect_equal(res$n_mc, 2500L)
    expect_equal(res$pf, 0)
    expect_identical(res$log_pf, -Inf)
    expect_true(is.infinite(res$cov_mc))
})


test_that("MC_IS_single: n_batch größer als n_max überschreitet n_max nicht", {
    x_var <- PROB_BASEVAR(
        Id = 1,
        Name = "X",
        DistributionType = "norm",
        Mean = 0,
        Sd = 1
    )

    marginal <- x_var$getlDistr()[[1]]

    res <- with_test_future_plan(future::sequential, {
        MC_IS(
            lsf = function(x) 1,
            lDistr = list(marginal),
            dps = 3,
            cov_user = 0.05,
            n_batch = 5000L,
            n_max = 2000L,
            backend = "future",
            seed = 123L,
            dataRecord = FALSE
        )
    })

    expect_equal(res$n_mc, 2000L)
})
