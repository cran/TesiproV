library(testthat)

test_that("PROB_BASEVAR: prepare and getlDistr produce funlist", {
    v <- PROB_BASEVAR(Id = 1, Name = "f_ck", DistributionType = "norm", Mean = 30, Sd = 1)

    # Accept both S3 and S4
    expect_true(inherits(v, "PROB_BASEVAR"))

    dd <- v$getlDistr()
    # dd kann wrapper (list(funlist, params)) oder funlist sein -> normalize
    funlist <- if (is.list(dd) && length(dd) >= 1 && is.list(dd[[1]]) && !is.null(dd[[1]]$p)) dd[[1]] else dd

    expect_true(is.list(funlist))
    expect_true(is.function(funlist$d))
    expect_true(is.function(funlist$p))
    expect_true(is.function(funlist$q))
    expect_true(is.function(funlist$r))
})


test_that("SYS_LSF getLSF evaluates correctly", {
    v1 <- PROB_BASEVAR(Id = 1, Name = "X1", DistributionType = "norm", Mean = 0, Sd = 1)

    lsf_obj <- SYS_LSF(name = "g1", vars = list(v1))
    lsf_obj$func <- function(X1) 3 - X1

    f <- lsf_obj$getLSF()
    expect_true(is.function(f))

    expect_equal(f(c(2)), 1)
    expect_equal(f(c(4)), -1)
})


library(testthat)

test_that("MC_IS_system: deterministic with fixed seed and dps (quick smoke test)", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }

    v1 <- PROB_BASEVAR(Id = 1, Name = "X1", DistributionType = "norm", Mean = 0, Sd = 1)
    v2 <- PROB_BASEVAR(Id = 2, Name = "X2", DistributionType = "norm", Mean = 0, Sd = 1)

    lsf1 <- function(x) 3 - x[1]
    lsf2 <- function(x) 3 - x[1]

    # get funlist (normalize)
    d1 <- v1$getlDistr()
    d1 <- if (is.list(d1) && length(d1) >= 1 && is.list(d1[[1]]) && !is.null(d1[[1]]$p)) d1[[1]] else d1
    d2 <- v2$getlDistr()
    d2 <- if (is.list(d2) && length(d2) >= 1 && is.list(d2[[1]]) && !is.null(d2[[1]]$p)) d2[[1]] else d2

    lDistr <- list(list(d1), list(d2))

    run_once <- function(seed_val) {
        MC_IS(
            lsf = list(lsf1, lsf2),
            lDistr = lDistr,
            cov_user = 0.5,
            n_batch = 1000,
            n_max = 2000,
            use_threads = 1,
            backend = backend,
            seed = seed_val,
            dps = c(3, 3),
            dataRecord = TRUE,
            debug.level = 0
        )
    }

    r1 <- run_once(1234)
    r2 <- run_once(1234)

    expect_equal(r1$pf, r2$pf)
    expect_equal(r1$n_mc, r2$n_mc)

    expect_true(is.numeric(r1$pf))
    expect_true(r1$pf >= 0 && r1$pf <= 1)

    # Test SYS_PROB wrapper
    lsf1obj <- SYS_LSF(name = "g1", vars = list(v1))
    lsf1obj$func <- function(X1) 3 - X1
    lsf2obj <- SYS_LSF(name = "g2", vars = list(v2))
    lsf2obj$func <- function(X2) 3 - X2
    ps <- SYS_PROB(sys_input = list(lsf1obj, lsf2obj), probMachines = list(), sys_type = "serial")
    ps$calculateSystemProbability(calcType = "MCIS", params = list(
        cov_user = 0.5, n_batch = 1000, n_max = 2000, use_threads = 1,
        seed = 1234, backend = "future", dps = c(3, 3), dataRecord = TRUE
    ))
    expect_true(length(ps$res_sys) >= 1)
    last <- ps$res_sys[[length(ps$res_sys)]]
    expect_true(!is.null(last$pf))
})


test_that("MC_IS_single is reproducible for same seed", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }
    v <- PROB_BASEVAR(
        Id = 1, Name = "X",
        DistributionType = "norm",
        Mean = 0, Sd = 1
    )

    lsf <- SYS_LSF(name = "g", vars = list(v))
    lsf$func <- function(X) 3 - X

    ps1 <- SYS_PROB(
        sys_input = list(lsf),
        probMachines = list(),
        sys_type = "serial"
    )

    ps1$calculateSystemProbability(
        calcType = "MCIS",
        params = list(
            cov_user = 0.01,
            n_batch = 5000,
            n_max = 1e5,
            seed = 1234,
            backend = backend
        )
    )

    beta1 <- ps1$beta_sys[1, 1]

    ps2 <- SYS_PROB(
        sys_input = list(lsf),
        probMachines = list(),
        sys_type = "serial"
    )

    ps2$calculateSystemProbability(
        calcType = "MCIS",
        params = list(
            cov_user = 0.01,
            n_batch = 5000,
            n_max = 1e5,
            seed = 1234,
            backend = backend
        )
    )

    beta2 <- ps2$beta_sys[1, 1]

    expect_equal(beta1, beta2)
})


test_that("MC_IS_system is reproducible for same seed", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }

    v1 <- PROB_BASEVAR(
        Id = 1, Name = "X1",
        DistributionType = "norm",
        Mean = 0, Sd = 1
    )

    v2 <- PROB_BASEVAR(
        Id = 2, Name = "X2",
        DistributionType = "norm",
        Mean = 0, Sd = 1
    )

    lsf1 <- function(x) 3 - x[1]
    lsf2 <- function(x) 3 - x[1]

    d1 <- v1$getlDistr()[[1]]
    d2 <- v2$getlDistr()[[1]]

    lDistr <- list(list(d1), list(d2))

    run_once <- function(seed_val) {
        MC_IS(
            lsf = list(lsf1, lsf2),
            lDistr = lDistr,
            cov_user = 0.5,
            n_batch = 1000,
            n_max = 2000,
            use_threads = 1,
            backend = backend,
            seed = seed_val,
            dps = c(3, 3),
            dataRecord = TRUE
        )
    }

    r1 <- run_once(1234)
    r2 <- run_once(1234)

    expect_equal(r1$pf, r2$pf)
    expect_equal(r1$n_mc, r2$n_mc)
})


test_that("ESS is positive and bounded", {

    v <- PROB_BASEVAR(
        Id = 1, Name = "X",
        DistributionType = "norm",
        Mean = 0, Sd = 1
    )

    lsf <- function(x) 3 - x[1]

    d1 <- v$getlDistr()[[1]]

    res <- MC_IS(
        lsf = list(lsf),
        lDistr = list(list(d1)),
        cov_user = 0.2,
        n_batch = 2000,
        n_max = 5000,
        use_threads = 1,
        backend = "future",
        seed = 1234
    )

    expect_true(res$ESS > 0)
    expect_true(res$ESS <= res$n_mc)
})


test_that("MC_IS_system handles intersecting limit states", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }

    v1 <- PROB_BASEVAR(
        Id = 1, Name = "X1",
        DistributionType = "norm",
        Mean = 0, Sd = 1
    )

    v2 <- PROB_BASEVAR(
        Id = 2, Name = "X2",
        DistributionType = "norm",
        Mean = 0, Sd = 1
    )

    lsf1 <- function(x) 3 - x[1]
    lsf2 <- function(x) 3 - (0.7 * x[1] + 0.7 * x[2])

    d1 <- v1$getlDistr()[[1]]
    d2 <- v2$getlDistr()[[1]]

    res <- MC_IS(
        lsf = list(lsf1, lsf2),
        lDistr = list(list(d1), list(d1, d2)),
        cov_user = 0.05,
        n_batch = 5000,
        n_max = 1e5,
        use_threads = 1,
        backend = backend,
        seed = 1234
    )

    expect_true(res$pf > 0)
    expect_true(res$pf < 1)
})


test_that("MC_IS serial is reproducible for fixed seed", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }

    v <- PROB_BASEVAR(
        Id = 1, Name = "X",
        DistributionType = "norm",
        Mean = 0, Sd = 1
    )

    lsf <- SYS_LSF(name = "g", vars = list(v))
    lsf$func <- function(X) 3 - X

    ps1 <- SYS_PROB(
        sys_input = list(lsf),
        probMachines = list(),
        sys_type = "serial"
    )

    ps2 <- SYS_PROB(
        sys_input = list(lsf),
        probMachines = list(),
        sys_type = "serial"
    )

    ps1$calculateSystemProbability(
        calcType = "MCIS",
        params = list(
            cov_user = 0.02,
            n_batch = 5000,
            n_max = 1e5,
            seed = 1234,
            backend = backend
        )
    )

    ps2$calculateSystemProbability(
        calcType = "MCIS",
        params = list(
            cov_user = 0.02,
            n_batch = 5000,
            n_max = 1e5,
            seed = 1234,
            backend = backend
        )
    )

    expect_equal(
        ps1$res_sys[[1]]$pf,
        ps2$res_sys[[1]]$pf
    )
})


test_that("adaptive alpha is reproducible for fixed seed", {
  skip_if_not_installed("future")
  old_plan <- future::plan()
  future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)

    v1 <- PROB_BASEVAR(
        Id = 1, Name = "X1",
        DistributionType = "norm",
        Mean = 0, Sd = 1
    )

    v2 <- PROB_BASEVAR(
        Id = 2, Name = "X2",
        DistributionType = "norm",
        Mean = 0, Sd = 1
    )

    lsf1 <- function(x) 3 - x[1]
    lsf2 <- function(x) 3 - (0.7 * x[1] + 0.7 * x[2])

    d1 <- v1$getlDistr()[[1]]
    d2 <- v2$getlDistr()[[1]]

    lDistr <- list(list(d1), list(d1, d2))

    run_once <- function(seed_val) {
        MC_IS(
            lsf = list(lsf1, lsf2),
            lDistr = lDistr,
            cov_user = 0.02,
            n_batch = 10000,
            n_max = 2e5,
            use_threads = 1,
            backend = "future",
            seed = seed_val,
            adaptive_alpha = TRUE
        )
    }

    r1 <- run_once(1234)
    r2 <- run_once(1234)

    expect_equal(r1$pf, r2$pf)
})


test_that("adaptive alpha does not catastrophically worsen Pf", {

    v1 <- PROB_BASEVAR(
        Id = 1, Name = "X1",
        DistributionType = "norm",
        Mean = 0, Sd = 1
    )

    v2 <- PROB_BASEVAR(
        Id = 2, Name = "X2",
        DistributionType = "norm",
        Mean = 0, Sd = 1
    )

    lsf1 <- function(x) 3 - x[1]
    lsf2 <- function(x) 3 - (0.7 * x[1] + 0.7 * x[2])

    d1 <- v1$getlDistr()[[1]]
    d2 <- v2$getlDistr()[[1]]

    lDistr <- list(list(d1), list(d1, d2))

    static_run <- MC_IS(
        lsf = list(lsf1, lsf2),
        lDistr = lDistr,
        cov_user = 0.02,
        n_batch = 10000,
        n_max = 2e5,
        use_threads = 1,
        backend = "future",
        seed = 1234,
        adaptive_alpha = FALSE
    )

    adaptive_run <- MC_IS(
        lsf = list(lsf1, lsf2),
        lDistr = lDistr,
        cov_user = 0.02,
        n_batch = 10000,
        n_max = 2e5,
        use_threads = 1,
        backend = "future",
        seed = 1234,
        adaptive_alpha = TRUE
    )

    expect_lt(
        abs(adaptive_run$pf - static_run$pf),
        0.5 * static_run$pf
    )
})


test_that("ESS is finite and positive for adaptive alpha", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }

    v1 <- PROB_BASEVAR(
        Id = 1, Name = "X1",
        DistributionType = "norm",
        Mean = 0, Sd = 1
    )

    v2 <- PROB_BASEVAR(
        Id = 2, Name = "X2",
        DistributionType = "norm",
        Mean = 0, Sd = 1
    )

    lsf1 <- function(x) 3 - x[1]
    lsf2 <- function(x) 3 - (0.7 * x[1] + 0.7 * x[2])

    d1 <- v1$getlDistr()[[1]]
    d2 <- v2$getlDistr()[[1]]

    lDistr <- list(list(d1), list(d1, d2))

    res <- MC_IS(
        lsf = list(lsf1, lsf2),
        lDistr = lDistr,
        cov_user = 0.02,
        n_batch = 10000,
        n_max = 2e5,
        use_threads = 1,
        backend = backend,
        seed = 1234,
        adaptive_alpha = TRUE
    )

    expect_true(is.finite(res$ESS))
    expect_gt(res$ESS, 0)
})


test_that("compute_weight_update fast mode accumulates correctly", {
    state <- list(
        total_sum_w = 0,
        total_sum_Iw = 0,
        total_sum_w2 = 0,
        log_total_sum_w = NULL,
        log_total_sum_Iw = NULL,
        log_total_sum_w2 = NULL,
        n_sim = 0
    )

    all_log_w <- c(log(2), log(1), log(3))
    all_I_vals <- c(1, 0, 1)

    res <- compute_weight_update(
        all_log_w = all_log_w,
        all_I_vals = all_I_vals,
        stability_mode = "fast",
        sys_type = "serial",
        state = state
    )

    expect_true(res$total_sum_w > 0)
    expect_equal(res$n_sim, 3)
    expect_true(res$total_sum_Iw <= res$total_sum_w)
    expect_true(res$total_sum_w2 > 0)
})


test_that("compute_weight_update robust mode accumulates log-sums", {
    state <- list(
        total_sum_w = 0,
        total_sum_Iw = 0,
        total_sum_w2 = 0,
        log_total_sum_w = -Inf,
        log_total_sum_Iw = -Inf,
        log_total_sum_w2 = -Inf,
        n_sim = 0
    )

    all_log_w <- c(-1, -2, -0.5)
    all_I_vals <- c(1, 0, 1)

    res <- compute_weight_update(
        all_log_w = all_log_w,
        all_I_vals = all_I_vals,
        stability_mode = "robust",
        sys_type = "parallel",
        state = state
    )

    expect_true(is.finite(res$log_total_sum_w))
    expect_equal(res$n_sim, 3)
})
