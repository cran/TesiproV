library(testthat)

test_that("MC_IS_single reproduces FORM for linear LSF", {
    skip_on_cran()

    # -----------------------------
    # Basisvariablen
    # -----------------------------
    R <- PROB_BASEVAR(
        Id = 1, Name = "R",
        DistributionType = "norm",
        Mean = 10, Sd = 2
    )

    S <- PROB_BASEVAR(
        Id = 2, Name = "S",
        DistributionType = "norm",
        Mean = 8, Sd = 2
    )

    # -----------------------------
    # LSF
    # -----------------------------
    lsf <- SYS_LSF(vars = list(R, S), name = "Linear test")

    lsf$func <- function(R, S) {
        R - S
    }

    lsf$check()

    # -----------------------------
    # MC_IS Maschine
    # -----------------------------
    machine_MC_IS <- PROB_MACHINE(
        name = "MC_IS",
        fCall = "MC_IS",
        options = list(
            backend = "future",
            n_batch = 2000,
            n_max = 5e4,
            cov_user = 0.05,
            seed = 123
        )
    )

    # -----------------------------
    # FORM Maschine
    # -----------------------------
    machine_FORM <- PROB_MACHINE(
        name = "FORM",
        fCall = "FORM"
    )

    # -----------------------------
    # Problem aufsetzen
    # -----------------------------
    ps <- SYS_PROB(
        sys_input = list(lsf),
        probMachines = list(machine_FORM, machine_MC_IS),
        debug.level = 0
    )

    with_test_future_plan(future::sequential, {
        ps$runMachines()
    })

    beta_FORM <- ps$beta_single["FORM", 1]
    beta_MCIS <- ps$beta_single["MC_IS", 1]

    # -----------------------------
    # Test
    # -----------------------------
    expect_lt(abs(beta_MCIS - beta_FORM), 0.05)
})


test_that("MC_IS_single works with non-normal distributions", {
    skip_on_cran()

    # -----------------------------
    # Basisvariablen
    # -----------------------------
    Z <- PROB_BASEVAR(
        Id = 1, Name = "Z",
        DistributionType = "norm",
        Mean = 100, Cov = 0.04
    )

    Fy <- PROB_BASEVAR(
        Id = 2, Name = "Fy",
        DistributionType = "lnorm",
        Mean = 40, Cov = 0.1
    )

    M <- PROB_BASEVAR(
        Id = 3, Name = "M",
        DistributionType = "gumbel",
        Package = "evd",
        Mean = 2000, Cov = 0.1
    )

    # -----------------------------
    # LSF
    # -----------------------------
    lsf <- SYS_LSF(vars = list(Z, Fy, M), name = "Non-normal test")

    lsf$func <- function(Z, Fy, M) {
        Z * Fy - M
    }

    lsf$check()

    # -----------------------------
    # Maschinen
    # -----------------------------
    machine_FORM <- PROB_MACHINE(
        name = "FORM",
        fCall = "FORM"
    )

    machine_MC_IS <- PROB_MACHINE(
        name = "MC_IS",
        fCall = "MC_IS",
        options = list(
            backend = "future",
            n_batch = 5000,
            n_max = 1e5,
            cov_user = 0.05,
            seed = 123
        )
    )

    # -----------------------------
    # Problem aufsetzen
    # -----------------------------
    ps <- SYS_PROB(
        sys_input = list(lsf),
        probMachines = list(machine_FORM, machine_MC_IS),
        debug.level = 0
    )

    with_test_future_plan(future::sequential, {
        ps$runMachines()
    })

    beta_FORM <- ps$beta_single["FORM", 1]
    beta_MCIS <- ps$beta_single["MC_IS", 1]

    # -----------------------------
    # Test
    # -----------------------------
    expect_lt(abs(beta_MCIS - beta_FORM), 0.1)
})


test_that("MC_IS_single estimates rare-event probability correctly", {
    skip_on_cran()

    # -----------------------------
    # Basisvariablen
    # -----------------------------
    R <- PROB_BASEVAR(
        Id = 1, Name = "R",
        DistributionType = "norm",
        Mean = 100, Sd = 10
    )

    S <- PROB_BASEVAR(
        Id = 2, Name = "S",
        DistributionType = "norm",
        Mean = 70, Sd = 10
    )

    # -----------------------------
    # LSF
    # -----------------------------
    lsf <- SYS_LSF(vars = list(R, S), name = "Rare event test")

    lsf$func <- function(R, S) {
        R - S
    }

    lsf$check()

    # -----------------------------
    # Maschinen
    # -----------------------------
    machine_FORM <- PROB_MACHINE(
        name = "FORM",
        fCall = "FORM"
    )

    machine_MC_IS <- PROB_MACHINE(
        name = "MC_IS",
        fCall = "MC_IS",
        options = list(
            backend = "future",
            n_batch = 5000,
            n_max = 5e5,
            cov_user = 0.05,
            seed = 123
        )
    )

    # -----------------------------
    # Problem
    # -----------------------------
    ps <- SYS_PROB(
        sys_input = list(lsf),
        probMachines = list(machine_FORM, machine_MC_IS),
        debug.level = 0
    )

    with_test_future_plan(future::sequential, {
        ps$runMachines()
    })

    beta_FORM <- ps$beta_single["FORM", 1]
    beta_MCIS <- ps$beta_single["MC_IS", 1]

    # -----------------------------
    # Test
    # -----------------------------
    expect_lt(abs(beta_MCIS - beta_FORM), 0.1)
})


test_that("MC_IS_single is reproducible with Unix parallel workers", {
    skip_on_os("windows")

    x <- PROB_BASEVAR(
        Id = 1,
        Name = "X",
        DistributionType = "norm",
        Mean = 0,
        Sd = 1
    )

    lsf <- SYS_LSF(
        name = "Parallel single-LSF test",
        vars = list(x)
    )

    lsf$func <- function(X) {
        3 - X
    }

    run_once <- function() {
        ps <- SYS_PROB(
            sys_input = list(lsf),
            probMachines = list(
                PROB_MACHINE(
                    name = "MC-IS",
                    fCall = "MC_IS",
                    options = list(
                        backend = "parallel",
                        use_threads = 2L,
                        n_batch = 1000L,
                        n_max = 2000L,
                        cov_user = 0.5,
                        seed = 1234L,
                        dataRecord = FALSE
                    )
                )
            ),
            sys_type = "serial",
            debug.level = 0
        )

        ps$runMachines()
        ps$res_single[[1]][[1]]
    }

    res1 <- run_once()
    res2 <- run_once()

    expect_equal(res1$n_mc, res2$n_mc)
    expect_equal(res1$beta, res2$beta)
    expect_equal(res1$pf, res2$pf)
})


test_that("MC_IS_single keeps ultra-rare-event weights in the log domain", {
  x_var <- PROB_BASEVAR(
    Id = 1,
    Name = "X",
    DistributionType = "norm",
    Mean = 0,
    Sd = 1
  )

  # use_threads = 1 creates no cluster, so no warning must occur
  expect_no_warning(
    res <- MC_IS(
      lsf = function(x) 40 - x[1],
      lDistr = list(x_var$getlDistr()[[1]]),
      dps = 40,
      cov_user = 1.0,
      n_batch = 2000L,
      n_max = 2000L,
      use_threads = 1L,
      backend = "parallel",
      seed = 12345L,
      dataRecord = FALSE,
      stability_mode = "robust"
    )
  )

  expect_equal(res$n_mc, 2000L)

  expect_true(
    "log_pf" %in% names(res),
    info = "MC_IS result must expose log_pf for ultra-rare events."
  )

  expect_true(is.finite(res$log_pf))

  # beta must remain available even if exp(log_pf) underflows to zero.
  expect_true(is.finite(res$beta))
  expect_gt(res$beta, 35)
  expect_lt(res$beta, 45)
})


test_that("MC_IS_single is reproducible with two Future workers", {
    x_var <- PROB_BASEVAR(
        Id = 1,
        Name = "X",
        DistributionType = "norm",
        Mean = 0,
        Sd = 1
    )

    run_once <- function() {
        MC_IS(
            lsf = function(x) 4 - x[1],
            lDistr = list(x_var$getlDistr()[[1]]),
            dps = 4,
            backend = "future",
            cov_user = 0,
            n_batch = 1000L,
            n_max = 2000L,
            seed = 12345L,
            dataRecord = FALSE
        )
    }

    results <- with_test_future_plan(
        future::multisession,
        {
            expect_equal(
                as.integer(future::nbrOfWorkers()),
                2L
            )

            list(run_once(), run_once())
        },
        workers = 2
    )

    expect_equal(results[[1]]$n_mc, 2000L)
    expect_equal(results[[1]]$pf, results[[2]]$pf)
    expect_equal(results[[1]]$beta, results[[2]]$beta)
})
