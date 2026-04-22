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
            n_batch = 2000,
            n_max = 5e4,
            cov_user = 0.05,
            seed = 123,
            use_threads = 1
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

    ps$runMachines()

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
            n_batch = 5000,
            n_max = 1e5,
            cov_user = 0.05,
            seed = 123,
            use_threads = 1
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

    ps$runMachines()

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
            n_batch = 5000,
            n_max = 1e5,
            cov_user = 0.05,
            seed = 123,
            use_threads = 1
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

    ps$runMachines()

    beta_FORM <- ps$beta_single["FORM", 1]
    beta_MCIS <- ps$beta_single["MC_IS", 1]

    # -----------------------------
    # Test
    # -----------------------------
    expect_lt(abs(beta_MCIS - beta_FORM), 0.1)
})
