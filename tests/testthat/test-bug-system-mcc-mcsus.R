# System analysis with crude Monte Carlo ("MCC") and subset simulation
# ("MCSUS") on the system limit-state function.
#
# Reference solutions:
#  - independent functions g_i = R_i - S_i (normal): p_i = Phi(-beta_i),
#    serial p = 1 - prod(1 - p_i), parallel p = prod(p_i)
#  - common load S in both functions: given S = s the resistances are
#    independent, p follows by numerical integration over s.

test_that("MCC-1: serial system with independent functions", {
  skip_on_cran()

  R1 <- PROB_BASEVAR$new(Id = 1, Name = "R1", DistributionType = "norm", Mean = 2, Sd = 0.5)
  S1 <- PROB_BASEVAR$new(Id = 2, Name = "S1", DistributionType = "norm", Mean = 1, Sd = 0.5)
  R2 <- PROB_BASEVAR$new(Id = 3, Name = "R2", DistributionType = "norm", Mean = 2.3, Sd = 0.5)
  S2 <- PROB_BASEVAR$new(Id = 4, Name = "S2", DistributionType = "norm", Mean = 1, Sd = 0.5)
  g1 <- SYS_LSF$new(name = "g1", vars = list(R1, S1))
  g1$func <- function(x) x[1] - x[2]
  g2 <- SYS_LSF$new(name = "g2", vars = list(R2, S2))
  g2$func <- function(x) x[1] - x[2]

  p_exact <- 1 - (1 - pnorm(-1 / sqrt(0.5))) * (1 - pnorm(-1.3 / sqrt(0.5)))
  beta_exact <- -qnorm(p_exact)

  set.seed(1)
  ps <- SYS_PROB$new(sys_input = list(g1, g2), probMachines = list(), sys_type = "serial")

  backend_mcc <- if (.Platform$OS.type == "windows") {
    "future"
  } else {
    "parallel"
  }

  mcc_params <- list(
    cov_user = 0.05,
    n_batch = 10000,
    n_max = 2e6,
    backend = backend_mcc
  )

  if (identical(backend_mcc, "parallel")) {
    mcc_params$use_threads <- 2L
  }

  if (identical(backend_mcc, "future")) {
    suppressMessages(
      with_test_future_plan(
        future::multisession,
        ps$calculateSystemProbability(
          calcType = "MCC",
          params = mcc_params
        ),
        workers = 2L
      )
    )
  } else {
    suppressMessages(
      ps$calculateSystemProbability(
        calcType = "MCC",
        params = mcc_params
      )
    )
  }

  res <- ps$res_sys[[length(ps$res_sys)]]

  expect_lt(abs(res$beta - beta_exact), 0.05)
  expect_equal(rownames(ps$beta_sys)[nrow(ps$beta_sys)], "MC Crude")
})

test_that("MCC-2: parallel system with independent functions", {
  skip_on_cran()

  R1 <- PROB_BASEVAR$new(Id = 1, Name = "R1", DistributionType = "norm", Mean = 2, Sd = 0.5)
  S1 <- PROB_BASEVAR$new(Id = 2, Name = "S1", DistributionType = "norm", Mean = 1, Sd = 0.5)
  R2 <- PROB_BASEVAR$new(Id = 3, Name = "R2", DistributionType = "norm", Mean = 2.3, Sd = 0.5)
  S2 <- PROB_BASEVAR$new(Id = 4, Name = "S2", DistributionType = "norm", Mean = 1, Sd = 0.5)
  g1 <- SYS_LSF$new(name = "g1", vars = list(R1, S1))
  g1$func <- function(x) x[1] - x[2]
  g2 <- SYS_LSF$new(name = "g2", vars = list(R2, S2))
  g2$func <- function(x) x[1] - x[2]

  p_exact <- pnorm(-1 / sqrt(0.5)) * pnorm(-1.3 / sqrt(0.5))
  beta_exact <- -qnorm(p_exact)

  set.seed(1)
  ps <- SYS_PROB$new(sys_input = list(g1, g2), probMachines = list(), sys_type = "parallel")
  backend_mcc <- if (.Platform$OS.type == "windows") {
    "future"
  } else {
    "parallel"
  }

  mcc_params <- list(
    cov_user = 0.05,
    n_batch = 10000,
    n_max = 2e6,
    backend = backend_mcc
  )

  if (identical(backend_mcc, "parallel")) {
    mcc_params$use_threads <- 2L
  }

  if (identical(backend_mcc, "future")) {
    suppressMessages(
      with_test_future_plan(
        future::multisession,
        ps$calculateSystemProbability(
          calcType = "MCC",
          params = mcc_params
        ),
        workers = 2L
      )
    )
  } else {
    suppressMessages(
      ps$calculateSystemProbability(
        calcType = "MCC",
        params = mcc_params
      )
    )
  }
  res <- ps$res_sys[[length(ps$res_sys)]]

  expect_lt(abs(res$beta - beta_exact), 0.08)
})

test_that("MCC-3: a variable shared by two functions enters the system only once", {
  skip_on_cran()

  R1 <- PROB_BASEVAR$new(Id = 1, Name = "R1", DistributionType = "norm", Mean = 2, Sd = 0.5)
  R2 <- PROB_BASEVAR$new(Id = 2, Name = "R2", DistributionType = "norm", Mean = 2.3, Sd = 0.5)
  S <- PROB_BASEVAR$new(Id = 3, Name = "S", DistributionType = "norm", Mean = 1, Sd = 0.5)
  g1 <- SYS_LSF$new(name = "g1", vars = list(R1, S))
  g1$func <- function(x) x[1] - x[2]
  g2 <- SYS_LSF$new(name = "g2", vars = list(R2, S))
  g2$func <- function(x) x[1] - x[2]

  p_serial <- integrate(
    function(s) dnorm(s, 1, 0.5) * (1 - (1 - pnorm((s - 2) / 0.5)) * (1 - pnorm((s - 2.3) / 0.5))),
    lower = -10, upper = 10, rel.tol = 1e-10
  )$value
  p_parallel <- integrate(
    function(s) dnorm(s, 1, 0.5) * pnorm((s - 2) / 0.5) * pnorm((s - 2.3) / 0.5),
    lower = -10, upper = 10, rel.tol = 1e-10
  )$value

  # Create the serial and parallel system objects before running either analysis.
  ps_s <- SYS_PROB$new(
    sys_input = list(g1, g2),
    probMachines = list(),
    sys_type = "serial"
  )

  ps_p <- SYS_PROB$new(
    sys_input = list(g1, g2),
    probMachines = list(),
    sys_type = "parallel"
  )

  # Select the backend for this test.
  backend_mcc <- if (.Platform$OS.type == "windows") {
    "future"
  } else {
    "parallel"
  }

  mcc_params <- list(
    cov_user = 0.05,
    n_batch = 10000,
    n_max = 2e6,
    backend = backend_mcc
  )

  if (identical(backend_mcc, "parallel")) {
    mcc_params$use_threads <- 2L
  }

  set.seed(1)

  # Run both system analyses with the selected backend.
  if (identical(backend_mcc, "future")) {
    suppressMessages(
      with_test_future_plan(
        future::multisession,
        {
          ps_s$calculateSystemProbability(
            calcType = "MCC",
            params = mcc_params
          )

          ps_p$calculateSystemProbability(
            calcType = "MCC",
            params = mcc_params
          )
        },
        workers = 2L
      )
    )
  } else {
    suppressMessages(
      ps_s$calculateSystemProbability(
        calcType = "MCC",
        params = mcc_params
      )
    )

    suppressMessages(
      ps_p$calculateSystemProbability(
        calcType = "MCC",
        params = mcc_params
      )
    )
  }

  expect_lt(abs(ps_s$res_sys[[1]]$beta - (-qnorm(p_serial))), 0.05)
  # Treating S as two independent copies would give beta of about 2.8 here,
  # the exact value is about 2.25.
  expect_lt(abs(ps_p$res_sys[[1]]$beta - (-qnorm(p_parallel))), 0.08)
})

test_that("MCSUS-1: subset simulation for serial and parallel systems with a common load", {
  skip_on_cran()

  R1 <- PROB_BASEVAR$new(Id = 1, Name = "R1", DistributionType = "norm", Mean = 2, Sd = 0.5)
  R2 <- PROB_BASEVAR$new(Id = 2, Name = "R2", DistributionType = "norm", Mean = 2.3, Sd = 0.5)
  S <- PROB_BASEVAR$new(Id = 3, Name = "S", DistributionType = "norm", Mean = 1, Sd = 0.5)
  g1 <- SYS_LSF$new(name = "g1", vars = list(R1, S))
  g1$func <- function(x) x[1] - x[2]
  g2 <- SYS_LSF$new(name = "g2", vars = list(R2, S))
  g2$func <- function(x) x[1] - x[2]

  p_serial <- integrate(
    function(s) dnorm(s, 1, 0.5) * (1 - (1 - pnorm((s - 2) / 0.5)) * (1 - pnorm((s - 2.3) / 0.5))),
    lower = -10, upper = 10, rel.tol = 1e-10
  )$value
  p_parallel <- integrate(
    function(s) dnorm(s, 1, 0.5) * pnorm((s - 2) / 0.5) * pnorm((s - 2.3) / 0.5),
    lower = -10, upper = 10, rel.tol = 1e-10
  )$value

  set.seed(1)
  ps_s <- SYS_PROB$new(sys_input = list(g1, g2), probMachines = list(), sys_type = "serial")
  suppressMessages(ps_s$calculateSystemProbability(
    calcType = "MCSUS", params = list(Nsubset = 5000, p0 = 0.1)
  ))
  ps_p <- SYS_PROB$new(sys_input = list(g1, g2), probMachines = list(), sys_type = "parallel")
  suppressMessages(ps_p$calculateSystemProbability(
    calcType = "MCSUS", params = list(Nsubset = 5000, p0 = 0.1)
  ))

  expect_lt(abs(ps_s$res_sys[[1]]$beta - (-qnorm(p_serial))), 0.15)
  expect_lt(abs(ps_p$res_sys[[1]]$beta - (-qnorm(p_parallel))), 0.15)
  expect_equal(rownames(ps_s$beta_sys)[nrow(ps_s$beta_sys)], "MC SubSam")
})
