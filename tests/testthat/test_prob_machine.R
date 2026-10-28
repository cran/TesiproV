test_that("PROB_MACHINE does not add options on its own", {
  machine <- PROB_MACHINE$new(name = "FORM", fCall = "FORM")

  expect_length(machine$options, 0)
})

test_that("PROB_MACHINE keeps the options given by the user", {
  machine <- PROB_MACHINE$new(name = "FORM", fCall = "FORM",
                              options = list(n_optim = 50L, loctol = 1e-6))

  expect_equal(machine$options$n_optim, 50L)
  expect_equal(machine$options$loctol, 1e-6)
})

test_that("PROB_MACHINE preserves explicit FORM options", {
  machine <- PROB_MACHINE(
    name = "FORM test",
    fCall = "FORM",
    options = list(
      n_optim = 50L,
      loctol = 1e-6
    )
  )

  expect_equal(machine$options$n_optim, 50L)
  expect_equal(machine$options$loctol, 1e-6)
})


test_that("PROB_MACHINE rejects unavailable methods", {
  machine <- PROB_MACHINE(
    name = "Unknown method",
    fCall = "NOT_A_REAL_METHOD"
  )

  expect_error(
    machine$checkMethodExists(),
    "Method 'NOT_A_REAL_METHOD' is not loaded",
    fixed = TRUE
  )
})


test_that("SYS_PROB forwards valid machine options to the algorithm", {
  received <- new.env(parent = emptyenv())

  mock_algorithm <- function(lsf,
                             lDistr,
                             custom_option = NULL,
                             debug.level = NULL) {
    received$custom_option <- custom_option
    received$debug.level <- debug.level

    list(
      method = "mock_algorithm",
      beta = 3.0,
      pf = stats::pnorm(-3.0),
      runtime = c(0, 0, 0, 0, 0)
    )
  }

  assign("mock_algorithm", mock_algorithm, envir = .GlobalEnv)

  on.exit(
    rm("mock_algorithm", envir = .GlobalEnv),
    add = TRUE
  )

  variable <- PROB_BASEVAR(
    Id = 1,
    Name = "X",
    DistributionType = "norm",
    Mean = 10,
    Sd = 1
  )

  lsf <- SYS_LSF(
    name = "Simple LSF",
    vars = list(variable)
  )

  lsf$func <- function(x) x[1] - 5

  machine <- PROB_MACHINE(
    name = "Mock method",
    fCall = "mock_algorithm",
    options = list(custom_option = 123)
  )

  ps <- SYS_PROB(
    sys_input = list(lsf),
    probMachines = list(machine),
    sys_type = "serial",
    debug.level = 0
  )

  ps$runMachines()

  expect_equal(received$custom_option, 123)
  expect_equal(received$debug.level, 0)

  expect_equal(ps$beta_single[1, 1], 3.0)
  expect_equal(ps$res_single[[1]][[1]]$method, "mock_algorithm")
})
