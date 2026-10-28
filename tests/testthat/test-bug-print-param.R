# Verification of the suspected issues in printResults() and SYS_PARAM
# Run with: devtools::test(filter = "bug-print-param")

test_that("PRT-1: SYS_PROB$printResults() reports the number of methods, not of problems", {
  R <- PROB_BASEVAR$new(Id = 1, Name = "R", DistributionType = "norm", Mean = 10, Sd = 1)
  S <- PROB_BASEVAR$new(Id = 2, Name = "S", DistributionType = "norm", Mean = 5, Sd = 2)
  lsf <- SYS_LSF$new(name = "g", func = function(R, S) R - S, vars = list(R, S))
  m_form <- PROB_MACHINE$new(name = "FORM", fCall = "FORM")
  m_mvfosm <- PROB_MACHINE$new(name = "MVFOSM", fCall = "MVFOSM")
  ps <- SYS_PROB$new(sys_input = list(lsf), sys_type = "serial",
                     probMachines = list(m_form, m_mvfosm))
  ps$runMachines()

  out <- capture.output(ps$printResults())

  expect_true(any(grepl("by 2 methods", out)))
})

test_that("PRT-2: SYS_PROB$printResults(path) closes the file even after an error (sink)", {
  R <- PROB_BASEVAR$new(Id = 1, Name = "R", DistributionType = "norm", Mean = 10, Sd = 1)
  S <- PROB_BASEVAR$new(Id = 2, Name = "S", DistributionType = "norm", Mean = 5, Sd = 2)
  lsf <- SYS_LSF$new(name = "g", func = function(R, S) R - S, vars = list(R, S))
  m_form <- PROB_MACHINE$new(name = "FORM", fCall = "FORM")
  ps <- SYS_PROB$new(sys_input = list(lsf), sys_type = "serial", probMachines = list(m_form))
  ps$runMachines()
  ps$res_single[[1]][[1]]$beta <- "not a number"  # provokes an error in sprintf()
  tmp <- tempfile()

  expect_error(ps$printResults(tmp))
  n_sink <- sink.number()
  while (sink.number() > 0) sink()  # clean up so that the test output is not swallowed

  expect_equal(n_sink, 0L)
})

test_that("PAR-1: SYS_PARAM$printResults() lists the runs of all parameter variables", {
  P1 <- PARAM_BASEVAR$new(Id = 1, Name = "R", DistributionType = "norm",
                          ParamType = "Mean", Sd = 1, ParamValues = c(10, 11))
  P2 <- PARAM_BASEVAR$new(Id = 2, Name = "S", DistributionType = "norm",
                          ParamType = "Mean", Sd = 2, ParamValues = c(4, 5, 6))
  lsf <- SYS_LSF$new(name = "g", func = function(R, S) R - S, vars = list(P1, P2))
  m_form <- PROB_MACHINE$new(name = "FORM", fCall = "FORM")
  ps <- SYS_PARAM$new(sys_input = list(lsf), probMachines = list(m_form))
  ps$runMachines()

  expect_equal(length(ps$beta_params), 5)
  expect_no_warning(out <- capture.output(ps$printResults()))
  expect_true(any(grepl("Run 5", out)))
})

test_that("PAR-2: SYS_PARAM$printResults() describes a PARAM_DETVAR variable", {
  P <- PARAM_BASEVAR$new(Id = 1, Name = "R", DistributionType = "norm",
                         ParamType = "Mean", Sd = 1, ParamValues = c(10, 11))
  D <- PARAM_DETVAR$new(Id = 2, Name = "D", ParamValues = c(1, 2))
  lsf <- SYS_LSF$new(name = "g", func = function(R, D) R - D, vars = list(P, D))
  m_form <- PROB_MACHINE$new(name = "FORM", fCall = "FORM")
  ps <- SYS_PARAM$new(sys_input = list(lsf), probMachines = list(m_form))
  ps$runMachines()

  out <- capture.output(ps$printResults())

  expect_true(any(grepl("Deterministic", out)))
})

test_that("PAR-3: SYS_PARAM with MVFOSM machine and options without isExpression runs", {
  P1 <- PARAM_BASEVAR$new(Id = 1, Name = "R", DistributionType = "norm",
                          ParamType = "Mean", Sd = 1, ParamValues = c(10, 11))
  S <- PROB_BASEVAR$new(Id = 2, Name = "S", DistributionType = "norm", Mean = 5, Sd = 2)
  lsf <- SYS_LSF$new(name = "g", func = function(R, S) R - S, vars = list(P1, S))
  m_mvfosm <- PROB_MACHINE$new(name = "MVFOSM", fCall = "MVFOSM",
                               options = list(debug.level = 0))
  ps <- SYS_PARAM$new(sys_input = list(lsf), probMachines = list(m_mvfosm))

  expect_no_error(ps$runMachines())
})

test_that("PAR-4: SYS_PARAM$printResults() reports the number of methods", {
  P1 <- PARAM_BASEVAR$new(Id = 1, Name = "R", DistributionType = "norm",
                          ParamType = "Mean", Sd = 1, ParamValues = c(10, 11))
  S <- PROB_BASEVAR$new(Id = 2, Name = "S", DistributionType = "norm", Mean = 5, Sd = 2)
  lsf <- SYS_LSF$new(name = "g", func = function(R, S) R - S, vars = list(P1, S))
  m_form <- PROB_MACHINE$new(name = "FORM", fCall = "FORM")
  m_mvfosm <- PROB_MACHINE$new(name = "MVFOSM", fCall = "MVFOSM")
  ps <- SYS_PARAM$new(sys_input = list(lsf), probMachines = list(m_form, m_mvfosm))
  ps$runMachines()

  out <- capture.output(ps$printResults())

  expect_true(any(grepl("by 2 methods", out)))
})
