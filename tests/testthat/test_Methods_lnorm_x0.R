library(testthat)
library(microbenchmark)

testthat::local_edition(3)


test_that("Spaethe (1991) - Example 2", {

  Force <- PROB_BASEVAR(
    Name = "Force", DistributionType = "gumbel", Package = "evd", Mean = 18, Sd = 2)

  # yield strength in kN
  fy <- PROB_BASEVAR(
    Name = "fy", DistributionType = "slnorm", Package = "brms", Mean = 265000, Sd = 25000, x0 = 160000)

  # plastic moment of resistance in m³ - considered deterministic in this example
  W_pl <- 2.5e-4
  # length of cantilever beam in m - considered deterministic in this example
  l    <- 2

  lsf_param <- SYS_LSF(vars = list(Force, fy), name = "Spaethe (1992) - Beispiel 2")
  lsf_param$func <- function(Force,fy) {
    return(W_pl*fy-l*Force)
  }
  lsf_param$check()

  machine_form <- PROB_MACHINE(name = "FORM", fCall = "FORM")
  #machine_form <- PROB_MACHINE(name = "SORM", fCall = "SORM")
  #machine_form <- PROB_MACHINE(name = "MC CoV 0.05", fCall = "MC_CRUDE", options = list("n_max" = 1e7, "cov_user" = 0.01, "use_threads" = 1L, "seed" = 1234, backend = "future"))


  ps_param <- SYS_PROB(
    sys_input = list(lsf_param),
    probMachines = list(machine_form)
  )

  ps_param$runMachines()

  pf_expected   <- 1.398e-4
  beta_expected <- -qnorm(pf_expected)
  rel_error <-abs(ps_param$beta_single[1] - beta_expected) / beta_expected
  expect_lt(rel_error, 0.02) # Abweichung 2% erlaubt
})


test_that("lnorm with x0 reduces to standard case when x0 = 0", {
  fy1 <- PROB_BASEVAR(Name="fy", DistributionType="slnorm", Package = "brms", Mean=100, Sd=10, x0=0)
  fy2 <- PROB_BASEVAR(Name="fy", DistributionType="lnorm", Mean=100, Sd=10)

  expect_equal(fy1$DistributionParameters[1], fy2$DistributionParameters[1])
  expect_equal(fy1$DistributionParameters[2], fy2$DistributionParameters[2])
  expect_equal(fy1$DistributionParameters[3], 0)
})
