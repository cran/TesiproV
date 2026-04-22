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
  machine_sorm <- PROB_MACHINE(name = "SORM", fCall = "SORM")
  machine_MCIS <- PROB_MACHINE("MCIS", "MC_IS",
                               options = list(
                                 cov_user = 0.025,
                                 n_max = 200000, # statt Millionen
                                 n_batch = 20000,
                                 use_threads = 1L,
                                 backend = "future",
                                 seed = 1234
                               ))
  machine_SubSam <- PROB_MACHINE("MC_SubSam", "MC_SubSam")
  machine_MCCrude <- PROB_MACHINE("MC_Crude", "MC_CRUDE")


  ps_param <- SYS_PROB(
    sys_input = list(lsf_param),
    probMachines = list(machine_form,machine_sorm,machine_MCIS,machine_SubSam,machine_MCCrude)
  )

  ps_param$runMachines()

  pf_expected   <- 1.398e-4
  beta_expected <- -qnorm(pf_expected)

  rel_error <- c(1,2,3,4,5)
  for (i in rel_error) {
    rel_error[i] <-abs(ps_param$beta_single[i] - beta_expected) / beta_expected
    expect_lt(rel_error[i], 0.025) # Abweichung 2.5% erlaubt
  }
})


test_that("lnorm with x0 reduces to standard case when x0 = 0", {
  fy1 <- PROB_BASEVAR(Name="fy", DistributionType="slnorm", Package = "brms", Mean=100, Sd=10, x0=0)
  fy2 <- PROB_BASEVAR(Name="fy", DistributionType="lnorm", Mean=100, Sd=10)

  expect_equal(fy1$DistributionParameters[1], fy2$DistributionParameters[1])
  expect_equal(fy1$DistributionParameters[2], fy2$DistributionParameters[2])
  expect_equal(fy1$DistributionParameters[3], 0)
})
