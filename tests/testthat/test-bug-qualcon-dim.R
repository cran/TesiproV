# Verification of the attribute plan "dim_k||n" in the Bayesian update
# The tolerance values are stored in PROB_QUALCON$DefectLimitLower / $DefectLimitUpper
# and can be chosen freely, ConformityParam only contains c(n, c_accept).
# Run with: devtools::test(filter = "bug-qualcon-dim")

test_that("QCD-1: one-sided lower limit matches the Monte-Carlo reference", {
  th <- c(0.001, 0.005, 0.007, 0.01, 0.012, 0.015, 0.02, 0.03, 0.05, 0.07,
          0.09, 0.12, 0.15, 0.2, 0.3, 0.4, 0.5, 0.9)
  hp <- c(30, 2, 3, 10)
  cp <- c(10, 1)
  pa <- .oc_compute_pa(th, "dim_k||n", cp, "st", hp)
  r  <- .oc_get_outgoing_pdf(th, pa, "dim_k||n", cp, "st", hp, lower = 27, upper = Inf)

  # Reference: draw the prior (M, S), weight with P_a(theta(M, S)), X | M, S normal
  set.seed(1)
  N  <- 2e6
  v  <- hp[4]
  S  <- hp[2] * sqrt(v / stats::rchisq(N, v))
  M  <- stats::rnorm(N, hp[1], S / sqrt(hp[3]))
  U  <- (27 - M) / S / sqrt(v / (v - 2))
  w  <- stats::approx(th, pa, xout = stats::pt(U, v), rule = 2, ties = "ordered")$y
  x  <- stats::rnorm(N, M, S)
  m_ref <- sum(w * x) / sum(w)
  s_ref <- sqrt(sum(w * (x - m_ref)^2) / sum(w))

  expect_true(all(is.finite(r$f_o)))
  expect_lt(abs(r$m_o - m_ref), 0.02)
  expect_lt(abs(r$s_o - s_ref), 0.05)
})

test_that("QCD-2: acceptance with few defects raises the mean compared to the prior", {
  th <- c(0.001, 0.005, 0.007, 0.01, 0.012, 0.015, 0.02, 0.03, 0.05, 0.07,
          0.09, 0.12, 0.15, 0.2, 0.3, 0.4, 0.5, 0.9)
  hp <- c(30, 2, 3, 10)
  cp <- c(10, 1)
  pa_plan  <- .oc_compute_pa(th, "dim_k||n", cp, "st", hp)
  pa_all <- rep(1, length(th))
  r_plan  <- .oc_get_outgoing_pdf(th, pa_plan,  "dim_k||n", cp, "st", hp, lower = 27)
  r_all <- .oc_get_outgoing_pdf(th, pa_all, "dim_k||n", cp, "st", hp, lower = 27)

  expect_gt(r_plan$m_o, r_all$m_o)
  expect_lt(r_plan$s_o, r_all$s_o)
})

test_that("QCD-3: stricter plan (c_accept = 0) shifts more than a lenient plan (c_accept = 3)", {
  th <- c(0.001, 0.005, 0.007, 0.01, 0.012, 0.015, 0.02, 0.03, 0.05, 0.07,
          0.09, 0.12, 0.15, 0.2, 0.3, 0.4, 0.5, 0.9)
  hp <- c(30, 2, 3, 10)
  pa_strict <- .oc_compute_pa(th, "dim_k||n", c(10, 0), "st", hp)
  pa_lenient   <- .oc_compute_pa(th, "dim_k||n", c(10, 3), "st", hp)
  r_strict  <- .oc_get_outgoing_pdf(th, pa_strict, "dim_k||n", c(10, 0), "st", hp, lower = 27)
  r_lenient    <- .oc_get_outgoing_pdf(th, pa_lenient,   "dim_k||n", c(10, 3), "st", hp, lower = 27)

  expect_gt(r_strict$m_o, r_lenient$m_o)
})

test_that("QCD-4: without a finite limit the update stops with a clear message", {
  th <- c(0.001, 0.005, 0.007, 0.01, 0.012, 0.015, 0.02, 0.03, 0.05, 0.07,
          0.09, 0.12, 0.15, 0.2, 0.3, 0.4, 0.5, 0.9)
  hp <- c(30, 2, 3, 10)
  pa <- .oc_compute_pa(th, "dim_k||n", c(10, 1), "st", hp)

  expect_error(
    .oc_get_outgoing_pdf(th, pa, "dim_k||n", c(10, 1), "st", hp),
    "finite"
  )
})

test_that("QCD-5: updateBasevar() with an attribute plan returns a posterior with a higher mean", {
  v <- PROB_BASEVAR$new(Name = "cover", DistributionType = "st",
                        DistributionParameters = c(30, 2, 3, 10))
  qc <- PROB_QUALCON$new(Name = "Cover", ConformityType = "dim_k||n",
                         ConformityParam = c(10, 1), DefectLimitLower = 27, var = v)
  qc$getOC()

  post <- NULL
  expect_no_error(suppressMessages(suppressWarnings(
    post <- v$updateBasevar(qc)
  )))

  expect_false(all(is.na(qc$Pa_grid)))
  expect_gt(post$Mean, v$Mean)
  expect_lt(post$Sd, v$Sd)
  expect_true(post$DistributionType %in% c("st", "emp"))
})

test_that("QCD-6: two-sided plan matches the Monte-Carlo reference", {
  th <- c(0.001, 0.005, 0.007, 0.01, 0.012, 0.015, 0.02, 0.03, 0.05, 0.07,
          0.09, 0.12, 0.15, 0.2, 0.3, 0.4, 0.5, 0.9)
  hp <- c(30, 2, 3, 10)
  cp <- c(10, 1)
  pa <- .oc_compute_pa(th, "dim_k||n", cp, "st", hp)
  r  <- .oc_get_outgoing_pdf(th, pa, "dim_k||n", cp, "st", hp, lower = 27, upper = 33)

  # theta is the share below the lower limit plus the share above the upper limit
  set.seed(2)
  N  <- 2e6
  v  <- hp[4]
  S  <- hp[2] * sqrt(v / stats::rchisq(N, v))
  M  <- stats::rnorm(N, hp[1], S / sqrt(hp[3]))
  k  <- sqrt(v / (v - 2))
  tt <- stats::pt((27 - M) / S / k, v) + 1 - stats::pt((33 - M) / S / k, v)
  w  <- stats::approx(th, pa, xout = tt, rule = 2, ties = "ordered")$y
  x  <- stats::rnorm(N, M, S)
  m_ref <- sum(w * x) / sum(w)
  s_ref <- sqrt(sum(w * (x - m_ref)^2) / sum(w))

  expect_true(all(is.finite(r$f_o)))
  expect_lt(abs(r$m_o - m_ref), 0.02)
  expect_lt(abs(r$s_o - s_ref), 0.05)
})

test_that("QCD-7: an upper limit alone lowers the mean, a lower limit alone raises it", {
  th <- c(0.001, 0.005, 0.007, 0.01, 0.012, 0.015, 0.02, 0.03, 0.05, 0.07,
          0.09, 0.12, 0.15, 0.2, 0.3, 0.4, 0.5, 0.9)
  hp <- c(30, 2, 3, 10)
  pa <- .oc_compute_pa(th, "dim_k||n", c(10, 1), "st", hp)
  r_lo <- .oc_get_outgoing_pdf(th, pa, "dim_k||n", c(10, 1), "st", hp, lower = 27, upper = Inf)
  r_hi  <- .oc_get_outgoing_pdf(th, pa, "dim_k||n", c(10, 1), "st", hp, lower = -Inf, upper = 33)

  expect_gt(r_lo$m_o, 30)
  expect_lt(r_hi$m_o, 30)
})

test_that("QCD-8: two-sided plan reduces the spread more than either side alone", {
  th <- c(0.001, 0.005, 0.007, 0.01, 0.012, 0.015, 0.02, 0.03, 0.05, 0.07,
          0.09, 0.12, 0.15, 0.2, 0.3, 0.4, 0.5, 0.9)
  hp <- c(30, 2, 3, 10)
  pa <- .oc_compute_pa(th, "dim_k||n", c(10, 1), "st", hp)
  r_both <- .oc_get_outgoing_pdf(th, pa, "dim_k||n", c(10, 1), "st", hp, lower = 27, upper = 33)
  r_lo <- .oc_get_outgoing_pdf(th, pa, "dim_k||n", c(10, 1), "st", hp, lower = 27, upper = Inf)
  r_hi  <- .oc_get_outgoing_pdf(th, pa, "dim_k||n", c(10, 1), "st", hp, lower = -Inf, upper = 33)

  expect_lt(r_both$s_o, r_lo$s_o)
  expect_lt(r_both$s_o, r_hi$s_o)
})

test_that("QCD-9: invalid limits are rejected", {
  th <- c(0.001, 0.005, 0.007, 0.01, 0.012, 0.015, 0.02, 0.03, 0.05, 0.07,
          0.09, 0.12, 0.15, 0.2, 0.3, 0.4, 0.5, 0.9)
  hp <- c(30, 2, 3, 10)
  pa <- .oc_compute_pa(th, "dim_k||n", c(10, 1), "st", hp)

  expect_error(
    .oc_get_outgoing_pdf(th, pa, "dim_k||n", c(10, 1), "st", hp, lower = 33, upper = 27),
    "DefectLimitLower < DefectLimitUpper"
  )
  expect_error(
    .oc_get_outgoing_pdf(th, pa, "dim_k||n", c(10, 1), "st", hp, lower = -Inf, upper = Inf),
    "finite"
  )
})

test_that("QCD-10: PROB_QUALCON stores limits with defaults -Inf / Inf and validates them", {
  qc <- PROB_QUALCON$new(Name = "Default", ConformityType = "dim_k||n",
                         ConformityParam = c(10, 1))
  expect_equal(qc$DefectLimitLower, -Inf)
  expect_equal(qc$DefectLimitUpper, Inf)

  qc2 <- PROB_QUALCON$new(Name = "Limits", ConformityType = "dim_k||n",
                          ConformityParam = c(10, 1), DefectLimitLower = 27, DefectLimitUpper = 33)
  expect_equal(qc2$DefectLimitLower, 27)
  expect_equal(qc2$DefectLimitUpper, 33)

  expect_error(
    PROB_QUALCON$new(Name = "invalid", ConformityType = "dim_k||n",
                     ConformityParam = c(10, 1), DefectLimitLower = 33, DefectLimitUpper = 27),
    "smaller"
  )
})
