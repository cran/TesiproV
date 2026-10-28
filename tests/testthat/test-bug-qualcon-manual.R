# Plan fck_EN-206 and manually specified OC curve
# f_ck is stored in DefectLimitLower and must always be set, ConformityParam
# only contains c(n) or c(n, N)
# Run with: devtools::test(filter = "bug-qualcon-manual")

test_that("QCM-1: fck_EN-206 without a finite DefectLimitLower stops with a clear message", {
  th <- c(0, 1)
  pa <- c(1, 1)
  hp <- c(3.65, 0.12, 3, 10)

  expect_error(
    .oc_get_outgoing_pdf(th, pa, "fck_EN-206", numeric(0), "lt", hp),
    "DefectLimitLower"
  )
  expect_error(
    .oc_get_outgoing_pdf(th, pa, "fck_EN-206", c(15), "lt", hp, lower = -Inf),
    "f_ck"
  )
})

test_that("QCM-2: fck_EN-206 rejects an upper defect limit", {
  th <- c(0, 1)
  pa <- c(1, 1)
  hp <- c(3.65, 0.12, 3, 10)

  expect_error(
    .oc_get_outgoing_pdf(th, pa, "fck_EN-206", c(15), "lt", hp, lower = 25, upper = 40),
    "DefectLimitUpper"
  )
})

test_that("QCM-3: manual OC curve without selection (P_a = 1) works without ConformityParam (lt)", {
  th <- c(0, 1)
  pa <- c(1, 1)
  hp <- c(3.65, 0.12, 3, 10)

  r <- .oc_get_outgoing_pdf(th, pa, "fck_EN-206", numeric(0), "lt", hp, lower = 25)

  expect_true(all(is.finite(r$f_o)))
  expect_true(is.finite(r$m_o))
  expect_gt(r$s_o, 0)
})

test_that("QCM-4: manual OC curve without selection does not depend on the chosen f_ck", {
  th <- c(0, 1)
  pa <- c(1, 1)
  hp <- c(30, 2, 3, 10)

  r_a <- .oc_get_outgoing_pdf(th, pa, "fck_EN-206", numeric(0), "st", hp, lower = 25)
  r_b <- .oc_get_outgoing_pdf(th, pa, "fck_EN-206", numeric(0), "st", hp, lower = 20)

  expect_equal(r_a$m_o, r_b$m_o, tolerance = 1e-10)
  expect_equal(r_a$s_o, r_b$s_o, tolerance = 1e-10)
})

test_that("QCM-5: updateBasevar() with a manual OC curve and DefectLimitLower runs", {
  v <- PROB_BASEVAR$new(Name = "fc", DistributionType = "lt",
                        DistributionParameters = c(3.65, 0.12, 3, 10))
  qc <- PROB_QUALCON$new(Name = "manual", ConformityType = "fck_EN-206",
                         DefectLimitLower = 25,
                         theta_grid = c(0, 1), Pa_grid = c(1, 1))

  post <- NULL
  expect_no_error(suppressMessages(suppressWarnings(
    post <- v$updateBasevar(PROB_QUALCON = qc, min_r2 = 0.9)
  )))
  expect_true(is.finite(post$Mean))
  expect_gt(post$Sd, 0)
})

test_that("QCM-6: updateBasevar() with a manual OC curve and without DefectLimitLower stops", {
  v <- PROB_BASEVAR$new(Name = "fc", DistributionType = "lt",
                        DistributionParameters = c(3.65, 0.12, 3, 10))
  qc <- PROB_QUALCON$new(Name = "manual", ConformityType = "fck_EN-206",
                         theta_grid = c(0, 1), Pa_grid = c(1, 1))

  expect_error(
    suppressMessages(suppressWarnings(v$updateBasevar(PROB_QUALCON = qc))),
    "DefectLimitLower"
  )
})

test_that("QCM-7: OC curve of fck_EN-206 uses DefectLimitLower as f_ck", {
  th <- c(0.05, 0.2)
  hp <- c(3.75, 0.105, 3, 10)

  set.seed(1)
  pa_25 <- .oc_compute_pa(th, "fck_EN-206", 5, "lt", hp, lower = 25, N = 5000)
  set.seed(1)
  pa_40 <- .oc_compute_pa(th, "fck_EN-206", 5, "lt", hp, lower = 40, N = 5000)

  expect_false(isTRUE(all.equal(pa_25, pa_40)))
  expect_true(all(pa_25 >= 0 & pa_25 <= 1))
})

test_that("QCM-8: OC curve of fck_EN-206 without n or without f_ck stops", {
  th <- c(0.05, 0.2)
  hp <- c(30, 3)

  expect_error(.oc_compute_pa(th, "fck_EN-206", numeric(0), "norm", hp, lower = 25), "must contain n")
  expect_error(.oc_compute_pa(th, "fck_EN-206", c(5), "norm", hp), "DefectLimitLower")
})
