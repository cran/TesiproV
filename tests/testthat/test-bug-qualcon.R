# Verification of findings on PROB_QUALCON, the OC helpers and updateBasevar()
# Run with: devtools::test(filter = "bug-qualcon")

th <- c(0.001, 0.005, 0.007, 0.01, 0.012, 0.015, 0.02, 0.03, 0.05, 0.07,
        0.09, 0.12, 0.15, 0.2, 0.3, 0.4, 0.5, 0.9)
pa_fake <- (1 - th)^15  # monotonically decreasing acceptance probability

test_that("QC-0 (control): outgoing PDF (lt) is finite for JCSS C30", {
  r <- suppressWarnings(.oc_get_outgoing_pdf(
    th, pa_fake, "fck_EN-206", c(15), "lt", c(3.75, 0.105, 3, 10), lower = 30
  ))
  expect_true(all(is.finite(r$f_o)))
})

test_that("QC-1: outgoing PDF (lt) is also finite for a large CoV (JCSS C15)", {
  r <- suppressWarnings(.oc_get_outgoing_pdf(
    th, pa_fake, "fck_EN-206", c(15), "lt", c(3.40, 0.14, 3, 10), lower = 30
  ))
  expect_true(all(is.finite(r$f_o)))
})

test_that("QC-2: outgoing PDF depends on f_ck (DefectLimitLower)", {
  hp <- c(3.75, 0.105, 3, 10)
  r1 <- suppressWarnings(.oc_get_outgoing_pdf(th, pa_fake, "fck_EN-206", c(15), "lt", hp, lower = 25))
  r2 <- suppressWarnings(.oc_get_outgoing_pdf(th, pa_fake, "fck_EN-206", c(15), "lt", hp, lower = 45))
  expect_false(isTRUE(all.equal(r1$m_o, r2$m_o)))
})

test_that("QC-3: n' and s' are not separately identifiable in the st density", {
  x <- seq(30, 55, length.out = 50)
  hp1 <- c(42, 3 * sqrt(5 / 6), 5, 10)
  hp2 <- c(42, 3 * sqrt(50 / 51), 50, 10)
  expect_equal(dst(x, hp1), dst(x, hp2))  # documents the ambiguity
})

test_that("QC-4: parametric fit to an exact st density converges", {
  X <- seq(30, 55, length.out = 300)
  f_o <- dst(X, c(42, 3, 5, 10))
  dx <- X[2] - X[1]
  m_o <- sum(X * f_o) * dx
  s_o <- sqrt(sum((X - m_o)^2 * f_o) * dx)
  fit <- suppressMessages(.oc_fit_posterior_hyper(
    list(X = X, f_o = f_o, m_o = m_o, s_o = s_o), "st", min_r2 = 0.99
  ))
  expect_false(is.null(fit))
})

test_that("QC-5: updateBasevar() computes a missing OC curve itself", {
  skip_on_cran()
  v <- PROB_BASEVAR$new(Name = "fc", DistributionType = "lt",
                        DistributionParameters = c(3.75, 0.105, 3, 10))
  qc <- PROB_QUALCON$new(Name = "EN206", ConformityType = "fck_EN-206",
                         ConformityParam = 15, NSim = 2000, DefectLimitLower = 30, var = v)
  expect_no_error(suppressMessages(suppressWarnings(
    post <- v$updateBasevar(qc)
  )))
  expect_false(all(is.na(qc$Pa_grid)))
})
