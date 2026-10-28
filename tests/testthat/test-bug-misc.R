# Verifikation der Vermutungen zu ExpressionToFunction(), PROB_DETVAR und updateBasevar()
# Ausfuehren: devtools::test(filter = "bug-misc")

test_that("EXPR-1: Variablennamen, die Teilstring eines anderen Namens sind (R und Rd)", {
  R <- PROB_BASEVAR$new(Id = 1, Name = "R", DistributionType = "norm", Mean = 10, Sd = 1)
  Rd <- PROB_BASEVAR$new(Id = 2, Name = "Rd", DistributionType = "norm", Mean = 3, Sd = 1)
  lsf <- SYS_LSF$new(name = "g", expr = expression(R - Rd), vars = list(R, Rd))

  expect_no_error(lsf$ExpressionToFunction())
  expect_equal(lsf$func(c(10, 3)), 7)
})

test_that("EXPR-2: Variablenname 'e' zerstoert nicht die Funktion exp()", {
  R <- PROB_BASEVAR$new(Id = 1, Name = "R", DistributionType = "norm", Mean = 10, Sd = 1)
  e <- PROB_BASEVAR$new(Id = 2, Name = "e", DistributionType = "norm", Mean = 3, Sd = 1)
  lsf <- SYS_LSF$new(name = "g", expr = expression(R - 2 * e + exp(0)), vars = list(R, e))

  expect_no_error(lsf$ExpressionToFunction())
  expect_equal(lsf$func(c(10, 3)), 5)
})

test_that("EXPR-3: lange Ausdruecke (mehrzeiliges deparse) werden umgewandelt", {
  a <- PROB_BASEVAR$new(Id = 1, Name = "a", DistributionType = "norm", Mean = 10, Sd = 1)
  b <- PROB_BASEVAR$new(Id = 2, Name = "b", DistributionType = "norm", Mean = 3, Sd = 1)
  lsf <- SYS_LSF$new(
    name = "g",
    expr = expression(a * b + a^2 + b^3 + sqrt(a) + sqrt(b) + log(a) + log(b) + exp(a / 10) + exp(b / 10) - 100),
    vars = list(a, b)
  )
  expected <- 10 * 3 + 100 + 27 + sqrt(10) + sqrt(3) + log(10) + log(3) + exp(1) + exp(0.3) - 100

  expect_no_error(lsf$ExpressionToFunction())
  expect_equal(lsf$func(c(10, 3)), expected)
})

test_that("DET-1: PROB_DETVAR mit negativem Wert hat eine nicht negative Standardabweichung", {
  d <- PROB_DETVAR$new(Name = "d", Value = -5)

  expect_equal(d$Mean, -5)
  expect_gte(d$Sd, 0)
})

test_that("DET-2: PROB_DETVAR mit Wert 0 liefert eine nutzbare Verteilung", {
  d <- PROB_DETVAR$new(Name = "d", Value = 0)
  ld <- d$getlDistr()[[1]]

  expect_no_warning(q <- ld$q(0.5))
  expect_equal(q, 0)
})

test_that("UPD-1: updateBasevar meldet den fehlgeschlagenen Fit mit lesbarem Text", {
  skip_on_cran()
  v <- PROB_BASEVAR$new(Name = "fc", DistributionType = "lt",
                        DistributionParameters = c(3.75, 0.105, 3, 10))
  qc <- PROB_QUALCON$new(Name = "EN206", ConformityType = "fck_EN-206",
                         ConformityParam = 15, NSim = 2000, DefectLimitLower = 30, var = v)
  msgs <- character()

  withCallingHandlers(
    suppressWarnings(v$updateBasevar(qc, debug.level = 1, min_r2 = 1.1)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  # exactly one message names the reason for the fallback
  expect_equal(sum(grepl("using the empirical distribution", msgs)), 1)
  expect_true(any(grepl("Parametric fit insufficient", msgs)))
  expect_true(all(nzchar(trimws(msgs))))
})

test_that("UPD-2: updateBasevar liefert bei gleichem globalen Seed dasselbe Ergebnis", {
  skip_on_cran()
  v <- PROB_BASEVAR$new(Name = "fc", DistributionType = "lt",
                        DistributionParameters = c(3.75, 0.105, 3, 10))
  qc_a <- PROB_QUALCON$new(Name = "EN206", ConformityType = "fck_EN-206",
                           ConformityParam = 15, NSim = 2000, DefectLimitLower = 30, var = v)
  qc_b <- PROB_QUALCON$new(Name = "EN206", ConformityType = "fck_EN-206",
                           ConformityParam = 15, NSim = 2000, DefectLimitLower = 30, var = v)

  set.seed(1)
  post_a <- suppressMessages(suppressWarnings(v$updateBasevar(qc_a, min_r2 = 1.1)))
  set.seed(1)
  post_b <- suppressMessages(suppressWarnings(v$updateBasevar(qc_b, min_r2 = 1.1)))

  expect_equal(post_a$DistributionParameters, post_b$DistributionParameters)
})

test_that("PLOT-1: PROB_BASEVAR$plot() runs with Name, Id and Description and accepts own labels", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  v <- PROB_BASEVAR(Id = 1, Name = "Abmessung", Description = "Concrete cover",
                    DistributionType = "st", DistributionParameters = c(30, 2, 3, 10))
  expect_no_error(v$plot())

  # without Id and Description
  v2 <- PROB_BASEVAR(Name = "X", DistributionType = "norm", Mean = 0, Sd = 1)
  expect_no_error(v2$plot())

  # own title and axis labels override the automatic ones
  expect_no_error(v$plot(main = "own title", xlab = "own x", ylab = "own y"))

  # prior and posterior after updateBasevar()
  qc <- PROB_QUALCON(Id = 1, ConformityType = "dim_k||n", ConformityParam = c(25, 1),
                     DefectLimitLower = 27, DefectLimitUpper = 33, var = v)
  post <- suppressMessages(v$updateBasevar(qc))
  expect_no_error(post$plot(prior = TRUE))

  # unit is optional and works together with prior
  expect_no_error(v$plot(unit = "mm"))
  expect_no_error(post$plot(prior = TRUE, unit = "mm"))
  expect_no_error(v$plot(unit = NULL))
})
