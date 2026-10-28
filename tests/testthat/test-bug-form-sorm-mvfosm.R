# Verifikation der Befunde zu FORM, SORM und MVFOSM
# Ausfuehren: devtools::test(filter = "bug-form")

# Linearer Grenzzustand g = R - S mit R ~ N(10, 1), S ~ N(5, 2): beta = 5 / sqrt(5)
R <- PROB_BASEVAR$new(Name = "R", DistributionType = "norm", Mean = 10, Sd = 1)
S <- PROB_BASEVAR$new(Name = "S", DistributionType = "norm", Mean = 5, Sd = 2)
lDistr_RS <- list(R$getlDistr(), S$getlDistr())
lsf_lin <- function(x) x[1] - x[2]

test_that("FORM-0 (Kontrolle): Beta stimmt bei linearem Normalfall", {
  res_ihlrf <- FORM(lsf_lin, lDistr_RS, optim_type = "ihlrf")
  res_rackfies <- FORM(lsf_lin, lDistr_RS, optim_type = "rackfies")
  res_auto <- FORM(lsf_lin, lDistr_RS)
  res_slsqp <- FORM(lsf_lin, lDistr_RS, optim_type = "slsqp")

  expect_equal(res_ihlrf$beta, 5 / sqrt(5), tolerance = 1e-3)
  expect_equal(res_rackfies$beta, 5 / sqrt(5), tolerance = 1e-3)
  expect_equal(res_auto$beta, 5 / sqrt(5), tolerance = 1e-3)
  expect_equal(res_slsqp$beta, 5 / sqrt(5), tolerance = 1e-3)
})

test_that("FORM-1: grosses Beta (> 7) wird erreicht", {
  R2 <- PROB_BASEVAR$new(Name = "R", DistributionType = "norm", Mean = 20, Sd = 1)
  S2 <- PROB_BASEVAR$new(Name = "S", DistributionType = "norm", Mean = 5, Sd = 1)
  lDistr_big <- list(R2$getlDistr(), S2$getlDistr())
  expected <- 15 / sqrt(2)  # 10.607

  res_ihlrf <- FORM(lsf_lin, lDistr_big, optim_type = "ihlrf")
  res_rackfies <- FORM(lsf_lin, lDistr_big, optim_type = "rackfies")

  expect_equal(res_ihlrf$beta, expected, tolerance = 1e-2)
  expect_equal(res_rackfies$beta, expected, tolerance = 1e-2)
})

test_that("FORM-2: FORM liefert u_points und ein sichtbares Ergebnis", {
  vis <- withVisible(FORM(lsf_lin, lDistr_RS, optim_type = "ihlrf"))
  expect_true(vis$visible)
  expect_false(is.null(vis$value$u_points))
  expect_true(vis$value$converged)
})

test_that("FORM-3: fehlende Konvergenz wird gemeldet (Warnung und Flag)", {
  expect_warning(
    res <- FORM(lsf_lin, lDistr_RS, n_optim = 0L, loctol = 1e-12, optim_type = "ihlrf"),
    "no convergence"
  )
  expect_false(res$converged)
})

test_that("SORM-1: Kruemmung wird beruecksichtigt (pf_SORM < pf_FORM, nahe am exakten Wert)", {
  # Versagensbereich {x2 > 3 + 0.2 x1^2} ist Teilmenge des FORM-Halbraums
  A <- PROB_BASEVAR$new(Name = "A", DistributionType = "norm", Mean = 0, Sd = 1)
  B <- PROB_BASEVAR$new(Name = "B", DistributionType = "norm", Mean = 0, Sd = 1)
  lDistr_AB <- list(A$getlDistr(), B$getlDistr())
  lsf <- function(x) 3 - x[2] + 0.2 * x[1]^2

  res <- suppressWarnings(SORM(lsf, lDistr_AB, debug.level = 0))
  exact <- integrate(function(u) dnorm(u) * pnorm(-(3 + 0.2 * u^2)), -Inf, Inf)$value

  expect_lt(res$pf, res$FORM_pf)
  expect_lt(abs(res$pf / exact - 1), 0.2)
})

test_that("SORM-2: Ergebnis ist invariant gegen lineare Skalierung der Variablen", {
  g_u <- function(u) 3 - u[2] + 0.2 * u[1]^2 + 0.05 * u[1]^3

  A <- PROB_BASEVAR$new(Name = "A", DistributionType = "norm", Mean = 0, Sd = 1)
  B <- PROB_BASEVAR$new(Name = "B", DistributionType = "norm", Mean = 0, Sd = 1)
  std <- suppressWarnings(SORM(g_u, list(A$getlDistr(), B$getlDistr())))

  A2 <- PROB_BASEVAR$new(Name = "A", DistributionType = "norm", Mean = 10, Sd = 2)
  B2 <- PROB_BASEVAR$new(Name = "B", DistributionType = "norm", Mean = 5, Sd = 1)
  lsf_x <- function(x) g_u(c((x[1] - 10) / 2, (x[2] - 5) / 1))
  scl <- suppressWarnings(SORM(lsf_x, list(A2$getlDistr(), B2$getlDistr())))

  expect_lt(scl$pf, scl$FORM_pf)
  expect_lt(abs(scl$pf / std$pf - 1), 0.05)
})

test_that("MVFOSM-1: Designpunkt hat das richtige Vorzeichen", {
  ld_R <- R$getlDistr()[[1]]
  ld_S <- S$getlDistr()[[1]]
  d <- list(list(ld_R, c(ld_R$mean, ld_R$Sd)),
            list(ld_S, c(ld_S$mean, ld_S$Sd)))

  res <- MVFOSM(lsf_lin, d, debug.level = 0)

  expect_equal(res$beta, 5 / sqrt(5), tolerance = 1e-6)
  expect_equal(unname(res$design.point.x), c(9, 9), tolerance = 1e-3)
  expect_equal(lsf_lin(res$design.point.x), 0, tolerance = 1e-6)
})

test_that("MVFOSM-2: debug.level hat einen Default", {
  ld_R <- R$getlDistr()[[1]]
  ld_S <- S$getlDistr()[[1]]
  d <- list(list(ld_R, c(ld_R$mean, ld_R$Sd)),
            list(ld_S, c(ld_S$mean, ld_S$Sd)))

  expect_no_error(MVFOSM(lsf_lin, d))
})

test_that("MVFOSM-3: isExpression = TRUE funktioniert bei nichtlinearem g", {
  Rres <- PROB_BASEVAR$new(Name = "Rres", DistributionType = "norm", Mean = 10, Sd = 1)
  Sact <- PROB_BASEVAR$new(Name = "Sact", DistributionType = "norm", Mean = 5, Sd = 2)
  ld_R <- Rres$getlDistr()[[1]]
  ld_S <- Sact$getlDistr()[[1]]
  d <- list(list(ld_R, c(ld_R$mean, ld_R$Sd)),
            list(ld_S, c(ld_S$mean, ld_S$Sd)))

  b_fun <- MVFOSM(function(x) x[1] * x[2] - 20, d, debug.level = 0)$beta
  b_expr <- MVFOSM(expression(Rres * Sact - 20), d, isExpression = TRUE, debug.level = 0)$beta

  expect_equal(b_expr, b_fun, tolerance = 1e-4)
})
