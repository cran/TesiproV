# Rackwitz (2006), Zuverlässigkeit und Lasten im konstruktiven Ingenieurbau,
# Teil I: Zuverlässigkeitstheoretische Grundlagen, TU Muenchen 1993-2006, S. 286: Beispiele 1-7, 9 und 11.

# Basisvariablen ---------------------------------------------------------
u1 <- PROB_BASEVAR(Id = 1, Name = "u1", DistributionType = "norm", Mean = 0, Sd = 1)
u2 <- PROB_BASEVAR(Id = 2, Name = "u2", DistributionType = "norm", Mean = 0, Sd = 1)

# Grenzzustandsfunktionen ---------------------------------------------------
lsf1 <- SYS_LSF(vars = list(u1, u2))
lsf1$func <- function(u1, u2) {
  (0.1 * (u1 + 0.5 * (u1^3 + u1)) + 12) / (0.18 * u2 + 3) - 3
}
lsf1$check()

lsf2 <- SYS_LSF(vars = list(u1, u2))
lsf2$func <- function(u1, u2) {
  2 * (0.5 * u1 + 0.87 * u2 + 3) + 5 * (0.5 * u1 - 0.87 * u2 - 1)^2 - 0.5
}
lsf2$check()

lsf3 <- SYS_LSF(vars = list(u1, u2))
lsf3$func <- function(u1, u2) {
  2 * ((-u1) + 0.5 * u2 + 4) + exp(5 * (0.67 * u1 - 0.67 * u2) - 1) - 6
}
lsf3$check()

lsf4 <- SYS_LSF(vars = list(u1, u2))
lsf4$func <- function(u1, u2) {
  u2 + 2 * (1 - 0.4 * u1 * (u1 + 0.8)) + exp(5 * (u1 - u2) - 2.5)
}
lsf4$check()

lsf6 <- SYS_LSF(vars = list(u1, u2))
lsf6$func <- function(u1, u2) 3 - u1 * u2
lsf6$check()

# LSF 5: X1, X2 ~ N(10, 5), Literatur beta 2.29 (NLPQL), hier 2.3121
p1 <- PROB_BASEVAR(Id = 11, Name = "p1", DistributionType = "norm", Mean = 10, Sd = 5)
p2 <- PROB_BASEVAR(Id = 12, Name = "p2", DistributionType = "norm", Mean = 10, Sd = 5)
lsf5 <- SYS_LSF(vars = list(p1, p2))
lsf5$func <- function(p1, p2) p1^3 + p1^2 * p2 + p2^3 - 18
lsf5$check()

# LSF 7: numerisches Rauschen, X1..X4 ~ LN(120, 12), X5 ~ LN(50, 15), X6 ~ LN(40, 12)
q1 <- PROB_BASEVAR(Id = 21, Name = "q1", DistributionType = "lnorm", Mean = 120, Sd = 12)
q2 <- PROB_BASEVAR(Id = 22, Name = "q2", DistributionType = "lnorm", Mean = 120, Sd = 12)
q3 <- PROB_BASEVAR(Id = 23, Name = "q3", DistributionType = "lnorm", Mean = 120, Sd = 12)
q4 <- PROB_BASEVAR(Id = 24, Name = "q4", DistributionType = "lnorm", Mean = 120, Sd = 12)
q5 <- PROB_BASEVAR(Id = 25, Name = "q5", DistributionType = "lnorm", Mean = 50, Sd = 15)
q6 <- PROB_BASEVAR(Id = 26, Name = "q6", DistributionType = "lnorm", Mean = 40, Sd = 12)
lsf7 <- SYS_LSF(vars = list(q1, q2, q3, q4, q5, q6))
lsf7$func <- function(q1, q2, q3, q4, q5, q6) {
  q1 + 2 * q2 + 2 * q3 + q4 - 5 * q5 - 5 * q6 +
    0.001 * (sin(100 * q1) + sin(100 * q2) + sin(100 * q3) +
               sin(100 * q4) + sin(100 * q5) + sin(100 * q6))
}
lsf7$check()

# LSF 9: Ellipsoid, 10 Variablen U_i ~ N(0, 1), beta = sqrt(22) exakt
w1 <- PROB_BASEVAR(Id = 31, Name = "w1", DistributionType = "norm", Mean = 0, Sd = 1)
w2 <- PROB_BASEVAR(Id = 32, Name = "w2", DistributionType = "norm", Mean = 0, Sd = 1)
w3 <- PROB_BASEVAR(Id = 33, Name = "w3", DistributionType = "norm", Mean = 0, Sd = 1)
w4 <- PROB_BASEVAR(Id = 34, Name = "w4", DistributionType = "norm", Mean = 0, Sd = 1)
w5 <- PROB_BASEVAR(Id = 35, Name = "w5", DistributionType = "norm", Mean = 0, Sd = 1)
w6 <- PROB_BASEVAR(Id = 36, Name = "w6", DistributionType = "norm", Mean = 0, Sd = 1)
w7 <- PROB_BASEVAR(Id = 37, Name = "w7", DistributionType = "norm", Mean = 0, Sd = 1)
w8 <- PROB_BASEVAR(Id = 38, Name = "w8", DistributionType = "norm", Mean = 0, Sd = 1)
w9 <- PROB_BASEVAR(Id = 39, Name = "w9", DistributionType = "norm", Mean = 0, Sd = 1)
w10 <- PROB_BASEVAR(Id = 40, Name = "w10", DistributionType = "norm", Mean = 0, Sd = 1)
lsf9 <- SYS_LSF(vars = list(w1, w2, w3, w4, w5, w6, w7, w8, w9, w10))
lsf9$func <- function(w1, w2, w3, w4, w5, w6, w7, w8, w9, w10) {
  w <- c(w1, w2, w3, w4, w5, w6, w7, w8, w9, w10)
  20 - sum(w^2 / (1 + (1:10) / 10))
}
lsf9$check()

# LSF 11: 26 Variablen, X0 ~ N(0.5, 0.1), Xi ~ N(0.2, 0.1), i = 1..25
r0 <- PROB_BASEVAR(Id = 70, Name = "r0", DistributionType = "norm", Mean = 0.5, Sd = 0.1)
r_list <- lapply(1:25, function(i) {
  PROB_BASEVAR(Id = 70 + i, Name = paste0("r", i), DistributionType = "norm", Mean = 0.2, Sd = 0.1)
})
lsf11 <- SYS_LSF(vars = c(list(r0), r_list))
lsf11$func <- function(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13,
                       r14, r15, r16, r17, r18, r19, r20, r21, r22, r23, r24, r25) {
  x <- c(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13,
         r14, r15, r16, r17, r18, r19, r20, r21, r22, r23, r24, r25)
  r0 - sum(x^2 / (1:25))
}
lsf11$check()

# Referenzwerte
beta1 <- 3.536169
beta2 <- 2.773283
beta3 <- 3.337477
beta4 <- 1.868724
beta5 <- 2.312095
beta6 <- 2.449490
beta7 <- 2.348131   # Literatur 2.35, Rauschen: nur mit auto und tolerance 1e-3
beta9 <- sqrt(22)
beta11 <- 2.918088

# Literaturwerte (Rackwitz, auf zwei Stellen gerundet)
beta_lit <- c(3.54, 2.76, 3.32, 1.87, 2.29, 2.45, 2.35, 4.69, 2.92)   # LSF 1-7, 9, 11


# ------------------------------------------------------------------------
# "auto" (Default): SLSQP, bei Bedarf iHLRF
# ------------------------------------------------------------------------
test_that("RACK: FORM auto loest LSF 1, 2, 3, 4, 6 mit Machine-Defaults (n_optim 20, loctol 1e-4)", {
  form_auto <- PROB_MACHINE(
    name = "FORM auto", fCall = "FORM",
    options = list(n_optim = 20, loctol = 1e-4, optim_type = "auto")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf1, lsf2, lsf3, lsf4, lsf6),
    probMachines = list(form_auto),
    debug.level = 0
  )

  expect_no_warning(ps$runMachines())

  expect_equal(ps$beta_single[1], beta1, tolerance = 1e-4)
  expect_equal(ps$beta_single[2], beta2, tolerance = 1e-4)
  expect_equal(ps$beta_single[3], beta3, tolerance = 1e-4)
  expect_equal(ps$beta_single[4], beta4, tolerance = 1e-4)
  expect_equal(ps$beta_single[5], beta6, tolerance = 1e-4)
})

test_that("RACK: FORM auto nutzt zuerst SLSQP und braucht hier keinen Fallback", {
  form_auto <- PROB_MACHINE(
    name = "FORM auto", fCall = "FORM",
    options = list(n_optim = 20, loctol = 1e-4, optim_type = "auto")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf1, lsf2, lsf3, lsf4, lsf6),
    probMachines = list(form_auto),
    debug.level = 0
  )
  ps$runMachines()

  for (i in 1:5) {
    expect_equal(ps$res_single[[i]][[1]]$solver, "SLSQP")
    expect_false(ps$res_single[[i]][[1]]$fallback)
    expect_true(ps$res_single[[i]][[1]]$converged)
  }
})

test_that("RACK: FORM auto loest LSF 1, 2, 3, 4, 6 mit engen Einstellungen (n_optim 200, loctol 1e-6)", {
  form_auto <- PROB_MACHINE(
    name = "FORM auto", fCall = "FORM",
    options = list(n_optim = 200, loctol = 1e-6, optim_type = "auto")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf1, lsf2, lsf3, lsf4, lsf6),
    probMachines = list(form_auto),
    debug.level = 0
  )

  expect_no_warning(ps$runMachines())

  expect_equal(ps$beta_single[1], beta1, tolerance = 1e-4)
  expect_equal(ps$beta_single[2], beta2, tolerance = 1e-4)
  expect_equal(ps$beta_single[3], beta3, tolerance = 1e-4)
  expect_equal(ps$beta_single[4], beta4, tolerance = 1e-4)
  expect_equal(ps$beta_single[5], beta6, tolerance = 1e-4)
})


# ------------------------------------------------------------------------
# "slsqp" explizit
# ------------------------------------------------------------------------
test_that("RACK: FORM slsqp loest LSF 1, 2, 3, 4, 6 (n_optim 20, loctol 1e-4)", {
  form_slsqp <- PROB_MACHINE(
    name = "FORM slsqp", fCall = "FORM",
    options = list(n_optim = 20, loctol = 1e-4, optim_type = "slsqp")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf1, lsf2, lsf3, lsf4, lsf6),
    probMachines = list(form_slsqp),
    debug.level = 0
  )

  expect_no_warning(ps$runMachines())

  expect_equal(ps$beta_single[1], beta1, tolerance = 1e-4)
  expect_equal(ps$beta_single[2], beta2, tolerance = 1e-4)
  expect_equal(ps$beta_single[3], beta3, tolerance = 1e-4)
  expect_equal(ps$beta_single[4], beta4, tolerance = 1e-4)
  expect_equal(ps$beta_single[5], beta6, tolerance = 1e-4)
})


# ------------------------------------------------------------------------
# iHLRF und rackfies explizit: nur dieses Verfahren, kein stiller Fallback
# ------------------------------------------------------------------------
test_that("RACK: FORM ihlrf loest die schwach gekruemmten LSF 1, 4 und 6", {
  form_ihlrf <- PROB_MACHINE(
    name = "FORM iHLRF", fCall = "FORM",
    options = list(n_optim = 20, loctol = 1e-4, optim_type = "ihlrf")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf1, lsf4, lsf6),
    probMachines = list(form_ihlrf),
    debug.level = 0
  )

  expect_no_warning(ps$runMachines())

  expect_equal(ps$beta_single[1], beta1, tolerance = 1e-4)
  expect_equal(ps$beta_single[2], beta4, tolerance = 1e-4)
  expect_equal(ps$beta_single[3], beta6, tolerance = 1e-4)
})

test_that("RACK: FORM rackfies loest LSF 6", {
  form_rackfies <- PROB_MACHINE(
    name = "FORM rackfies", fCall = "FORM",
    options = list(n_optim = 20, loctol = 1e-4, optim_type = "rackfies")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf6),
    probMachines = list(form_rackfies),
    debug.level = 0
  )

  expect_no_warning(ps$runMachines())

  expect_equal(ps$beta_single[1], beta6, tolerance = 1e-4)
})

test_that("RACK-2: FORM ihlrf konvergiert mit n_optim = 200", {
  form_ihlrf <- PROB_MACHINE(
    name = "FORM iHLRF", fCall = "FORM",
    options = list(n_optim = 200, loctol = 1e-6, optim_type = "ihlrf")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf2),
    probMachines = list(form_ihlrf),
    debug.level = 0
  )

  expect_no_warning(ps$runMachines())

  expect_equal(ps$beta_single[1], beta2, tolerance = 1e-4)
})

test_that("RACK-2: FORM ihlrf warnt bei n_optim = 20 (kein stiller Fallback)", {
  form_ihlrf <- PROB_MACHINE(
    name = "FORM iHLRF", fCall = "FORM",
    options = list(n_optim = 20, loctol = 1e-4, optim_type = "ihlrf")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf2),
    probMachines = list(form_ihlrf),
    debug.level = 0
  )

  expect_warning(ps$runMachines(), "no convergence")
})

test_that("RACK-3: FORM ihlrf warnt auch mit n_optim = 200", {
  form_ihlrf <- PROB_MACHINE(
    name = "FORM iHLRF", fCall = "FORM",
    options = list(n_optim = 200, loctol = 1e-6, optim_type = "ihlrf")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf3),
    probMachines = list(form_ihlrf),
    debug.level = 0
  )

  expect_warning(ps$runMachines(), "different optim_type")
})

test_that("RACK-2: FORM rackfies pendelt und warnt", {
  form_rackfies <- PROB_MACHINE(
    name = "FORM rackfies", fCall = "FORM",
    options = list(n_optim = 50, loctol = 1e-4, optim_type = "rackfies")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf2),
    probMachines = list(form_rackfies),
    debug.level = 0
  )

  expect_warning(ps$runMachines(), "different optim_type")
})

test_that("RACK-3: FORM rackfies pendelt und warnt", {
  form_rackfies <- PROB_MACHINE(
    name = "FORM rackfies", fCall = "FORM",
    options = list(n_optim = 50, loctol = 1e-4, optim_type = "rackfies")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf3),
    probMachines = list(form_rackfies),
    debug.level = 0
  )

  expect_warning(ps$runMachines(), "different optim_type")
})


# ------------------------------------------------------------------------
# Weitere Beispiele (5, 7, 9, 11), hoehere Dimension und Rauschen
# ------------------------------------------------------------------------
test_that("RACK: FORM auto loest die Beispiele 5, 7, 9 und 11", {
  form_auto <- PROB_MACHINE(
    name = "FORM auto", fCall = "FORM",
    options = list(n_optim = 100, loctol = 1e-5, optim_type = "auto")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf5, lsf7, lsf9, lsf11),
    probMachines = list(form_auto),
    debug.level = 0
  )

  expect_no_warning(ps$runMachines())

  expect_equal(ps$beta_single[1], beta5, tolerance = 1e-4)
  expect_equal(ps$beta_single[2], beta7, tolerance = 1e-3)
  expect_equal(ps$beta_single[3], beta9, tolerance = 1e-4)
  expect_equal(ps$beta_single[4], beta11, tolerance = 1e-4)
})

test_that("RACK: FORM auto loest die Beispiele 5, 9 und 11", {
  form_auto <- PROB_MACHINE(
    name = "FORM auto", fCall = "FORM",
    options = list(n_optim = 100, loctol = 1e-5, optim_type = "auto")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf5, lsf9, lsf11),
    probMachines = list(form_auto),
    debug.level = 0
  )

  expect_no_warning(ps$runMachines())

  expect_equal(ps$beta_single[1], beta5, tolerance = 1e-4)
  expect_equal(ps$beta_single[2], beta9, tolerance = 1e-4)
  expect_equal(ps$beta_single[3], beta11, tolerance = 1e-4)
})

test_that("RACK: FORM auto loest Beispiel 7 (numerisches Rauschen)", {
  form_auto <- PROB_MACHINE(
    name = "FORM auto", fCall = "FORM",
    options = list(n_optim = 100, loctol = 1e-5, optim_type = "auto")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf7),
    probMachines = list(form_auto),
    debug.level = 0
  )

  # sin() of large arguments differs in the last digits between platforms,
  # the convergence criterion may therefore fail. Only beta is checked.
  suppressWarnings(ps$runMachines())

  expect_equal(ps$beta_single[1], beta7, tolerance = 1e-3)
})

test_that("RACK: FORM ihlrf loest die Beispiele 5 und 11", {
  form_ihlrf <- PROB_MACHINE(
    name = "FORM iHLRF", fCall = "FORM",
    options = list(n_optim = 100, loctol = 1e-5, optim_type = "ihlrf")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf5, lsf11),
    probMachines = list(form_ihlrf),
    debug.level = 0
  )

  expect_no_warning(ps$runMachines())

  expect_equal(ps$beta_single[1], beta5, tolerance = 1e-4)
  expect_equal(ps$beta_single[2], beta11, tolerance = 1e-4)
})

test_that("RACK: FORM rackfies loest die Beispiel 11", {
  form_rackfies <- PROB_MACHINE(
    name = "FORM rackfies", fCall = "FORM",
    options = list(n_optim = 50, loctol = 1e-5, optim_type = "rackfies")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf11),
    probMachines = list(form_rackfies),
    debug.level = 0
  )

  expect_no_warning(ps$runMachines())

  expect_equal(ps$beta_single[1], beta11, tolerance = 1e-4)
})


# ------------------------------------------------------------------------
# Vergleich mit den Literaturwerten (Rackwitz), Toleranz 2 %
# Die Literaturwerte sind gerundet und stammen teils aus nicht voll
# konvergierten Laeufen (LSF 3: 3.32, LSF 5: 2.29), daher keine engere Toleranz.
# ------------------------------------------------------------------------
test_that("RACK: FORM auto stimmt fuer alle Beispiele mit den Literaturwerten ueberein", {
  form_auto <- PROB_MACHINE(
    name = "FORM auto", fCall = "FORM",
    options = list(n_optim = 100, loctol = 1e-5, optim_type = "auto")
  )
  ps <- SYS_PROB(
    sys_input = list(lsf1, lsf2, lsf3, lsf4, lsf5, lsf6, lsf7, lsf9, lsf11),
    probMachines = list(form_auto),
    debug.level = 0
  )

  # convergence warnings are checked in the tests above (example 7 is
  # platform dependent), here only the beta values are compared
  suppressWarnings(ps$runMachines())

  expect_equal(ps$beta_single[1], beta_lit[1], tolerance = 0.02)
  expect_equal(ps$beta_single[2], beta_lit[2], tolerance = 0.02)
  expect_equal(ps$beta_single[3], beta_lit[3], tolerance = 0.02)
  expect_equal(ps$beta_single[4], beta_lit[4], tolerance = 0.02)
  expect_equal(ps$beta_single[5], beta_lit[5], tolerance = 0.02)
  expect_equal(ps$beta_single[6], beta_lit[6], tolerance = 0.02)
  expect_equal(ps$beta_single[7], beta_lit[7], tolerance = 0.02)
  expect_equal(ps$beta_single[8], beta_lit[8], tolerance = 0.02)
  expect_equal(ps$beta_single[9], beta_lit[9], tolerance = 0.02)
})
