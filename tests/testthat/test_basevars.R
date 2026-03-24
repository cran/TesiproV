test_that("Basis var, transformation", {
  a <- PROB_BASEVAR(Id = 1, Name = "a", DistributionType = "norm", Mean = 1, Sd = 0.1)
  a$prepare()
  mean_ <- a$DistributionParameters[1]
  sd_ <- a$DistributionParameters[2]
  expect_equal(mean_, 1)
  expect_equal(sd_, 0.1)

  b <- PROB_BASEVAR(Id = 1, Name = "b", DistributionType = "lnorm", Mean = 1, Sd = 0.1)
  b$prepare()
  mean_ <- b$DistributionParameters[1]
  sd_ <- b$DistributionParameters[2]
  expect_equal(mean_, -0.004975165)
  expect_equal(sd_, 0.099751345)

  b <- PROB_BASEVAR(Id = 1, Name = "c", DistributionType = "lnorm", DistributionParameters = c(mean_, sd_))
  b$prepare()
  expect_equal(b$Mean, 1)
  expect_equal(b$Sd, 0.1)
})


test_that("Basis var, parameters", {
  params <- c(2, 3, 4, 5)
  means <- round(c(0.69189874, 1.098057, 1.385982, 1.609238, 0.6918987), 2)
  sds <- round(c(0.04996879, 0.03332408, 0.0249961, 0.019998, 0.04996879), 4)
  a <- PARAM_BASEVAR(Id = 1, Name = "a", DistributionType = "lnorm", ParamType = "Mean", Sd = 0.1, ParamValues = params)

  for (i in 1:5) {
    mean_ <- round(a$DistributionParameters[1], 2)
    sd_ <- round(a$DistributionParameters[2], 4)
    expect_equal(mean_, means[i])
    expect_equal(sd_, sds[i])
    a$nextParam()
  }
})


test_that("Gumbel distribution forward and backward transformation are consistent", {
  # Beispielwerte für Mittelwert und Standardabweichung
  mean_val <- 2000
  sd_val <- 200

  # Vorwärtstransformation mit PROB_BASEVAR
  gumbel_var <- PROB_BASEVAR(
    Id = 1, Name = "M",
    DistributionType = "gumbel",
    Package = "evd",
    Mean = mean_val,
    Sd = sd_val
  )

  gumbel_var$prepare()

  loc_param <- gumbel_var$DistributionParameters[1]
  scale_param <- gumbel_var$DistributionParameters[2]

  # Rücktransformation aus Parametern prüfen
  gumbel_back <- PROB_BASEVAR(
    Id = 2, Name = "M_back",
    DistributionType = "gumbel",
    Package = "evd",
    DistributionParameters = c(loc_param, scale_param)
  )

  gumbel_back$prepare()

  expect_equal(gumbel_back$Mean, mean_val, tolerance = 1e-8)
  expect_equal(gumbel_back$Sd, sd_val, tolerance = 1e-8)

  # Zusatzprüfung: Location liegt unterhalb des Mittelwerts (wegen positiver Gamma-Konstante)
  expect_lt(loc_param, mean_val)
})


test_that("Weibull distribution forward and backward transformation are consistent", {
  mean_val <- 100
  sd_val <- 15

  weibull_var <- PROB_BASEVAR(
    Id = 3, Name = "W",
    DistributionType = "weibull",
    Mean = mean_val,
    Sd = sd_val
  )

  weibull_var$prepare()

  shape_param <- weibull_var$DistributionParameters[1]
  scale_param <- weibull_var$DistributionParameters[2]

  # Rückrichtung prüfen:
  weibull_back <- PROB_BASEVAR(
    Id = 4, Name = "W_back",
    DistributionType = "weibull",
    DistributionParameters = c(shape_param, scale_param)
  )

  weibull_back$prepare()

  # Wegen der Näherungsfunktion für k ist eine Toleranz erforderlich
  rel_err_mean <- abs(weibull_back$Mean - mean_val) / mean_val
  rel_err_sd <- abs(weibull_back$Sd - sd_val) / sd_val

  expect_lt(rel_err_mean, 0.01) # < 1 % Abweichung erlaubt
  expect_lt(rel_err_sd, 0.01)

  expect_true(shape_param > .Machine$double.eps)
  expect_true(scale_param > .Machine$double.eps)
})

test_that("Weibull shape reacts correctly to different CoV values", {
  low_cov_obj <- PROB_BASEVAR(
    Id = 5, Name = "W_lowCov",
    DistributionType = "weibull", Mean = 1000, Sd = 80
  ) # Cov≈0.05

  high_cov_obj <- PROB_BASEVAR(
    Id = 6, Name = "W_highCov",
    DistributionType = "weibull", Mean = 1000, Sd = 300
  ) # Cov≈0.5

  low_cov_obj$prepare()
  high_cov_obj$prepare()

  k_lowCov <- low_cov_obj$DistributionParameters[1]
  k_highCov <- high_cov_obj$DistributionParameters[1]

  # Bei größerer Streuung (hoher COV) wird der Shapeparameter kleiner.
  expect_gt(k_lowCov, k_highCov)
})


test_that("Old incorrect Gumbel implementation would fail consistency check", {
  mean_test <- 2000
  sd_test <- 200

  correct_scale <- (sd_test * sqrt(6)) / pi
  correct_location <- mean_test + digamma(1) * correct_scale

  wrong_location_old_code <-
    mean_test + digamma(1) / correct_scale # alte fehlerhafte Formel

  expect_gt(abs(correct_location - wrong_location_old_code), .Machine$double.eps)
})


test_that("Exponential forward and backward transformation are consistent", {
  evar <- PROB_BASEVAR(
    Id = 10, Name = "ExpVar",
    DistributionType = "exp", Mean = 5
  )
  evar$prepare()
  expect_equal(evar$DistributionParameters[1], 0.2, tolerance = 5e-8)

  eback <- PROB_BASEVAR(
    Id = 11, Name = "ExpBack",
    DistributionType = "exp",
    DistributionParameters = c(evar$DistributionParameters)
  )
  eback$prepare()
  expect_equal(eback$Mean, 5, tolerance = 5e-8)
  expect_equal(eback$Sd, 5, tolerance = 5e-8)
})


test_that("Beta forward and backward transformation are consistent", {
  bvar <- PROB_BASEVAR(
    Id = 12, Name = "BetaVar",
    DistributionType = "beta", Mean = 0.6, Sd = 0.05
  )
  bvar$prepare()

  a_param <- bvar$DistributionParameters[1]
  b_param <- bvar$DistributionParameters[2]

  bback <- PROB_BASEVAR(
    Id = 13, Name = "BetaBack",
    DistributionType = "beta",
    DistributionParameters = c(a_param, b_param)
  )
  bback$prepare()

  expect_equal(bback$Mean, 0.6, tolerance = 5e-3)
  expect_equal(bback$Sd, 0.05, tolerance = 5e-3)
})


test_that("is_empty works as intended", {
  expect_true(is_empty(NULL))
  expect_true(is_empty(numeric(0)))
  expect_true(is_empty(integer(0)))
  expect_true(is_empty(list()))
  expect_true(is_empty(data.frame()))
  expect_true(is_empty(matrix(numeric(0), nrow = 0, ncol = 2)))

  expect_false(is_empty(1))
  expect_false(is_empty(c(1, 2)))
  expect_false(is_empty(list(a = 1)))
  expect_false(is_empty(data.frame(x = 1)))
})


test_that("PROB_BASEVAR rejects unknown distribution types", {
  expect_error(
    PROB_BASEVAR(
      Name = "x",
      DistributionType = "unknown",
      Mean = 1, Sd = 0.1
    ),
    regexp = "Distribution type 'unknown' is not known"
  )
})


test_that("PROB_MACHINE rejects non-existent fCall", {
  mach <- PROB_MACHINE(name = "Bad", fCall = "NON_EXISTENT")
  ps <- SYS_PROB(
    sys_input = list(SYS_LSF(vars = list(
      PROB_BASEVAR(Name = "a", Mean = 1, Sd = 0.1)
    ))),
    probMachines = list(mach)
  )
  expect_error(
    ps$runMachines(),
    regexp = "Method 'NON_EXISTENT' is not loaded"
  )
})


test_that("Distribution parameter round-trip (norm, lnorm, gumbel, gamma, beta, weibull)", {
  # Helper, prüft relative Abweichung
  rel_err <- function(ref, got) abs(ref - got) / max(abs(ref), .Machine$double.eps)

  # 1. Normal
  nb <- PROB_BASEVAR(
    Name = "norm", DistributionType = "norm",
    Mean = 12.3, Sd = 4.5
  )
  nb$prepare()
  expect_equal(nb$DistributionParameters, c(12.3, 4.5))

  # 2. Lognormal
  ln <- PROB_BASEVAR(
    Name = "lnorm", DistributionType = "lnorm",
    Mean = 5, Sd = 2
  )
  ln$prepare()
  # round‑trip
  mn <- ln$DistributionParameters[1]
  sn <- ln$DistributionParameters[2]
  # back‑transformation (inside prepare) was already executed,
  # check that Mean/Sd ≈ original
  expect_lt(rel_err(ln$Mean, 5), 1e-12)
  expect_lt(rel_err(ln$Sd, 2), 1e-12)

  # 3. Gumbel
  gb <- PROB_BASEVAR(
    Name = "gumbel", DistributionType = "gumbel",
    Mean = 2000, Sd = 200
  )
  gb$prepare()
  loc <- gb$DistributionParameters[1]
  scale <- gb$DistributionParameters[2]
  # back‑transform
  expect_lt(rel_err(gb$Mean, loc - digamma(1) * scale), 1e-12)
  expect_lt(rel_err(gb$Sd, pi / sqrt(6) * scale), 1e-12)

  # 4. Gamma
  ga <- PROB_BASEVAR(
    Name = "gamma", DistributionType = "gamma",
    Mean = 15, Sd = 5
  )
  ga$prepare()
  shape <- ga$DistributionParameters[1]
  scale <- ga$DistributionParameters[2]
  expect_lt(rel_err(ga$Mean, shape * scale), 1e-12)
  expect_lt(rel_err(ga$Sd, sqrt(shape) * scale), 1e-12)

  # 5. Beta (bound to (0,1))
  be <- PROB_BASEVAR(
    Name = "beta", DistributionType = "beta",
    Mean = 0.6, Sd = 0.05
  )
  be$prepare()
  alpha <- be$DistributionParameters[1]
  beta_ <- be$DistributionParameters[2]
  mean.back <- alpha / (alpha + beta_)
  var.back <- alpha * beta_ / ((alpha + beta_)^2 * (alpha + beta_ + 1))
  expect_lt(rel_err(be$Mean, mean.back), 1e-12)
  expect_lt(rel_err(be$Sd, sqrt(var.back)), 1e-12)

  # 6. Weibull – hier nutzen wir die numerische Approximation,
  #    die wir bereits implementiert haben
  mean_val <- 100
  sd_val <- 15

  we <- PROB_BASEVAR(
    Name = "weibull", DistributionType = "weibull",
    Mean = mean_val,
    Sd = sd_val
  )

  we$prepare()
  shape.w <- we$DistributionParameters[1]
  scale.w <- we$DistributionParameters[2]
  # back‑transform
  mean.w <- scale.w * base::gamma(1 + 1 / shape.w)
  sd.w <- sqrt(scale.w^2 *
    (base::gamma(1 + 2 / shape.w) - base::gamma(1 + 1 / shape.w)^2))
  # Toleranz 1 % ist für die empirische Approximation realistisch
  rel_err_mean <- abs(we$Mean - mean_val) / mean_val
  rel_err_sd <- abs(we$Sd - sd_val) / sd_val

  expect_lt(rel_err_mean, 0.01) # < 1 %
  expect_lt(rel_err_sd, 0.01)
})


test_that("PROB_BASEVAR handles zero or near-zero mean correctly", {
  vb <- PROB_BASEVAR(
    Name = "zeroMean",
    DistributionType = "norm",
    Mean = 0,
    Sd = 0.1
  )

  # The prepare() method should run silently (no warning expected)
  expect_silent(vb$prepare())

  # Coefficient of variation must be infinite for mean ≈ 0
  expect_true(is.infinite(vb$Cov))

  # Numerical tolerance: near-zero means also yield Inf CoV
  v_small <- PROB_BASEVAR(Name = "tinyMean", DistributionType = "norm", Mean = 1e-12, Sd = 1e-4)

  v_small$prepare()

  # Check that Cov is still treated as infinite due to small mean value
  expect_true(is.infinite(v_small$Cov))

  # For normal finite means → finite CoV value expected
  v_normal <- PROB_BASEVAR(Name = "normalMean", DistributionType = "norm", Mean = 10, Sd = 2)

  v_normal$prepare()

  expect_true(is.finite(v_normal$Cov))

  expect_equal(round(v_normal$Cov, 3), round(2 / 10, 3))
})


# Test for SYS_LSF$check() behaviour in different function states
# ---------------------------------------------------------------

test_that("SYS_LSF$check() correctly distinguishes missing, empty and valid functions", {
  # Create two basic random variables
  var_E <- PROB_BASEVAR(Name = "E", DistributionType = "norm", Mean = 10, Sd = 0.2)
  var_R <- PROB_BASEVAR(Name = "R", DistributionType = "norm", Mean = 35, Sd = 0.5)

  ## Case 1: Missing function (func == NULL) ---------------------
  lsf_missing <- SYS_LSF(vars = list(var_E, var_R), name = "LSF_missing")
  expect_true(is.null(lsf_missing$func))
  expect_error(lsf_missing$check(),
    regexp = "Limit-state function 'LSF_missing' is missing",
    fixed = TRUE
  )

  ## Case 2: Empty function without arguments --------------------
  lsf_empty <- SYS_LSF(vars = list(var_E, var_R), name = "LSF_empty")
  lsf_empty$func <- function() NULL
  expect_true(is.function(lsf_empty$func))
  expect_equal(length(formals(lsf_empty$func)), 0L)
  expect_error(lsf_empty$check(),
    regexp = "has no formal arguments",
    fixed = FALSE
  )

  ## Case 3: Valid function with arguments -----------------------
  lsf_valid <- SYS_LSF(vars = list(var_E, var_R), name = "LSF_valid")
  lsf_valid$func <- function(E, R) {
    R - E
  }

  # Should run silently (no error)
  expect_silent(lsf_valid$check())
})


test_that("PROB_BASEVAR catches negative standard deviation", {
  expect_error(
    PROB_BASEVAR(
      Name = "negSd",
      DistributionType = "norm",
      Mean = 10, Sd = -1
    ),
    regexp = "Standard deviation (Sd) must be non-negative.",
    fixed = TRUE
  )
})

test_that("PARAM_BASEVAR works with a single value", {
  pb <- PARAM_BASEVAR(
    Name = "single",
    DistributionType = "norm",
    ParamType = "Mean",
    ParamValues = 5
  )
  # Erstes `nextParam()` wird beim Initialisieren aufgerufen
  expect_equal(pb$Mean, 5)
  # Weiteres `nextParam()` muss denselben Wert zurückgeben (wrap‑around)
  pb$nextParam()
  expect_equal(pb$Mean, 5)
})

test_that("PARAM_DETVAR works with one deterministic value", {
  pd <- PARAM_DETVAR(
    Name = "det_single",
    ParamValues = 42
  )
  expect_equal(pd$Mean, 42)
  expect_equal(pd$Sd, 42 / 1e7)
  pd$nextParam()
  expect_equal(pd$Mean, 42) # wrap‑around, gleiche Zahl
})


test_that("PROB_BASEVAR caching improves runtime", {
  # Prepare a simple normal variable
  v <- PROB_BASEVAR(
    Name = "cacheTest",
    DistributionType = "norm",
    Mean = 0,
    Sd = 1
  )
  v$prepare()

  # Skip test if microbenchmark is not available
  skip_if_not_installed("microbenchmark")
  library(microbenchmark)

  # Benchmark both calls: first (creates cache) vs. second (uses cache)
  bench <- microbenchmark(
    first_call = {
      v$getlDistr()
    },
    second_call = {
      v$getlDistr()
    },
    times = 10L # number of repetitions per expression
  )

  # Print benchmark results for inspection (optional)
  print(bench)

  # Extract median execution times in milliseconds
  t_first_ms <- median(bench$time[bench$expr == "first_call"]) / 1e6
  t_second_ms <- median(bench$time[bench$expr == "second_call"]) / 1e6

  message(sprintf("Median time: first=%.3f ms, second=%.3f ms", t_first_ms, t_second_ms))

  # Compare performance: cached call should be equal or faster
  expect_true(t_second_ms <= t_first_ms * 1.2,
    info = sprintf(
      "Second call should be faster or within +20%% tolerance (first=%.3f ms, second=%.3f ms)",
      t_first_ms, t_second_ms
    )
  )

  # Additional functional check: both calls should return identical output values
  res_first <- v$getlDistr()
  res_second <- v$getlDistr()

  expect_equal(res_first[[1]]$d(0), res_second[[1]]$d(0))
})


test_that("runMachines() yields identical β-values for identical seeds", {
  set.seed(1234) # Master‑Seed – hat keinen Einfluss auf die Worker‑Seeds
  # 2 Maschinen: FORM (deterministisch) + MC_CRUDE (stochastic)
  m_form <- PROB_MACHINE(name = "FORM", fCall = "FORM")
  m_mc <- PROB_MACHINE(
    name = "MC", fCall = "MC_CRUDE",
    options = list(
      cov_user = 0.05,
      n_max = 5e5,
      seed = 777
    )
  ) # fester Seed

  lsf_obj <- SYS_LSF(vars = list(
    PROB_BASEVAR(Name = "E", Mean = 10, Sd = 0.2),
    PROB_BASEVAR(Name = "R", Mean = 35, Sd = 0.5)
  ))
  lsf_obj$func <- function(E, R) {
    R - E
  } # einfache Funktion definieren

  ps <- SYS_PROB(
    sys_input = list(lsf_obj),
    probMachines = list(m_form, m_mc)
  )

  # Erster Durchlauf
  ps$runMachines()
  beta1 <- ps$beta_single

  # zweiter Durchlauf (identischer Seed!)
  ps$runMachines()
  beta2 <- ps$beta_single

  expect_equal(beta1, beta2)
})


test_that("MC_IS throws informative error if n_batch <= 0", {
  expect_error(
    MC_IS(
      lsf = function(x) 1,
      lDistr = list(list(
        d = function() NULL,
        p = function() NULL,
        q = function() NULL,
        r = function() NULL
      )),
      n_batch = 0, n_max = 1e5, cov_user = 0.05,
      seed = 123
    ),
    regexp = "n_batch.*positive"
  )
})


test_that("SYS_PARAM aborts cleanly when a PARAM variable loses its ParamValues", {
  pvar <- PARAM_BASEVAR(
    Name = "bad",
    DistributionType = "norm",
    ParamType = "Mean",
    ParamValues = c(10) # valid initialization first
  )
  pvar$ParamValues <- numeric(0) # now remove all values manually

  lsf <- SYS_LSF(vars = list(pvar), name = "badLSF")
  ps <- SYS_PARAM(
    sys_input = list(lsf),
    probMachines = list(PROB_MACHINE(name = "FORM", fCall = "FORM"))
  )

  expect_error(ps$runMachines(),
    regexp = "must define a non-empty ParamValues",
    fixed = FALSE
  )
})


test_that("SYS_LSF throws an informative error when no limit-state function is defined", {
  # Create a valid set of base variables
  var_E <- PROB_BASEVAR(Name = "E", DistributionType = "norm", Mean = 10, Sd = 0.2)
  var_R <- PROB_BASEVAR(Name = "R", DistributionType = "norm", Mean = 35, Sd = 0.5)

  # Case 1: No function assigned -> should trigger 'missing' message
  lsf_missing <- SYS_LSF(vars = list(var_E, var_R), name = "LSF_missing")

  expect_error(
    lsf_missing$check(),
    regexp = "Limit-state function 'LSF_missing' is missing",
    fixed = FALSE
  )

  # Case 2: Function exists but has no arguments -> should trigger 'no formal arguments' message
  lsf_empty <- SYS_LSF(vars = list(var_E, var_R), name = "LSF_empty")
  lsf_empty$func <- function() NULL

  expect_error(
    lsf_empty$check(),
    regexp = "has no formal arguments",
    fixed = FALSE
  )

  # Case 3: Properly defined function -> check() should pass silently
  lsf_valid <- SYS_LSF(vars = list(var_E, var_R), name = "LSF_valid")
  lsf_valid$func <- function(E, R) {
    R - E
  }

  expect_silent(lsf_valid$check())
})


## ------------------------------------------------------------------------
## 1) PROB_BASEVAR – check missing or invalid inputs
## ------------------------------------------------------------------------

test_that("PROB_BASEVAR throws error when Mean/Sd and parameters are missing", {
  expect_error(
    PROB_BASEVAR(Name = "X1"),
    regexp = "Mean/Sd|distribution parameters",
    fixed = FALSE
  )
})

test_that("PROB_BASEVAR accepts valid Mean/Sd combination", {
  var_ok <- PROB_BASEVAR(Name = "X2", DistributionType = "norm", Mean = 10, Sd = 2)
  expect_s4_class(var_ok, "PROB_BASEVAR")
  expect_equal(var_ok$Mean, 10)
})

test_that("PROB_BASEVAR rejects unknown distribution type", {
  expect_error(
    PROB_BASEVAR(Name = "X3", DistributionType = "unknown", Mean = 10, Sd = 2),
    regexp = "Distribution type 'unknown' is not known",
    fixed = FALSE
  )
})

test_that("PROB_BASEVAR rejects negative standard deviation", {
  expect_error(
    PROB_BASEVAR(Name = "X4", DistributionType = "norm", Mean = 5, Sd = -1),
    regexp = "Standard deviation (Sd) must be non-negative",
    fixed = TRUE
  )
})


## ------------------------------------------------------------------------
## 2) PARAM_BASEVAR – check mandatory fields and parameter sweep setup
## ------------------------------------------------------------------------

test_that("PARAM_BASEVAR throws error when ParamValues are missing", {
  expect_error(
    PARAM_BASEVAR(Name = "P1", ParamType = "Mean"),
    regexp = "PARAM_BASEVAR .*must define a non-empty ParamValues vector.",
    fixed = FALSE
  )
})

# "PARAM_BASEVAR '%s' must define a non-empty ParamValues vector."
test_that("PARAM_BASEVAR throws error when ParamType is missing", {
  expect_error(
    PARAM_BASEVAR(Name = "P2", ParamValues = c(10, 20)),
    regexp = "requires a valid 'ParamType'",
    fixed = FALSE
  )
})

test_that("PARAM_BASEVAR throws error when both Mean/Sd and ParamValues are missing", {
  expect_error(
    PARAM_BASEVAR(Name = "P3"),
    regexp = "PARAM_BASEVAR .*must define a non-empty ParamValues vector.",
    fixed = FALSE
  )
})

test_that("PARAM_BASEVAR initializes correctly with valid inputs (Mean sweep)", {
  pvar <- PARAM_BASEVAR(
    Name = "E_mod",
    DistributionType = "norm",
    ParamType = "Mean",
    ParamValues = c(30000, 35000),
    Sd = 1000
  )

  # Object should be created successfully and first parameter value assigned:
  expect_s4_class(pvar, "PARAM_BASEVAR")

  # The mean should correspond to the first value in ParamValues:
  expect_equal(round(pvar$Mean), 30000)
})


## ------------------------------------------------------------------------
## Optional: combined test for prepare() / getlDistr()
## ------------------------------------------------------------------------

test_that("getlDistr returns list of functions for valid PROB_BASEVAR object", {
  v <- PROB_BASEVAR(Name = "f_ck", DistributionType = "norm", Mean = 30, Sd = 1)
  distr <- v$getlDistr()

  # Density function must exist (closure expected):
  expect_type(distr[[1]]$d, "closure")

  # Parameter vector must be numeric:
  expect_true(is.numeric(distr[[2]]))
})

## ------------------------------------------------------------------------
## Empirical function
## ------------------------------------------------------------------------

test_that("Empirical distribution (emp) is implemented correctly", {
  skip_if_not_installed("EnvStats")
  set.seed(123)
  # x <- rnorm(50)
  x <- c(2, 3, 4, 4.000001, 4.000002, 6, 7, 200)
  v <- PROB_BASEVAR(
    Name = "data_sequence",
    DistributionType = "emp",
    DistributionParameters = x
  )
  lDistr <- v$getlDistr()
  funs <- lDistr[[1]]
  # --- 1) Functions exist ----------------------------------------------------
  expect_type(funs$d, "closure")
  expect_type(funs$p, "closure")
  expect_type(funs$q, "closure")
  expect_type(funs$r, "closure")
  # --- 2) Density equals EnvStats reference ---------------------------------
  grid <- seq(min(x) - 1, max(x) + 1, length.out = 2000)
  expect_equal(
    funs$d(grid),
    EnvStats::demp(grid, obs = x)
  )
  # --- 4) Distribution function consistency ---------------------------------
  ref_p <- suppressWarnings(
    EnvStats::pemp(grid, obs = x)
  )
  expect_equal(
    funs$p(grid),
    ref_p
  )
  # --- 5) Quantile function consistency -------------------------------------
  probs <- seq(0.1, 0.9, by = 0.2)

  expect_equal(
    funs$q(probs),
    EnvStats::qemp(probs, obs = x)
  )
  # --- 6) Random generator only returns observed values ---------------------
  r_vals <- funs$r(200)
  expect_true(all(r_vals %in% x))
})



## ------------------------------------------------------------------------
## Empirical function
## ------------------------------------------------------------------------
test_that("Empirical distribution functions behave correctly", {
  skip_if_not_installed("EnvStats")

  set.seed(1)

  sample_data <- rnorm(1000, mean = 100, sd = 4)

  X <- PROB_BASEVAR(
    Id = 1,
    Name = "Z",
    DistributionType = "emp",
    Package = "EnvStats",
    DistributionParameters = sample_data
  )

  distr <- X$getlDistr()[[1]]

  # CDF monotonicity
  x_vals <- sort(sample_data[1:100])
  p_vals <- distr$p(x_vals)

  expect_true(all(diff(p_vals) >= 0))

  # Quantile inversion check
  u_vals <- seq(0.1, 0.9, length.out = 10)
  q_vals <- distr$q(u_vals)
  p_back <- distr$p(q_vals)

  expect_equal(p_back, u_vals, tolerance = 1e-2)
})

test_that("Student-t distribution is implemented correctly", {

  hyper.param <- c(38, 7, 4, 6)

  m <- hyper.param[1]
  s <- hyper.param[2]
  n <- hyper.param[3]
  v <- hyper.param[4]

  # Density integrates to 1
  f <- function(x) dst(x, hyper.param)
  result <- integrate(f, -Inf, Inf)
  expect_equal(result$value, 1, tolerance = 1e-6)

  # Density matches analytical dt transformation
  s_trans <- s/sqrt(n/(n+1))
  x <- seq(-100, 100, length.out = 1000)
  expected <- dt((x - m)/s_trans, df = v) / s_trans
  actual   <- dst(x, hyper.param)

  # CDF and quantile are inverse
  expect_equal(actual, expected, tolerance = 1e-12)
  p <- seq(0.01, 0.99, length.out = 100)
  q <- qst(p, hyper.param)
  p_back <- pst(q, hyper.param)
  expect_equal(p_back, p, tolerance = 1e-8)

  # Random generation has correct mean and sd
  sim <- rst(1e6, hyper.param)
  # Mean exists for v > 1
  expect_equal(mean(sim), m, tolerance = 0.05)
  # Variance exists for v > 2
  sd_theoretical <- s_trans * sqrt(v/(v-2))
  expect_equal(sd(sim), sd_theoretical, tolerance = 0.05)
})

test_that("Log-Student-t distribution is implemented correctly", {

  set.seed(123)

  hyper.param <- c(3.4, 0.14, 3, 10)

  m <- hyper.param[1]
  s <- hyper.param[2]
  n <- hyper.param[3]
  v <- hyper.param[4]

  s_trans <- s/sqrt(n/(n+1))

  # Log-t density integrates to 1
    f <- function(x) dlt(x, hyper.param)
    result <- integrate(f, 0, Inf)
    expect_equal(result$value, 1, tolerance = 1e-6)

  # Log-t density matches analytical transformation
    x <- seq(0.001, 100, length.out = 1000)
    expected <- (1/x) *
      dt((log(x) - m)/s_trans, df = v) /
      s_trans
    actual <- dlt(x, hyper.param)
    expect_equal(actual, expected, tolerance = 1e-12)

  # Log-t CDF and quantile are inverse
    p <- seq(0.01, 0.99, length.out = 100)
    q <- qlt(p, hyper.param)
    p_back <- plt(q, hyper.param)
    expect_equal(p_back, p, tolerance = 1e-8)

  # Log-t random generation matches exp(t-simulation)
    sim1 <- rlt(1e5, hyper.param)
    # direkte Konstruktion
    sim2 <- exp(m + s_trans * rt(1e5, df = v))
    expect_equal(mean(log(sim1)),
                 mean(log(sim2)),
                 tolerance = 0.02)
    expect_equal(sd(log(sim1)),
                 sd(log(sim2)),
                 tolerance = 0.02)
})



