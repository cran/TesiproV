# =============================================================================
#  test_qualcon.R
#  Tests for PROB_QUALCON (constructor, getOC, isStale, show/summary),
#  .oc_compute_pa() and PROB_BASEVAR$updateBasevar().
#
#  State against which these tests are written:
#   - PROB_QUALCON$initialize(): var and plan (ConformityType/Param) are both
#     optional at construction time - they are only required by getOC().
#     The constructor only warns for an *incomplete* manual OC curve (only
#     one of Pa_grid and theta_grid given, or different lengths). The "ANY"
#     field var: absence is detected via .is_unset() and then normalised to a
#     real NULL.
#   - PROB_QUALCON$getOC(): validates var, DistributionType, ConformityType
#     and ConformityParam BEFORE calling .oc_compute_pa().
#   - .oc_compute_pa(): validates conf.param for "dim_k||n" (c(n, c_accept))
#     and for "fck_EN-206" (c(n[, N]), with f_ck given as lower limit) and
#     explicitly returns P_a in both branches.
#   - updateBasevar() has a parametric fit branch (min_r2) and falls back to
#     the empirical distribution if the fit is poor.
# =============================================================================

make_test_var <- function(DistributionType = "st",
                          DistributionParameters = c(30, 3, 3, 10),
                          Name = "fck_test") {
  PROB_BASEVAR(
    Id = 1, Name = Name,
    DistributionType = DistributionType,
    DistributionParameters = DistributionParameters
  )
}

# -----------------------------------------------------------------------------
#  .oc_compute_pa -- return(P_a) Regressionstests + conf.param-Validierung
# -----------------------------------------------------------------------------

test_that(".oc_compute_pa returns a numeric P_a vector for conf.type = 'dim_k||n'", {
  res <- .oc_compute_pa(
    theta = c(0.01, 0.1, 0.5),
    conf.type = "dim_k||n",
    conf.param = c(20, 1),
    dist.type = "st",
    dist.param = c(30, 3, 3, 10)
  )

  expect_false(is.null(res))
  expect_type(res, "double")
  expect_length(res, 3)
  expect_equal(res, stats::pbinom(1, size = 20, prob = c(0.01, 0.1, 0.5)))
})

test_that(".oc_compute_pa rejects incomplete conf.param for 'dim_k||n'", {
  expect_error(
    .oc_compute_pa(theta = c(0.1), conf.type = "dim_k||n",
                   conf.param = numeric(0), dist.type = "st", dist.param = c(30,3,3,10)),
    "c\\(n, c_accept\\)"
  )
  expect_error(
    .oc_compute_pa(theta = c(0.1), conf.type = "dim_k||n",
                   conf.param = c(20), dist.type = "st", dist.param = c(30,3,3,10)),
    "c\\(n, c_accept\\)"
  )
})

test_that(".oc_compute_pa rejects incomplete conf.param for 'fck_EN-206'", {
  expect_error(
    .oc_compute_pa(theta = c(0.1), conf.type = "fck_EN-206",
                   conf.param = c(15), dist.type = "st", dist.param = c(30,3,3,10)),
    "DefectLimitLower"
  )
  expect_error(
    .oc_compute_pa(theta = c(0.1), conf.type = "fck_EN-206",
                   conf.param = numeric(0), dist.type = "st", dist.param = c(30,3,3,10),
                   lower = 25),
    "must contain n"
  )
})

test_that(".oc_compute_pa errors clearly on unknown conf.type", {
  expect_error(
    .oc_compute_pa(theta = c(0.1), conf.type = "unbekannt",
                   conf.param = c(1), dist.type = "st", dist.param = c(30,3,3,10)),
    "not known"
  )
})

# -----------------------------------------------------------------------------
#  PROB_QUALCON -- Konstruktor
# -----------------------------------------------------------------------------

test_that("PROB_QUALCON constructor does not warn when nothing is set yet", {
  expect_no_warning(PROB_QUALCON(Id = 1))
})

test_that("PROB_QUALCON constructor does not warn when only a plan is given (var set later)", {
  expect_no_warning(
    PROB_QUALCON(Id = 1, ConformityType = "fck_EN-206", ConformityParam = c(15), DefectLimitLower = 3.65)
  )
})

test_that("PROB_QUALCON constructor does not warn when var is given without a plan yet", {
  var <- make_test_var()
  expect_no_warning(PROB_QUALCON(Id = 1, var = var))
})

test_that("PROB_QUALCON constructor warns on a mismatched manual OC curve", {
  expect_warning(
    PROB_QUALCON(Id = 1, Pa_grid = c(0.9, 0.5, 0.1)),
    "supplied together"
  )
  expect_warning(
    PROB_QUALCON(Id = 1, theta_grid = c(0.1, 0.2, 0.3)),
    "supplied together"
  )
  expect_warning(
    PROB_QUALCON(Id = 1, Pa_grid = c(0.9, 0.5), theta_grid = c(0.1, 0.2, 0.3)),
    "different lengths"
  )
})

test_that("PROB_QUALCON constructor accepts a matching manual OC curve without warning", {
  expect_no_warning(
    PROB_QUALCON(Id = 1, theta_grid = c(0.01, 0.1, 0.5), Pa_grid = c(0.95, 0.6, 0.1))
  )
})

test_that("PROB_QUALCON constructor normalizes an unset $var to real NULL", {
  qc <- PROB_QUALCON(Id = 1, ConformityType = "fck_EN-206", ConformityParam = c(15), DefectLimitLower = 3.65)
  expect_true(is.null(qc$var))
  expect_identical(class(qc$var), "NULL")
})

test_that("PROB_QUALCON constructor rejects a var that is not a PROB_BASEVAR", {
  expect_error(
    PROB_QUALCON(
      Id = 1, ConformityType = "dim_k||n", ConformityParam = c(20, 1),
      var = list(not = "a basevar")
    ),
    "PROB_BASEVAR"
  )
})

# -----------------------------------------------------------------------------
#  PROB_QUALCON$getOC() -- Validierung
# -----------------------------------------------------------------------------

test_that("PROB_QUALCON$getOC requires $var to be attached", {
  qc <- PROB_QUALCON(Id = 1, ConformityType = "dim_k||n", ConformityParam = c(20, 1))
  expect_error(qc$getOC(), "No PROB_BASEVAR attached")
})

test_that("PROB_QUALCON$getOC requires var$DistributionType 'st' or 'lt'", {
  var_norm <- PROB_BASEVAR(Id = 1, Name = "x", DistributionType = "norm", Mean = 10, Sd = 1)
  qc <- PROB_QUALCON(
    Id = 1, ConformityType = "dim_k||n", ConformityParam = c(20, 1),
    var = var_norm
  )
  expect_error(qc$getOC(), "'st' or 'lt'")
})

test_that("PROB_QUALCON$getOC requires ConformityType to be set", {
  var <- make_test_var()
  qc <- PROB_QUALCON(Id = 1, var = var)
  expect_error(qc$getOC(), "ConformityType not set")
})

test_that("PROB_QUALCON$getOC requires ConformityParam to be set", {
  var <- make_test_var()
  qc <- PROB_QUALCON(Id = 1, var = var, ConformityType = "dim_k||n")
  expect_error(qc$getOC(), "ConformityParam not set")
})

test_that("PROB_QUALCON$getOC computes Pa_grid for a dim_k||n plan", {
  var <- make_test_var()
  qc <- PROB_QUALCON(
    Id = 1, Name = "attribute_plan",
    ConformityType = "dim_k||n", ConformityParam = c(20, 1),
    var = var
  )

  expect_true(all(is.na(qc$Pa_grid)))

  qc$getOC()

  expected <- stats::pbinom(1, size = 20, prob = qc$theta_grid)
  expect_equal(qc$Pa_grid, expected)
})

test_that("PROB_QUALCON$getOC does not recompute without force", {
  var <- make_test_var()
  qc <- PROB_QUALCON(
    Id = 1, ConformityType = "dim_k||n", ConformityParam = c(20, 1), var = var
  )
  qc$getOC()
  first <- qc$Pa_grid

  expect_message(qc$getOC(), "already available")
  expect_equal(qc$Pa_grid, first)
})

test_that("PROB_QUALCON$getOC recomputes with force = TRUE", {
  var <- make_test_var()
  qc <- PROB_QUALCON(
    Id = 1, ConformityType = "dim_k||n", ConformityParam = c(20, 1), var = var
  )
  qc$getOC()

  qc$ConformityParam <- c(20, 2)
  expect_message(qc$getOC(force = TRUE), "OC curve computed")

  expected <- stats::pbinom(2, size = 20, prob = qc$theta_grid)
  expect_equal(qc$Pa_grid, expected)
})

test_that("PROB_QUALCON$getOC with force = TRUE uses the default theta_grid again", {
  var <- make_test_var()
  qc <- PROB_QUALCON(
    Id = 1, ConformityType = "dim_k||n", ConformityParam = c(25, 1),
    DefectLimitLower = 27, DefectLimitUpper = 33,
    theta_grid = c(0, 1), Pa_grid = c(1, 1), var = var
  )
  suppressMessages(qc$getOC(force = TRUE))

  default_grid <- c(0.001, 0.005, 0.007, 0.01, 0.012, 0.015, 0.02, 0.03, 0.05,
                    0.07, 0.09, 0.12, 0.15, 0.2, 0.3, 0.4, 0.5, 0.9)
  expect_equal(qc$theta_grid, default_grid)
  expect_equal(qc$Pa_grid, stats::pbinom(1, size = 25, prob = default_grid))
})

test_that("PROB_QUALCON has NSim with default 100000 and validates it", {
  var <- make_test_var()
  qc <- PROB_QUALCON(
    Id = 1, ConformityType = "fck_EN-206", ConformityParam = 15,
    DefectLimitLower = 30, var = var
  )
  expect_equal(qc$NSim, 100000)

  qc2 <- PROB_QUALCON(
    Id = 2, ConformityType = "fck_EN-206", ConformityParam = 15, NSim = 2000,
    DefectLimitLower = 30, var = var
  )
  expect_equal(qc2$NSim, 2000)

  expect_error(
    PROB_QUALCON(Id = 3, ConformityType = "fck_EN-206", ConformityParam = 15,
                 NSim = 0, DefectLimitLower = 30, var = var),
    "NSim"
  )
})

test_that("PROB_QUALCON rejects the sample size in ConformityParam for 'fck_EN-206'", {
  var <- make_test_var()
  expect_error(
    PROB_QUALCON(Id = 1, ConformityType = "fck_EN-206", ConformityParam = c(15, 2000),
                 DefectLimitLower = 30, var = var),
    "NSim"
  )
})

test_that(".oc_compute_pa rejects c(n, N) and an invalid N for 'fck_EN-206'", {
  expect_error(
    .oc_compute_pa(theta = c(0.1), conf.type = "fck_EN-206",
                   conf.param = c(15, 2000), dist.type = "st",
                   dist.param = c(30, 3, 3, 10), lower = 25),
    "argument N"
  )
  expect_error(
    .oc_compute_pa(theta = c(0.1), conf.type = "fck_EN-206",
                   conf.param = 15, dist.type = "st",
                   dist.param = c(30, 3, 3, 10), lower = 25, N = -1),
    "positive"
  )
})

# -----------------------------------------------------------------------------
#  PROB_QUALCON$isStale()
# -----------------------------------------------------------------------------

test_that("PROB_QUALCON$isStale is FALSE right after getOC()", {
  var <- make_test_var()
  qc <- PROB_QUALCON(
    Id = 1, ConformityType = "dim_k||n", ConformityParam = c(20, 1), var = var
  )
  qc$getOC()
  expect_false(qc$isStale())
})

test_that("PROB_QUALCON$isStale detects a changed $var after getOC()", {
  var <- make_test_var()
  qc <- PROB_QUALCON(
    Id = 1, ConformityType = "dim_k||n", ConformityParam = c(20, 1), var = var
  )
  qc$getOC()

  var$DistributionParameters <- c(35, 2, 3, 10)
  var$prepare()

  expect_true(qc$isStale())
})

test_that("PROB_QUALCON$isStale returns NA without a prior snapshot", {
  var <- make_test_var()
  qc_manual <- PROB_QUALCON(
    Id = 1, var = var,
    theta_grid = c(0.01, 0.1, 0.5),
    Pa_grid    = c(0.95, 0.6, 0.1)
  )
  expect_true(is.na(qc_manual$isStale()))
})

# -----------------------------------------------------------------------------
#  PROB_QUALCON$show() / $summary()
# -----------------------------------------------------------------------------

test_that("PROB_QUALCON print output is compact and shows only var's name", {
  var <- make_test_var(Name = "fck")
  qc <- PROB_QUALCON(
    Id = 1, Name = "attribute_plan",
    ConformityType = "dim_k||n", ConformityParam = c(20, 1), var = var
  )
  qc$getOC()

  out <- capture.output(print(qc))
  expect_true(any(grepl("var\\s*:\\s*fck", out)))
  expect_false(any(grepl("Field \"DistributionType\"", out)))
  expect_false(any(grepl("Field \"\\.cache\"", out)))
})

test_that("PROB_QUALCON summary reports point count once Pa_grid is available", {
  var <- make_test_var()
  qc <- PROB_QUALCON(
    Id = 1, ConformityType = "dim_k||n", ConformityParam = c(20, 1), var = var
  )

  out_before <- capture.output(qc$summary())
  expect_true(any(grepl("not yet computed", out_before)))

  qc$getOC()
  out_after <- capture.output(qc$summary())
  expect_true(any(grepl(sprintf("available \\(\\s*%d\\s*points\\)", length(qc$theta_grid)), out_after)))
})

# -----------------------------------------------------------------------------
#  PROB_BASEVAR$updateBasevar()
# -----------------------------------------------------------------------------
#  Both plan types are used here.

test_that("updateBasevar rejects a var with unsupported DistributionType", {
  var_norm <- PROB_BASEVAR(Id = 1, Name = "x", DistributionType = "norm", Mean = 10, Sd = 1)
  qc <- PROB_QUALCON(
    Id = 1,
    ConformityType  = "fck_EN-206",
    ConformityParam = c(15),
    DefectLimitLower = 10,
    var = var_norm,
    theta_grid = c(0.01, 0.5), Pa_grid = c(0.9, 0.1)
  )
  expect_error(var_norm$updateBasevar(qc), "'st' or 'lt'")
})

test_that("updateBasevar rejects a PROB_QUALCON of the wrong class", {
  var_prior <- PROB_BASEVAR(
    Id = 1, Name = "fck", DistributionType = "lt",
    DistributionParameters = c(3.65, 0.12, 3, 10)
  )
  expect_error(var_prior$updateBasevar(list()), "PROB_QUALCON")
})

test_that("updateBasevar returns a PROB_BASEVAR copy and leaves the prior unmodified", {
  var_prior <- PROB_BASEVAR(
    Id = 1, Name = "fck",
    DistributionType = "lt",
    DistributionParameters = c(3.65, 0.12, 3, 10)
  )

  qc <- PROB_QUALCON(
    Id = 1, Name = "manual_oc",
    ConformityType  = "fck_EN-206",
    ConformityParam = c(15),
    DefectLimitLower = 3.65,
    var             = var_prior,
    theta_grid = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 0.9),
    Pa_grid    = c(0.99,  0.95, 0.8,  0.6, 0.3, 0.05, 0.001)
  )

  prior_type   <- var_prior$DistributionType
  prior_params <- var_prior$DistributionParameters
  prior_mean   <- var_prior$Mean

  set.seed(1)
  var_post <- var_prior$updateBasevar(qc, min_r2 = 0)  # min_r2 = 0 -> the fit is almost always accepted

  expect_s4_class(var_post, "PROB_BASEVAR")
  expect_true(is.finite(var_post$Mean))
  expect_true(is.finite(var_post$Sd))

  expect_identical(var_prior$DistributionType, prior_type)
  expect_equal(var_prior$DistributionParameters, prior_params)
  expect_equal(var_prior$Mean, prior_mean)

  expect_false(identical(var_post$.cache, var_prior$.cache))
})

test_that("updateBasevar accepts the parametric fit when min_r2 is permissive", {
  var_prior <- PROB_BASEVAR(
    Id = 1, Name = "fck", DistributionType = "lt",
    DistributionParameters = c(3.65, 0.12, 3, 10)
  )
  qc <- PROB_QUALCON(
    Id = 1,
    ConformityType  = "fck_EN-206",
    ConformityParam = c(15),
    DefectLimitLower = 3.65,
    var = var_prior,
    theta_grid = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 0.9),
    Pa_grid    = c(0.99,  0.95, 0.8,  0.6, 0.3, 0.05, 0.001)
  )

  set.seed(1)
  var_post <- var_prior$updateBasevar(qc, min_r2 = 0)

  expect_identical(var_post$DistributionType, "lt")
  expect_length(var_post$DistributionParameters, 4)
})

test_that("updateBasevar falls back to the empirical distribution when min_r2 is unreachable", {
  var_prior <- PROB_BASEVAR(
    Id = 1, Name = "fck", DistributionType = "lt",
    DistributionParameters = c(3.65, 0.12, 3, 10)
  )
  qc <- PROB_QUALCON(
    Id = 1,
    ConformityType  = "fck_EN-206",
    ConformityParam = c(15),
    DefectLimitLower = 3.65,
    var = var_prior,
    theta_grid = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 0.9),
    Pa_grid    = c(0.99,  0.95, 0.8,  0.6, 0.3, 0.05, 0.001)
  )

  set.seed(1)
  var_post <- var_prior$updateBasevar(qc, min_r2 = 0.999999999)

  expect_identical(var_post$DistributionType, "emp")
  expect_length(var_post$DistributionParameters, 100000)
  expect_true(is.numeric(var_post$DistributionParameters))
})

test_that("updateBasevar warns when the PROB_QUALCON is stale relative to var", {
  var_prior <- PROB_BASEVAR(
    Id = 1, Name = "fck", DistributionType = "lt",
    DistributionParameters = c(3.65, 0.12, 3, 10)
  )

  qc <- PROB_QUALCON(
    Id = 1, ConformityType = "dim_k||n", ConformityParam = c(20, 1),
    var = var_prior
  )
  qc$getOC()
  expect_false(qc$isStale())

  var_prior$DistributionParameters <- c(3.98, 0.07, 3, 10)
  var_prior$prepare()
  expect_true(qc$isStale())

  qc$ConformityType  <- "fck_EN-206"
  qc$ConformityParam <- c(15)
  qc$DefectLimitLower <- 3.98
  qc$theta_grid <- c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 0.9)
  qc$Pa_grid    <- c(0.99,  0.95, 0.8,  0.6, 0.3, 0.05, 0.001)

  set.seed(1)
  expect_warning(var_prior$updateBasevar(qc, min_r2 = 0), "changed since getOC")
})
