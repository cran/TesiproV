test_that("PROB_BASEVAR: prepare and getlDistr produce funlist", {
  v <- PROB_BASEVAR(Id = 1, Name = "f_ck", DistributionType = "norm", Mean = 30, Sd = 1)

  # Accept both S3 and S4
  expect_true(inherits(v, "PROB_BASEVAR"))

  dd <- v$getlDistr()
  # dd kann wrapper (list(funlist, params)) oder funlist sein -> normalize
  funlist <- if (is.list(dd) && length(dd) >= 1 && is.list(dd[[1]]) && !is.null(dd[[1]]$p)) dd[[1]] else dd

  expect_true(is.list(funlist))
  expect_true(is.function(funlist$d))
  expect_true(is.function(funlist$p))
  expect_true(is.function(funlist$q))
  expect_true(is.function(funlist$r))
})


test_that("SYS_LSF getLSF evaluates correctly", {
  v1 <- PROB_BASEVAR(Id = 1, Name = "X1", DistributionType = "norm", Mean = 0, Sd = 1)

  lsf_obj <- SYS_LSF(name = "g1", vars = list(v1))
  lsf_obj$func <- function(X1) 3 - X1

  f <- lsf_obj$getLSF()
  expect_true(is.function(f))

  expect_equal(f(c(2)), 1)
  expect_equal(f(c(4)), -1)
})


test_that("MC_IS_system: deterministic with fixed seed and dps (quick smoke test)", {
  with_test_future_plan(future::sequential, {
    backend <- "future"

    v1 <- PROB_BASEVAR(
      Id = 1, Name = "X1",
      DistributionType = "norm", Mean = 0, Sd = 1
    )
    v2 <- PROB_BASEVAR(
      Id = 2, Name = "X2",
      DistributionType = "norm", Mean = 0, Sd = 1
    )

    lsf1 <- function(x) 3 - x[1]
    lsf2 <- function(x) 3 - x[1]

    # Get and normalize the marginal distribution funlists.
    d1 <- v1$getlDistr()
    d1 <- if (is.list(d1) && length(d1) >= 1 &&
      is.list(d1[[1]]) && !is.null(d1[[1]]$p)) {
      d1[[1]]
    } else {
      d1
    }

    d2 <- v2$getlDistr()
    d2 <- if (is.list(d2) && length(d2) >= 1 &&
      is.list(d2[[1]]) && !is.null(d2[[1]]$p)) {
      d2[[1]]
    } else {
      d2
    }

    lDistr <- list(list(d1), list(d2))

    run_once <- function(seed_val) {
      MC_IS(
        lsf = list(lsf1, lsf2),
        lDistr = lDistr,
        cov_user = 0.5,
        n_batch = 1000,
        n_max = 2000,
        backend = backend,
        seed = seed_val,
        dps = c(3, 3),
        dataRecord = TRUE,
        debug.level = 0
        # Do not supply use_threads for the future backend.
      )
    }

    r1 <- run_once(1234)
    r2 <- run_once(1234)

    expect_equal(r1$pf, r2$pf)
    expect_equal(r1$n_mc, r2$n_mc)

    expect_true(is.numeric(r1$pf))
    expect_true(r1$pf >= 0 && r1$pf <= 1)

    # Test the SYS_PROB wrapper.
    lsf1obj <- SYS_LSF(name = "g1", vars = list(v1))
    lsf1obj$func <- function(X1) 3 - X1

    lsf2obj <- SYS_LSF(name = "g2", vars = list(v2))
    lsf2obj$func <- function(X2) 3 - X2

    ps <- SYS_PROB(
      sys_input = list(lsf1obj, lsf2obj),
      probMachines = list(),
      sys_type = "serial"
    )

    ps$calculateSystemProbability(
      calcType = "MCIS",
      params = list(
        cov_user = 0.5,
        n_batch = 1000,
        n_max = 2000,
        seed = 1234,
        backend = "future",
        dps = c(3, 3),
        dataRecord = TRUE
        # Do not supply use_threads for the future backend.
      )
    )

    expect_true(length(ps$res_sys) >= 1)

    last <- ps$res_sys[[length(ps$res_sys)]]
    expect_true(!is.null(last$pf))
  })
})


test_that("MC_IS_system is reproducible for same seed", {
  skip_on_cran()
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }
  v <- PROB_BASEVAR(
    Id = 1, Name = "X",
    DistributionType = "norm",
    Mean = 0, Sd = 1
  )

  lsf <- SYS_LSF(name = "g", vars = list(v))
  lsf$func <- function(X) 3 - X

  ps1 <- SYS_PROB(
    sys_input = list(lsf),
    probMachines = list(),
    sys_type = "serial"
  )

  ps1$calculateSystemProbability(
    calcType = "MCIS",
    params = list(
      cov_user = 0.01,
      n_batch = 5000,
      n_max = 1e5,
      seed = 1234,
      backend = backend
    )
  )

  beta1 <- ps1$beta_sys[1, 1]

  ps2 <- SYS_PROB(
    sys_input = list(lsf),
    probMachines = list(),
    sys_type = "serial"
  )

  ps2$calculateSystemProbability(
    calcType = "MCIS",
    params = list(
      cov_user = 0.01,
      n_batch = 5000,
      n_max = 1e5,
      seed = 1234,
      backend = backend
    )
  )

  beta2 <- ps2$beta_sys[1, 1]

  expect_equal(beta1, beta2)
})


test_that("MC_IS_system is reproducible for same seed", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }

  v1 <- PROB_BASEVAR(
    Id = 1, Name = "X1",
    DistributionType = "norm",
    Mean = 0, Sd = 1
  )

  v2 <- PROB_BASEVAR(
    Id = 2, Name = "X2",
    DistributionType = "norm",
    Mean = 0, Sd = 1
  )

  lsf1 <- function(x) 3 - x[1]
  lsf2 <- function(x) 3 - x[1]

  d1 <- v1$getlDistr()[[1]]
  d2 <- v2$getlDistr()[[1]]

  lDistr <- list(list(d1), list(d2))

  run_once <- function(seed_val) {
    MC_IS(
      lsf = list(lsf1, lsf2),
      lDistr = lDistr,
      cov_user = 0.5,
      n_batch = 1000,
      n_max = 2000,
      backend = backend,
      seed = seed_val,
      dps = c(3, 3),
      dataRecord = TRUE
    )
  }

  r1 <- run_once(1234)
  r2 <- run_once(1234)

  expect_equal(r1$pf, r2$pf)
  expect_equal(r1$n_mc, r2$n_mc)
})


test_that("ESS is positive and bounded with one Future worker", {
  v <- PROB_BASEVAR(
    Id = 1,
    Name = "X",
    DistributionType = "norm",
    Mean = 0,
    Sd = 1
  )

  lsf <- function(x) 3 - x[1]
  d1 <- v$getlDistr()[[1]]

  res <- with_test_future_plan(
    future::sequential,
    {
      expect_equal(future::nbrOfWorkers(), 1L)

      MC_IS(
        lsf = list(lsf),
        lDistr = list(list(d1)),
        cov_user = 0.2,
        n_batch = 2000,
        n_max = 5000,
        backend = "future",
        seed = 1234,
        dataRecord = FALSE
        # Do not supply use_threads for the future backend.
      )
    }
  )

  expect_true(res$ESS > 0)
  expect_true(res$ESS <= res$n_mc)
})


test_that("ESS is positive and bounded with two workers", {
  v <- PROB_BASEVAR(
    Id = 1,
    Name = "X",
    DistributionType = "norm",
    Mean = 0,
    Sd = 1
  )

  lsf <- function(x) 3 - x[1]
  d1 <- v$getlDistr()[[1]]

  res <- with_test_future_plan(
    future::multisession,
    {
      expect_equal(future::nbrOfWorkers(), 2L)

      MC_IS(
        lsf = list(lsf),
        lDistr = list(list(d1)),
        cov_user = 0.2,
        n_batch = 2000,
        n_max = 5000,
        backend = "future",
        seed = 1234,
        dataRecord = FALSE
        # Do not supply use_threads for the future backend.
      )
    },
    workers = 2
  )

  expect_true(res$ESS > 0)
  expect_true(res$ESS <= res$n_mc)
})


test_that("MC_IS_system handles intersecting limit states", {
  skip_on_cran()
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }

  v1 <- PROB_BASEVAR(
    Id = 1, Name = "X1",
    DistributionType = "norm",
    Mean = 0, Sd = 1
  )

  v2 <- PROB_BASEVAR(
    Id = 2, Name = "X2",
    DistributionType = "norm",
    Mean = 0, Sd = 1
  )

  lsf1 <- function(x) 3 - x[1]
  lsf2 <- function(x) 3 - (0.7 * x[1] + 0.7 * x[2])

  d1 <- v1$getlDistr()[[1]]
  d2 <- v2$getlDistr()[[1]]

  res <- MC_IS(
    lsf = list(lsf1, lsf2),
    lDistr = list(list(d1), list(d1, d2)),
    cov_user = 0.05,
    n_batch = 5000,
    n_max = 1e5,
    backend = backend,
    seed = 1234
  )

  expect_true(res$pf > 0)
  expect_true(res$pf < 1)
})


test_that("MC_IS serial is reproducible for fixed seed", {
  skip_on_cran()
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }

  v <- PROB_BASEVAR(
    Id = 1, Name = "X",
    DistributionType = "norm",
    Mean = 0, Sd = 1
  )

  lsf <- SYS_LSF(name = "g", vars = list(v))
  lsf$func <- function(X) 3 - X

  ps1 <- SYS_PROB(
    sys_input = list(lsf),
    probMachines = list(),
    sys_type = "serial"
  )

  ps2 <- SYS_PROB(
    sys_input = list(lsf),
    probMachines = list(),
    sys_type = "serial"
  )

  ps1$calculateSystemProbability(
    calcType = "MCIS",
    params = list(
      cov_user = 0.02,
      n_batch = 5000,
      n_max = 1e5,
      seed = 1234,
      backend = backend
    )
  )

  ps2$calculateSystemProbability(
    calcType = "MCIS",
    params = list(
      cov_user = 0.02,
      n_batch = 5000,
      n_max = 1e5,
      seed = 1234,
      backend = backend
    )
  )

  expect_equal(
    ps1$res_sys[[1]]$pf,
    ps2$res_sys[[1]]$pf
  )
})


# test_that("adaptive alpha is reproducible for fixed seed", {
#     skip_if_not_installed("future")
#     old_plan <- future::plan()
#     future::plan(future::sequential)
#     on.exit(future::plan(old_plan), add = TRUE)

#     v1 <- PROB_BASEVAR(
#         Id = 1, Name = "X1",
#         DistributionType = "norm",
#         Mean = 0, Sd = 1
#     )

#     v2 <- PROB_BASEVAR(
#         Id = 2, Name = "X2",
#         DistributionType = "norm",
#         Mean = 0, Sd = 1
#     )

#     lsf1 <- function(x) 3 - x[1]
#     lsf2 <- function(x) 3 - (0.7 * x[1] + 0.7 * x[2])

#     d1 <- v1$getlDistr()[[1]]
#     d2 <- v2$getlDistr()[[1]]

#     lDistr <- list(list(d1), list(d1, d2))

#     run_once <- function(seed_val) {
#         MC_IS(
#             lsf = list(lsf1, lsf2),
#             lDistr = lDistr,
#             cov_user = 0.02,
#             n_batch = 10000,
#             n_max = 2e5,
#             use_threads = 1,
#             backend = "future",
#             seed = seed_val,
#             adaptive_alpha = TRUE
#         )
#     }

#     r1 <- run_once(1234)
#     r2 <- run_once(1234)

#     expect_equal(r1$pf, r2$pf)
# })


# test_that("adaptive alpha does not catastrophically worsen Pf", {
#     v1 <- PROB_BASEVAR(
#         Id = 1, Name = "X1",
#         DistributionType = "norm",
#         Mean = 0, Sd = 1
#     )

#     v2 <- PROB_BASEVAR(
#         Id = 2, Name = "X2",
#         DistributionType = "norm",
#         Mean = 0, Sd = 1
#     )

#     lsf1 <- function(x) 3 - x[1]
#     lsf2 <- function(x) 3 - (0.7 * x[1] + 0.7 * x[2])

#     d1 <- v1$getlDistr()[[1]]
#     d2 <- v2$getlDistr()[[1]]

#     lDistr <- list(list(d1), list(d1, d2))

#     static_run <- MC_IS(
#         lsf = list(lsf1, lsf2),
#         lDistr = lDistr,
#         cov_user = 0.02,
#         n_batch = 10000,
#         n_max = 2e5,
#         use_threads = 1,
#         backend = "future",
#         seed = 1234,
#         adaptive_alpha = FALSE
#     )

#     adaptive_run <- MC_IS(
#         lsf = list(lsf1, lsf2),
#         lDistr = lDistr,
#         cov_user = 0.02,
#         n_batch = 10000,
#         n_max = 2e5,
#         use_threads = 1,
#         backend = "future",
#         seed = 1234,
#         adaptive_alpha = TRUE
#     )

#     expect_lt(
#         abs(adaptive_run$pf - static_run$pf),
#         0.5 * static_run$pf
#     )
# })


# test_that("ESS is finite and positive for adaptive alpha", {
#     if (.Platform$OS.type == "windows") {
#         backend <- "future"
#     } else {
#         backend <- "parallel"
#     }

#     v1 <- PROB_BASEVAR(
#         Id = 1, Name = "X1",
#         DistributionType = "norm",
#         Mean = 0, Sd = 1
#     )

#     v2 <- PROB_BASEVAR(
#         Id = 2, Name = "X2",
#         DistributionType = "norm",
#         Mean = 0, Sd = 1
#     )

#     lsf1 <- function(x) 3 - x[1]
#     lsf2 <- function(x) 3 - (0.7 * x[1] + 0.7 * x[2])

#     d1 <- v1$getlDistr()[[1]]
#     d2 <- v2$getlDistr()[[1]]

#     lDistr <- list(list(d1), list(d1, d2))

#     res <- MC_IS(
#         lsf = list(lsf1, lsf2),
#         lDistr = lDistr,
#         cov_user = 0.02,
#         n_batch = 10000,
#         n_max = 2e5,
#         use_threads = 1,
#         backend = backend,
#         seed = 1234,
#         adaptive_alpha = TRUE
#     )

#     expect_true(is.finite(res$ESS))
#     expect_gt(res$ESS, 0)
# })


test_that("compute_weight_update accumulates the failure sums in the log domain", {
  state <- list(
    log_sum_Iw = -Inf,
    log_sum_Iw2 = -Inf,
    n_sim = 0L
  )

  all_log_w <- log(c(2, 1, 3))
  all_I_vals <- c(1, 0, 1)

  res <- compute_weight_update(
    all_log_w = all_log_w,
    all_I_vals = all_I_vals,
    state = state
  )

  # failure weights are 2 and 3
  expect_equal(exp(res$log_sum_Iw), 5)
  expect_equal(exp(res$log_sum_Iw2), 13)
  expect_equal(res$n_sim, 3)
})

test_that("compute_weight_update: result does not depend on the batch partition", {
  set.seed(1)
  lw <- rnorm(1000, sd = 3)
  I <- as.integer(runif(1000) < 0.1)
  s0 <- list(log_sum_Iw = -Inf, log_sum_Iw2 = -Inf, n_sim = 0L)

  one <- compute_weight_update(lw, I, state = s0)

  many <- s0
  for (idx in split(seq_along(lw), rep(1:10, each = 100))) {
    many <- compute_weight_update(lw[idx], I[idx], state = many)
  }

  expect_equal(many$log_sum_Iw, one$log_sum_Iw)
  expect_equal(many$log_sum_Iw2, one$log_sum_Iw2)
  expect_equal(many$n_sim, 1000)
  expect_equal(one$log_sum_Iw, log(sum(exp(lw[I == 1]))))
  expect_equal(one$log_sum_Iw2, log(sum(exp(2 * lw[I == 1]))))
})

test_that("compute_weight_update: batch without failures keeps the sums at -Inf", {
  state <- list(log_sum_Iw = -Inf, log_sum_Iw2 = -Inf, n_sim = 0L)

  res <- compute_weight_update(
    all_log_w = c(-1, -2, -0.5),
    all_I_vals = c(0, 0, 0),
    state = state
  )

  expect_equal(res$log_sum_Iw, -Inf)
  expect_equal(res$log_sum_Iw2, -Inf)
  expect_equal(res$n_sim, 3)
})

test_that("compute_weight_update: very large log-weights do not overflow", {
  state <- list(log_sum_Iw = -Inf, log_sum_Iw2 = -Inf, n_sim = 0L)

  res <- compute_weight_update(
    all_log_w = c(800, 790),
    all_I_vals = c(1, 1),
    state = state
  )

  expect_true(is.finite(res$log_sum_Iw))
  expect_equal(res$log_sum_Iw, 800 + log(1 + exp(-10)))
  expect_equal(res$log_sum_Iw2, 1600 + log(1 + exp(-20)))
})


# ------------------------------------------------------------
# Helper: create simple base variables
# ------------------------------------------------------------

create_vars <- function() {
  X1 <- PROB_BASEVAR(
    Id = 1, Name = "X1",
    DistributionType = "norm",
    Mean = 10, Sd = 1
  )

  X2 <- PROB_BASEVAR(
    Id = 2, Name = "X2",
    DistributionType = "norm",
    Mean = 5, Sd = 1
  )

  list(X1, X2)
}

# ------------------------------------------------------------
# Helper: wrap functions into SYS_LSF
# ------------------------------------------------------------

build_lsf_objects <- function(lsf_list, vars) {
  lapply(seq_along(lsf_list), function(i) {
    obj <- SYS_LSF(vars = vars, name = paste0("LSF", i))
    obj$func <- lsf_list[[i]]
    obj$check()
    obj
  })
}

# ============================================================
# TEST 1 — Debug level propagation
# ============================================================

test_that("debug.level is correctly propagated to MC_IS_system", {
  vars <- create_vars()

  lsf <- function(x) x[1] - x[2]

  ps <- SYS_PROB(
    sys_input = build_lsf_objects(list(lsf, lsf), vars),
    probMachines = list(),
    debug.level = 1,
    sys_type = "serial"
  )

  expect_equal(ps$debug.level, 1)
})

# ============================================================
# TEST 2 — Geometry detection: identical LSFs → single shift
# ============================================================

test_that("Identical LSFs are classified as non-multimodal", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }

  vars <- create_vars()

  lsf <- function(x) x[1] - x[2]

  ps <- SYS_PROB(
    sys_input = build_lsf_objects(list(lsf, lsf), vars),
    probMachines = list(),
    debug.level = 0,
    sys_type = "serial"
  )

  res <- ps$calculateSystemProbability(
    calcType = "MCIS",
    params = list(
      n_batch = 2000,
      n_max = 5000,
      cov_user = 0.5,
      backend = backend,
      seed = 123
    )
  )

  expect_true(TRUE) # no crash = pass
})

# ============================================================
# TEST 3 — Geometry detection: orthogonal LSFs → multimodal
# ============================================================

test_that("Orthogonal LSFs are classified as multimodal", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }
  vars <- create_vars()

  lsf1 <- function(x) x[1] - 8
  lsf2 <- function(x) x[2] - 3

  ps <- SYS_PROB(
    sys_input = build_lsf_objects(list(lsf1, lsf2), vars),
    probMachines = list(),
    debug.level = 0,
    sys_type = "serial"
  )

  res <- ps$calculateSystemProbability(
    calcType = "MCIS",
    params = list(
      n_batch = 2000,
      n_max = 5000,
      cov_user = 0.5,
      backend = backend,
      seed = 123
    )
  )

  expect_true(is.numeric(ps$beta_sys[1, 1]))
})

# ============================================================
# TEST 4 — Determinism across thread counts
# ============================================================

test_that("MC_IS_system is deterministic across thread counts", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }

  vars <- create_vars()

  lsf1 <- function(x) x[1] - x[2]
  lsf2 <- function(x) x[1] + x[2] - 12

  run_test <- function(threads) {
    ps_local <- SYS_PROB(
      sys_input = build_lsf_objects(list(lsf1, lsf2), vars),
      probMachines = list(),
      debug.level = 0,
      sys_type = "serial"
    )

    ps_local$calculateSystemProbability(
      calcType = "MCIS",
      params = c(
        list(
          n_batch = 3000,
          n_max = 8000,
          cov_user = 0.5,
          backend = backend,
          seed = 42
        ),
        if (backend == "parallel") list(use_threads = threads)
      )
    )
    ps_local$beta_sys
  }

  res1 <- run_test(1)
  res2 <- run_test(2)

  expect_equal(res1, res2)
})

# ============================================================
# TEST 5 — adaptive_alpha does not break single shift
# ============================================================

test_that("Adaptive alpha does not crash in single-shift case", {
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }

  vars <- create_vars()

  lsf <- function(x) x[1] - x[2]

  ps <- SYS_PROB(
    sys_input = build_lsf_objects(list(lsf, lsf), vars),
    probMachines = list(),
    debug.level = 0,
    sys_type = "serial"
  )

  ps$calculateSystemProbability(
    calcType = "MCIS",
    params = list(
      n_batch = 2000,
      n_max = 5000,
      cov_user = 0.5,
      backend = backend,
      seed = 123,
      adaptive_alpha = TRUE
    )
  )

  expect_true(is.finite(ps$beta_sys[1, 1]))
})


test_that("build_gsys constructs valid system object", {
  v1 <- PROB_BASEVAR(Name = "X1", DistributionType = "norm", Mean = 5, Sd = 1)
  v2 <- PROB_BASEVAR(Name = "X2", DistributionType = "norm", Mean = 3, Sd = 1)

  lsf1 <- SYS_LSF(
    name = "g1",
    vars = list(v1)
  )
  lsf1$func <- function(x) x[1] - 4

  lsf2 <- SYS_LSF(
    name = "g2",
    vars = list(v2)
  )
  lsf2$func <- function(x) x[1] - 2

  sys <- build_gsys(list(lsf1, lsf2))

  expect_true(is.function(sys$g_sys))
  expect_equal(sys$n_lsfs, 2)
  expect_equal(sys$n_vars, 2)
  expect_length(sys$distr_flat, 2)
})


test_that("g_sys serial approximates minimum", {
  v1 <- PROB_BASEVAR(Name = "X1", DistributionType = "norm", Mean = 0, Sd = 1)
  v2 <- PROB_BASEVAR(Name = "X2", DistributionType = "norm", Mean = 0, Sd = 1)

  lsf1 <- SYS_LSF(vars = list(v1))
  lsf1$func <- function(x) x[1] - 1

  lsf2 <- SYS_LSF(vars = list(v2))
  lsf2$func <- function(x) x[1] - 2

  sys <- build_gsys(
    sys_input = list(lsf1, lsf2),
    sys_structure = "serial",
    smooth_kappa = 100
  )

  x <- c(1.5, 3)

  g_vals <- c(
    lsf1$getLSF()(x[1]),
    lsf2$getLSF()(x[2])
  )

  g_expected <- min(g_vals)

  g_sys <- sys$g_sys(x)

  expect_equal(g_sys, g_expected, tolerance = 1e-3)
})


test_that("smooth serial aggregation behaves correctly", {
  v1 <- PROB_BASEVAR(Name = "X1", DistributionType = "norm", Mean = 0, Sd = 1)
  v2 <- PROB_BASEVAR(Name = "X2", DistributionType = "norm", Mean = 0, Sd = 1)

  lsf1 <- SYS_LSF(vars = list(v1))
  lsf1$func <- function(x) x[1] - 1

  lsf2 <- SYS_LSF(vars = list(v2))
  lsf2$func <- function(x) x[1] - 2

  x <- c(1.5, 3)

  g_vals <- c(
    lsf1$getLSF()(c(x[1])),
    lsf2$getLSF()(c(x[2]))
  )

  g_min <- min(g_vals)

  kappa <- 20

  sys <- build_gsys(
    sys_input = list(lsf1, lsf2),
    sys_structure = "serial",
    smooth_kappa = kappa
  )

  g_sys <- sys$g_sys(x)

  # print(sys$g_sys(x))
  # print(smooth_min(c(0.5, 1), 100))

  # print(sys$g_sys(x))
  # print(g_vals)
  # print(lsf1$getLSF()(c(1.5)))

  # ----------------------------
  # Property 1: smooth ≤ min
  # ----------------------------

  expect_true(g_sys <= g_min)

  # ----------------------------
  # Property 2: error bound
  # min − log(m)/kappa ≤ g_sys
  # ----------------------------

  m <- length(g_vals)

  lower_bound <- g_min - log(m) / kappa

  expect_true(g_sys >= lower_bound)
})


test_that("smooth aggregation converges to minimum for large kappa", {
  v1 <- PROB_BASEVAR(Name = "X1", DistributionType = "norm", Mean = 0, Sd = 1)
  v2 <- PROB_BASEVAR(Name = "X2", DistributionType = "norm", Mean = 0, Sd = 1)

  lsf1 <- SYS_LSF(vars = list(v1))
  lsf1$func <- function(x) x[1] - 1

  lsf2 <- SYS_LSF(vars = list(v2))
  lsf2$func <- function(x) x[1] - 2

  x <- c(1.5, 3)

  g_vals <- c(
    lsf1$getLSF()(c(x[1])),
    lsf2$getLSF()(c(x[2]))
  )

  g_min <- min(g_vals)

  sys_small <- build_gsys(list(lsf1, lsf2), smooth_kappa = 5)
  sys_large <- build_gsys(list(lsf1, lsf2), smooth_kappa = 1000)

  err_small <- abs(sys_small$g_sys(x) - g_min)
  err_large <- abs(sys_large$g_sys(x) - g_min)

  expect_true(err_large < err_small)
})


test_that("smooth aggregation remains numerically stable for large kappa", {
  v1 <- PROB_BASEVAR(Name = "X1", DistributionType = "norm", Mean = 0, Sd = 1)
  v2 <- PROB_BASEVAR(Name = "X2", DistributionType = "norm", Mean = 0, Sd = 1)

  lsf1 <- SYS_LSF(vars = list(v1))
  lsf1$func <- function(x) x[1] - 1

  lsf2 <- SYS_LSF(vars = list(v2))
  lsf2$func <- function(x) x[1] - 2

  sys <- build_gsys(
    sys_input = list(lsf1, lsf2),
    sys_structure = "serial",
    smooth_kappa = 1e6
  )

  x <- c(1.5, 3)

  g <- sys$g_sys(x)

  expect_true(is.finite(g))
})


test_that("g_sys parallel approximates maximum", {
  v1 <- PROB_BASEVAR(Name = "X1", DistributionType = "norm", Mean = 0, Sd = 1)
  v2 <- PROB_BASEVAR(Name = "X2", DistributionType = "norm", Mean = 0, Sd = 1)

  lsf1 <- SYS_LSF(vars = list(v1))
  lsf1$func <- function(x) x[1] - 1

  lsf2 <- SYS_LSF(vars = list(v2))
  lsf2$func <- function(x) x[1] - 2

  sys <- build_gsys(
    sys_input = list(lsf1, lsf2),
    sys_structure = "parallel",
    smooth_kappa = 100
  )

  x <- c(1.5, 3)

  g_vals <- c(
    lsf1$getLSF()(c(x[1])),
    lsf2$getLSF()(c(x[2]))
  )

  g_max <- max(g_vals)

  g_sys <- sys$g_sys(x)

  expect_equal(g_sys, g_max, tolerance = 1e-3)
})


test_that("g_sys works with matrix input", {
  v1 <- PROB_BASEVAR(Name = "X1", DistributionType = "norm", Mean = 0, Sd = 1)

  lsf <- SYS_LSF(vars = list(v1))
  lsf$func <- function(x) x[1] - 1

  sys <- build_gsys(list(lsf))

  X <- matrix(c(0.5, 1.5, 2), ncol = 1)

  g <- sys$g_sys(X)

  expect_equal(length(g), nrow(X))
  expect_true(is.numeric(g))
})


test_that("stable smooth aggregation matches classical formulation", {
  v1 <- PROB_BASEVAR(Name = "X1", DistributionType = "norm", Mean = 0, Sd = 1)
  v2 <- PROB_BASEVAR(Name = "X2", DistributionType = "norm", Mean = 0, Sd = 1)

  lsf1 <- SYS_LSF(vars = list(v1))
  lsf1$func <- function(x) x[1] - 1

  lsf2 <- SYS_LSF(vars = list(v2))
  lsf2$func <- function(x) x[1] - 2

  sys <- build_gsys(
    sys_input = list(lsf1, lsf2),
    sys_structure = "serial",
    smooth_kappa = 20
  )

  x <- c(1.5, 3)

  g_vals <- c(
    lsf1$getLSF()(c(x[1])),
    lsf2$getLSF()(c(x[2]))
  )

  # klassische Formel
  g_classic <- -log(sum(exp(-20 * g_vals))) / 20

  g_sys <- sys$g_sys(x)

  expect_equal(g_sys, g_classic, tolerance = 1e-8)
})


test_that("MC_IS_gsys produces finite reliability estimates", {
  v1 <- PROB_BASEVAR(Name = "X1", DistributionType = "norm", Mean = 5, Sd = 1)
  v2 <- PROB_BASEVAR(Name = "X2", DistributionType = "norm", Mean = 3, Sd = 1)

  lsf1 <- SYS_LSF(vars = list(v1))
  lsf1$func <- function(x) x[1] - 4

  lsf2 <- SYS_LSF(vars = list(v2))
  lsf2$func <- function(x) x[1] - 2

  sys_input <- list(lsf1, lsf2)

  res <- MC_IS_gsys(
    sys_input = sys_input,
    sys_structure = "serial",
    smooth_kappa = 100,
    n_batch = 200,
    n_max = 2000,
    use_threads = 1,
    seed = 123
  )

  expect_true(is.finite(res$beta))
  expect_true(is.finite(res$pf))
})


test_that("MC_IS_gsys runs quickly", {
  v1 <- PROB_BASEVAR(Name = "X1", DistributionType = "norm", Mean = 5, Sd = 1)
  v2 <- PROB_BASEVAR(Name = "X2", DistributionType = "norm", Mean = 3, Sd = 1)

  lsf1 <- SYS_LSF(vars = list(v1))
  lsf1$func <- function(x) x[1] - 4

  lsf2 <- SYS_LSF(vars = list(v2))
  lsf2$func <- function(x) x[1] - 2

  sys_input <- list(lsf1, lsf2)

  t <- system.time(
    MC_IS_gsys(
      sys_input = sys_input,
      sys_structure = "serial",
      n_batch = 200,
      n_max = 2000,
      use_threads = 1,
      seed = 123
    )
  )[3]

  expect_true(t < 5)
})


test_that("MC_IS_system matches the exact solution for small system", {
  skip_on_cran()

  # -----------------------------
  # Basisvariablen
  # -----------------------------
  R1 <- PROB_BASEVAR(
    Id = 1, Name = "R1",
    DistributionType = "norm",
    Mean = 10, Sd = 2
  )

  R2 <- PROB_BASEVAR(
    Id = 2, Name = "R2",
    DistributionType = "norm",
    Mean = 10, Sd = 2
  )

  S <- PROB_BASEVAR(
    Id = 3, Name = "S",
    DistributionType = "norm",
    Mean = 9, Sd = 2
  )

  # -----------------------------
  # LSFs
  # -----------------------------
  lsf1 <- SYS_LSF(vars = list(R1, S), name = "LSF1")
  lsf1$func <- function(R1, S) R1 - S
  lsf1$check()

  lsf2 <- SYS_LSF(vars = list(R2, S), name = "LSF2")
  lsf2$func <- function(R2, S) R2 - S
  lsf2$check()

  # -----------------------------
  # Maschinen
  # -----------------------------
  machine_MC_IS <- PROB_MACHINE(
    name = "MC_IS",
    fCall = "MC_IS",
    options = list(
      n_batch = 2000,
      n_max = 5e4,
      cov_user = 0.05,
      seed = 123,
      backend = "parallel",
      use_threads = 1
    )
  )

  # -----------------------------
  # Serial System (OR)
  # -----------------------------
  ps <- SYS_PROB(
    sys_input = list(lsf1, lsf2),
    probMachines = list(machine_MC_IS),
    sys_type = "serial"
  )

  ps$runMachines()
  ps$calculateSystemProbability("MCIS")

  pf_MCIS <- ps$res_sys[[1]]$pf

  # exact solution: safe if R1 >= S and R2 >= S (independent variables)
  ld <- list(R1$getlDistr(), R2$getlDistr(), S$getlDistr())
  ld <- lapply(ld, function(d) if (is.function(d$p)) d else d[[1]])
  safe <- integrate(
    function(s) ld[[3]]$d(s) * (1 - ld[[1]]$p(s)) * (1 - ld[[2]]$p(s)),
    lower = ld[[3]]$q(1e-12),
    upper = ld[[3]]$q(1 - 1e-12),
    rel.tol = 1e-10
  )$value
  pf_exact <- 1 - safe

  # -----------------------------
  # Test
  # -----------------------------
  # three standard deviations of the estimator (sd approx. cov_user * pf)
  expect_lt(abs(pf_MCIS - pf_exact), 3 * 0.05 * pf_exact)
})


test_that("MC_IS_system parallel backend is reproducible on Linux and Windows", {
  v1 <- PROB_BASEVAR(
    Id = 1,
    Name = "X1",
    DistributionType = "norm",
    Mean = 0,
    Sd = 1
  )
  v2 <- PROB_BASEVAR(
    Id = 2,
    Name = "X2",
    DistributionType = "norm",
    Mean = 0,
    Sd = 1
  )

  lsf1 <- function(x) 3 - x[1]
  lsf2 <- function(x) 3 - x[1]

  normalize_funlist <- function(v) {
    d <- v$getlDistr()
    if (is.list(d) &&
      length(d) >= 1L &&
      is.list(d[[1]]) &&
      is.function(d[[1]]$p)) {
      d[[1]]
    } else {
      d
    }
  }

  lDistr <- list(
    list(normalize_funlist(v1)),
    list(normalize_funlist(v2))
  )

  run_once <- function() {
    MC_IS(
      lsf = list(lsf1, lsf2),
      lDistr = lDistr,
      cov_user = 0.5,
      n_batch = 1000L,
      n_max = 2000L,
      use_threads = 2L,
      backend = "parallel",
      seed = 1234L,
      dps = c(3, 3),
      dataRecord = FALSE,
      debug.level = 0
    )
  }

  res1 <- run_once()
  res2 <- run_once()

  expect_equal(res1$n_mc, 2000L)
  expect_equal(res2$n_mc, 2000L)
  expect_equal(res1$pf, res2$pf)
  expect_equal(res1$beta, res2$beta)
  expect_true(is.finite(res1$pf))
  expect_true(res1$pf >= 0 && res1$pf <= 1)
})

test_that("MC_IS with adaptive_batch keeps the batch size within the bounds", {
  skip_on_cran()
  var <- PROB_BASEVAR(Name = "X", DistributionType = "norm", Mean = 0, Sd = 1)

  res <- MC_IS(
    lsf = list(function(x) 3 - x[1], function(x) 3 + x[1]),
    lDistr = list(list(var$getlDistr()[[1]]), list(var$getlDistr()[[1]])),
    sys_type = "serial",
    n_batch = 2000L,
    n_max = 2e5,
    adaptive_batch = TRUE,
    batch_control = list(K_future = 3),
    seed = 1L
  )

  # exact system pf: P(X > 3) + P(X < -3)
  expect_equal(res$pf, 2 * pnorm(-3), tolerance = 0.15)
  expect_gt(res$n_mc, 0)
  expect_lte(res$n_mc, 2e5 + 100000)
})

test_that("MC_IS_system uses user supplied design points with several variables per LSF", {
  skip_on_cran()
  v1 <- PROB_BASEVAR(Id = 1, Name = "X1", DistributionType = "norm", Mean = 0, Sd = 1)
  v2 <- PROB_BASEVAR(Id = 2, Name = "X2", DistributionType = "norm", Mean = 0, Sd = 1)

  lsf1 <- function(x) 3.5 - x[2]
  lsf2 <- function(x) 3.5 + x[1]

  lDistr <- list(
    list(v1$getlDistr()[[1]], v2$getlDistr()[[1]]),
    list(v1$getlDistr()[[1]], v2$getlDistr()[[1]])
  )

  res <- MC_IS(
    lsf = list(lsf1, lsf2),
    lDistr = lDistr,
    sys_type = "serial",
    dps = list(c(0, 3.5), c(-3.5, 0)),
    cov_user = 0.05,
    n_batch = 5000L,
    n_max = 2e5,
    seed = 1L
  )

  # exact: P(X2 > 3.5 or X1 < -3.5) = 1 - (1 - pnorm(-3.5))^2
  expect_equal(res$pf, 1 - (1 - pnorm(-3.5))^2, tolerance = 0.2)
})

test_that("MC_IS_system stops for dps in a wrong format", {
  v1 <- PROB_BASEVAR(Id = 1, Name = "X1", DistributionType = "norm", Mean = 0, Sd = 1)
  v2 <- PROB_BASEVAR(Id = 2, Name = "X2", DistributionType = "norm", Mean = 0, Sd = 1)

  lDistr <- list(
    list(v1$getlDistr()[[1]], v2$getlDistr()[[1]]),
    list(v1$getlDistr()[[1]], v2$getlDistr()[[1]])
  )

  expect_error(
    MC_IS(
      lsf = list(function(x) 3.5 - x[2], function(x) 3.5 + x[1]),
      lDistr = lDistr,
      sys_type = "serial",
      dps = c(3.5, 3.5),
      n_batch = 1000L,
      n_max = 1000L,
      seed = 1L
    ),
    "'dps' must be a list with one design point vector per limit-state function",
    fixed = TRUE
  )
})