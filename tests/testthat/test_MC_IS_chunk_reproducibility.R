test_that("MC_IS results are bit-identical across Future worker counts", {
    x_var <- PROB_BASEVAR(
        Id = 1,
        Name = "X",
        DistributionType = "norm",
        Mean = 0,
        Sd = 1
    )
    marginal <- x_var$getlDistr()[[1]]

    lsf <- function(x) {
        if (is.matrix(x)) {
            4 - x[, 1]
        } else {
            4 - x[1]
        }
    }

    run_case <- function() {
        MC_IS(
            lsf = lsf,
            lDistr = list(marginal),
            backend = "future",
            # Do not supply use_threads: the active Future plan controls workers.
            dps = 4,
            seed = 12345L,
            n_batch = 10000L,
            n_max = 10000L,
            cov_user = 0,
            dataRecord = FALSE,
            debug.level = 0
        )
    }

    res_seq <- with_test_future_plan(
        future::sequential,
        run_case()
    )

    res_multi <- with_test_future_plan(
        future::multisession,
        run_case(),
        workers = 2L
    )

    fields_to_compare <- c(
        "method", "beta", "pf", "log_pf", "var",
        "cov_mc", "ESS", "ESS_ratio", "n_mc", "design_points"
    )

    expect_identical(
        res_seq[fields_to_compare],
        res_multi[fields_to_compare]
    )
})


# test_that("MC_IS results are bit-identical across Future worker counts and chunk boundaries", {
#     # Batch sizes straddle the proposed fixed chunk size of 1000.
#     cases <- list(
#         list(name = "1var_below_chunk", weights = 1, threshold = 4, n_batch = 999L),
#         list(name = "1var_at_chunk", weights = 1, threshold = 4, n_batch = 1000L),
#         list(name = "1var_above_chunk", weights = 1, threshold = 4, n_batch = 1001L),
#         list(name = "1var_multiple_chunks", weights = 1, threshold = 4, n_batch = 4097L),
#         list(
#             name = "4var_multiple_chunks",
#             weights = c(1, -0.5, 0.8, 1.2),
#             threshold = 4,
#             n_batch = 4097L
#         )
#     )

#     seeds <- c(123L, 987L)

#     run_case <- function(case, seed) {
#         weights <- case$weights
#         n_vars <- length(weights)

#         vars <- lapply(seq_len(n_vars), function(i) {
#             PROB_BASEVAR(
#                 Id = i,
#                 Name = paste0("X", i),
#                 DistributionType = "norm",
#                 Mean = 0,
#                 Sd = 1
#             )
#         })

#         lDistr <- lapply(vars, function(v) v$getlDistr()[[1]])

#         # Linear LSF, vectorized for a matrix of samples.
#         lsf <- function(x) {
#             if (is.matrix(x)) {
#                 case$threshold - drop(x %*% weights)
#             } else {
#                 case$threshold - sum(x * weights)
#             }
#         }

#         # Design point on the linear limit-state surface.
#         dps <- case$threshold * weights / sum(weights^2)

#         MC_IS(
#             lsf = lsf,
#             lDistr = lDistr,
#             backend = "future",
#             # The active Future plan, not use_threads, controls the worker count.
#             dps = dps,
#             seed = seed,
#             n_batch = case$n_batch,
#             n_max = case$n_batch,
#             cov_user = 0,
#             dataRecord = FALSE,
#             debug.level = 0
#         )
#     }

#     run_suite <- function() {
#         results <- list()

#         for (case in cases) {
#             for (seed in seeds) {
#                 key <- paste(case$name, seed, sep = "_")
#                 results[[key]] <- run_case(case, seed)
#             }
#         }

#         results
#     }

#     fields_to_compare <- c(
#         "method",
#         "beta",
#         "pf",
#         "log_pf",
#         "var",
#         "cov_mc",
#         "ESS",
#         "ESS_ratio",
#         "n_mc",
#         "design_points"
#     )

#     reference <- suppressMessages(
#         with_test_future_plan(
#             future::sequential,
#             run_suite()
#         )
#     )

#     for (workers in c(2L, 4L)) {
#         actual <- suppressMessages(
#             with_test_future_plan(
#                 future::multisession,
#                 run_suite(),
#                 workers = workers
#             )
#         )

#         expect_identical(names(actual), names(reference))

#         for (key in names(reference)) {
#             expect_identical(
#                 actual[[key]][fields_to_compare],
#                 reference[[key]][fields_to_compare],
#                 info = sprintf(
#                     "Results differ between sequential and multisession (%d workers) for %s",
#                     workers,
#                     key
#                 )
#             )
#         }
#     }
# })
