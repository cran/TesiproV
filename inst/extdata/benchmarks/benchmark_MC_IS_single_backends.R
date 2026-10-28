# ============================================================
# Benchmark MC_IS_single:
# parallel vs future::multisession vs future::multicore
#
# Run from the TesiproV package root:
#   source("inst/extdata/benchmarks/benchmark_MC_IS_single_backends.R")
#
# Set benchmark_version to "before" or "after" for separate CSV files.
# ============================================================

# ------------------------------------------------------------
# Benchmark settings
# ------------------------------------------------------------

benchmark_version <- "after"

n_vars <- 16L
n_batch <- 250000L
n_max <- n_batch

n_repeats <- 10L
warmup_runs <- 1L
seed_value <- 12345L

# Requested parallel worker counts.
worker_counts <- c(1L, 2L, 4L, 8L)

# Workloads:
#   "cheap"       exposes transform / chunk / dispatch overhead
#   "compute_heavy" gives the LSF more work per sample
workloads <- c("cheap", "compute_heavy")

output_dir <- "inst/extdata/benchmarks/results"

# ------------------------------------------------------------
# Validation and package setup
# ------------------------------------------------------------

if (!file.exists("DESCRIPTION")) {
    stop("Run this script from the TesiproV package root.", call. = FALSE)
}

for (pkg in c("devtools", "future")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(
            sprintf("Package '%s' is required to run this benchmark.", pkg),
            call. = FALSE
        )
    }
}

devtools::load_all(".", quiet = TRUE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw_file <- file.path(
    output_dir,
    sprintf("mc_is_single_backends_%s.csv", benchmark_version)
)

summary_file <- file.path(
    output_dir,
    sprintf("mc_is_single_backends_%s_summary.csv", benchmark_version)
)

# ------------------------------------------------------------
# Base variables and marginal distribution funlists
# ------------------------------------------------------------

vars <- lapply(seq_len(n_vars), function(i) {
    PROB_BASEVAR(
        Id = i,
        Name = paste0("X", i),
        DistributionType = "norm",
        Mean = 0,
        Sd = 1
    )
})

# MC_IS accepts the marginal funlists returned by getlDistr()[[1]].
lDistr <- lapply(vars, function(v) v$getlDistr()[[1]])

# The event for the cheap workload is X1 >= 4.
# For standard normal variables, the physical-space design point is
# also the standard-normal-space point.
dps <- c(4, rep(0, n_vars - 1L))

# ------------------------------------------------------------
# Vectorized LSFs
# ------------------------------------------------------------

make_lsf <- function(workload) {
    force(workload)

    function(x) {
        if (is.matrix(x)) {
            z <- x[, 1L]

            if (workload == "cheap") {
                return(4 - z)
            }

            # A vectorized computational workload with a small perturbation
            # of the cheap limit state. The proposal remains centred near
            # the failure boundary.
            fy <- x[, 2L]
        } else {
            z <- x[1L]

            if (workload == "cheap") {
                return(4 - z)
            }

            fy <- x[2L]
        }

        val <- numeric(length(z))

        for (k in seq_len(100L)) {
            val <- val +
                sin(z * 0.001 * k) *
                    cos(fy * 0.002 * k)
        }

        4 - z + 1e-6 * val
    }
}

# ------------------------------------------------------------
# Scenario construction
# ------------------------------------------------------------

scenarios <- list()

add_scenario <- function(name,
                         backend,
                         plan_type = NA_character_,
                         workers = NA_integer_) {
    scenarios[[length(scenarios) + 1L]] <<- list(
        name = name,
        backend = backend,
        plan_type = plan_type,
        workers = as.integer(workers)
    )
}

# Base parallel backend: worker count is passed via use_threads.
for (w in worker_counts) {
    add_scenario(
        name = paste0("parallel_", w),
        backend = "parallel",
        workers = w
    )
}

# Future sequential plan.
add_scenario(
    name = "future_sequential",
    backend = "future",
    plan_type = "sequential",
    workers = 1L
)

# Future multisession plans.
for (w in worker_counts) {
    add_scenario(
        name = paste0("future_multisession_", w),
        backend = "future",
        plan_type = "multisession",
        workers = w
    )
}

# Future multicore is only included when supported by the current session.
multicore_supported <- isTRUE(
    tryCatch(
        future::supportsMulticore(),
        error = function(e) FALSE
    )
)

if (multicore_supported) {
    for (w in worker_counts[worker_counts >= 2L]) {
        add_scenario(
            name = paste0("future_multicore_", w),
            backend = "future",
            plan_type = "multicore",
            workers = w
        )
    }
} else {
    message(
        "Skipping future::multicore: it is not supported in this R session."
    )
}

# ------------------------------------------------------------
# Future-plan helper
# ------------------------------------------------------------

with_scenario_plan <- function(scenario, code) {
    if (scenario$backend != "future") {
        return(force(code))
    }

    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)

    if (scenario$plan_type == "sequential") {
        future::plan(future::sequential)
    } else if (scenario$plan_type == "multisession") {
        future::plan(
            future::multisession,
            workers = scenario$workers
        )
    } else if (scenario$plan_type == "multicore") {
        future::plan(
            future::multicore,
            workers = scenario$workers
        )
    } else {
        stop("Unknown Future plan type: ", scenario$plan_type)
    }

    force(code)
}

# ------------------------------------------------------------
# One fixed-budget MC_IS run
# ------------------------------------------------------------

run_one <- function(scenario, workload) {
    lsf <- make_lsf(workload)

    args <- list(
        lsf = lsf,
        lDistr = lDistr,
        backend = scenario$backend,
        dps = dps,
        seed = seed_value,
        n_batch = n_batch,
        n_max = n_max,
        cov_user = 0,
        dataRecord = FALSE,
        debug.level = 0
    )

    # Only the "parallel" backend uses use_threads.
    # For "future", the active Future plan defines the worker count.
    if (scenario$backend == "parallel") {
        args$use_threads <- scenario$workers
    }

    timing <- system.time({
        result <- do.call(MC_IS, args)
    })

    list(
        elapsed_sec = unname(timing[["elapsed"]]),
        result = result
    )
}

# ------------------------------------------------------------
# Run scenarios
# ------------------------------------------------------------

set.seed(2026)
scenario_order <- sample(seq_along(scenarios))

rows <- list()
row_id <- 0L

cat("MC_IS_single backend benchmark\n")
cat("Version:", benchmark_version, "\n")
cat("OS:", Sys.info()[["sysname"]], "\n")
cat("R:", R.version.string, "\n")
cat("Variables:", n_vars, "\n")
cat("Samples per run:", n_batch, "\n")
cat("Warm-up runs:", warmup_runs, "\n")
cat("Measured runs per scenario/workload:", n_repeats, "\n")
cat("Workloads:", paste(workloads, collapse = ", "), "\n")

for (scenario_idx in scenario_order) {
    scenario <- scenarios[[scenario_idx]]

    for (workload in workloads) {
        cat("\n------------------------------------------------------------\n")
        cat(
            "Scenario:", scenario$name,
            "| Workload:", workload,
            "| Backend:", scenario$backend,
            "| Requested workers:", scenario$workers,
            "\n"
        )
        cat("------------------------------------------------------------\n")

        scenario_results <- with_scenario_plan(
            scenario,
            {
                actual_workers <- if (scenario$backend == "future") {
                    as.integer(future::nbrOfWorkers())
                } else {
                    scenario$workers
                }

                cat("Actual workers:", actual_workers, "\n")

                # Warm-up runs are excluded from the recorded timings.
                if (warmup_runs > 0L) {
                    for (i in seq_len(warmup_runs)) {
                        invisible(run_one(scenario, workload))
                    }
                }

                measured <- vector("list", n_repeats)

                for (i in seq_len(n_repeats)) {
                    measured[[i]] <- run_one(scenario, workload)
                    r <- measured[[i]]$result

                    if (r$n_mc != n_batch) {
                        stop(
                            sprintf(
                                "Expected n_mc = %d but got %d in %s.",
                                n_batch, r$n_mc, scenario$name
                            ),
                            call. = FALSE
                        )
                    }

                    cat(sprintf(
                        paste0(
                            "Run %d/%d: %.3f s | n_mc = %d",
                            " | beta = %.12f | pf = %.12e\n"
                        ),
                        i,
                        n_repeats,
                        measured[[i]]$elapsed_sec,
                        r$n_mc,
                        r$beta,
                        r$pf
                    ))
                }

                list(
                    actual_workers = actual_workers,
                    measured = measured
                )
            }
        )

        ref <- scenario_results$measured[[1L]]$result

        for (i in seq_along(scenario_results$measured)) {
            out <- scenario_results$measured[[i]]
            r <- out$result

            row_id <- row_id + 1L

            rows[[row_id]] <- data.frame(
                version = benchmark_version,
                os = Sys.info()[["sysname"]],
                r_version = R.version.string,
                scenario = scenario$name,
                backend = scenario$backend,
                future_plan = scenario$plan_type,
                requested_workers = scenario$workers,
                actual_workers = scenario_results$actual_workers,
                workload = workload,
                run = i,
                elapsed_sec = out$elapsed_sec,
                samples_per_sec = r$n_mc / out$elapsed_sec,
                n_mc = r$n_mc,
                beta = r$beta,
                pf = r$pf,
                log_pf = r$log_pf,
                beta_diff_from_first_run = abs(r$beta - ref$beta),
                pf_diff_from_first_run = abs(r$pf - ref$pf),
                stringsAsFactors = FALSE
            )
        }
    }
}

benchmark_results <- do.call(rbind, rows)

# ------------------------------------------------------------
# Summaries and speed-up relative to Future sequential
# ------------------------------------------------------------

group_key <- interaction(
    benchmark_results$workload,
    benchmark_results$scenario,
    drop = TRUE,
    lex.order = TRUE
)

summary_results <- do.call(
    rbind,
    lapply(split(benchmark_results, group_key), function(d) {
        data.frame(
            version = d$version[1L],
            workload = d$workload[1L],
            scenario = d$scenario[1L],
            backend = d$backend[1L],
            future_plan = d$future_plan[1L],
            requested_workers = d$requested_workers[1L],
            actual_workers = d$actual_workers[1L],
            n_repeats = nrow(d),
            median_time_sec = median(d$elapsed_sec),
            mean_time_sec = mean(d$elapsed_sec),
            sd_time_sec = stats::sd(d$elapsed_sec),
            median_samples_per_sec = median(d$samples_per_sec),
            beta = d$beta[1L],
            pf = d$pf[1L],
            max_beta_diff = max(d$beta_diff_from_first_run),
            max_pf_diff = max(d$pf_diff_from_first_run),
            stringsAsFactors = FALSE
        )
    })
)

summary_results$speedup_vs_future_sequential <- NA_real_

for (workload in workloads) {
    idx_seq <- which(
        summary_results$workload == workload &
            summary_results$scenario == "future_sequential"
    )

    if (length(idx_seq) == 1L) {
        seq_median <- summary_results$median_time_sec[idx_seq]
        idx_workload <- summary_results$workload == workload

        summary_results$speedup_vs_future_sequential[idx_workload] <-
            seq_median / summary_results$median_time_sec[idx_workload]
    }
}

# ------------------------------------------------------------
# Save and print results
# ------------------------------------------------------------

utils::write.csv(
    benchmark_results,
    raw_file,
    row.names = FALSE
)

utils::write.csv(
    summary_results,
    summary_file,
    row.names = FALSE
)

cat("\n================ Raw results ================\n")
print(benchmark_results, row.names = FALSE)

cat("\n================ Summary ====================\n")
print(summary_results, row.names = FALSE)

cat("\nRaw results written to:\n", normalizePath(raw_file), "\n")
cat("Summary written to:\n", normalizePath(summary_file), "\n")
