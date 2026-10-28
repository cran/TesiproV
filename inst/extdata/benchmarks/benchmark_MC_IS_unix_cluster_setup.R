# ======================================================================
# benchmark_MC_IS_unix_cluster_setup.R
#
# Benchmark for MC_IS startup/setup overhead on Unix-like systems.
#
# Purpose:
# Compare runtime before and after removing the unused parallel cluster
# creation on Unix for backend = "parallel".
#
# Run from the package root, e.g.:
#   devtools::load_all()
#   source("inst/extdata/benchmarks/benchmark_MC_IS_unix_cluster_setup.R")
# ======================================================================

options(digits = 14)

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

# Change manually before / after the refactoring.
# This value is used in the output-file name.
benchmark_version <- "after"

# The benchmark targets the Unix mclapply path.
if (.Platform$OS.type == "windows") {
    stop(
        "This benchmark is intended for Unix-like systems because it measures ",
        "the unnecessary cluster creation before mclapply().",
        call. = FALSE
    )
}

backend <- "parallel"

# Number of workers tested.
worker_list <- c(1L, 4L, 8L)

# Number of measured repetitions per worker configuration.
n_runs <- 10L

# One unrecorded warm-up per worker configuration.
warmup_runs <- 1L

# Use exactly one MC batch in every run.
# In MC_IS_single(), n_batch is the total batch size and is split among
# workers. Setting n_max = n_batch therefore limits every run to one batch.
n_batch_setup <- 2000L
n_max_setup <- 2000L

# Fixed seed: beta and number of samples should be reproducible.
seed_value <- 12345L

# A permissive CoV target ensures that the benchmark does not need
# additional batches because of a strict convergence requirement.
cov_user_setup <- 1.0

# ----------------------------------------------------------------------
# Output directory
# ----------------------------------------------------------------------

output_dir <- file.path(
    "inst",
    "extdata",
    "benchmarks",
    "results"
)

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

output_csv <- file.path(
    output_dir,
    paste0(
        "mc_is_unix_cluster_setup_",
        benchmark_version,
        ".csv"
    )
)

output_summary_csv <- file.path(
    output_dir,
    paste0(
        "mc_is_unix_cluster_setup_",
        benchmark_version,
        "_summary.csv"
    )
)

# ----------------------------------------------------------------------
# Base variables
# ----------------------------------------------------------------------

Z <- PROB_BASEVAR(
    Id = 1,
    Name = "Z",
    DistributionType = "norm",
    Mean = 100,
    Cov = 0.04
)

Fy <- PROB_BASEVAR(
    Id = 2,
    Name = "Fy",
    DistributionType = "lnorm",
    Mean = 40,
    Cov = 0.10
)

M <- PROB_BASEVAR(
    Id = 3,
    Name = "M",
    DistributionType = "gumbel",
    Package = "evd",
    Mean = 2000,
    Cov = 0.10
)

# ----------------------------------------------------------------------
# Fast vectorizable limit-state function
# ----------------------------------------------------------------------

# Explicit scalar arguments are used so that SYS_LSF$getLSF() can also
# evaluate the LSF efficiently for matrix input in MC_IS.
lsf_setup_fun <- function(Z, Fy, M) {
    Z * Fy - M
}

lsf_setup <- SYS_LSF(
    vars = list(Z, Fy, M),
    name = "MC_IS Unix cluster setup benchmark"
)

lsf_setup$func <- lsf_setup_fun
lsf_setup$check()

# ----------------------------------------------------------------------
# Helper: execute one short MC_IS calculation
# ----------------------------------------------------------------------

run_one_setup_case <- function(workers, seed) {
    machine_mcis <- PROB_MACHINE(
        name = "MC IS setup benchmark",
        fCall = "MC_IS",
        options = list(
            n_batch = n_batch_setup,
            n_max = n_max_setup,
            cov_user = cov_user_setup,
            use_threads = workers,
            backend = backend,
            seed = seed,
            dataRecord = FALSE,
            densityType = "norm"
        )
    )

    ps <- SYS_PROB(
        sys_input = list(lsf_setup),
        probMachines = list(machine_mcis),
        sys_type = "serial",
        debug.level = 0
    )

    tic <- proc.time()
    ps$runMachines()
    elapsed <- (proc.time() - tic)[["elapsed"]]

    res <- ps$res_single[[1]][[1]]

    list(
        elapsed_sec = elapsed,
        n_mc = res$n_mc,
        beta = res$beta,
        pf = res$pf
    )
}

# ----------------------------------------------------------------------
# Benchmark
# ----------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("MC_IS Unix cluster setup benchmark\n")
cat("============================================================\n")
cat(sprintf("Version label: %s\n", benchmark_version))
cat(sprintf("Operating system: %s\n", Sys.info()[["sysname"]]))
cat(sprintf("Backend: %s\n", backend))
cat(sprintf("Workers: %s\n", paste(worker_list, collapse = ", ")))
cat(sprintf("Warm-up runs: %d\n", warmup_runs))
cat(sprintf("Measured runs per configuration: %d\n", n_runs))
cat(sprintf("n_batch = %d, n_max = %d\n", n_batch_setup, n_max_setup))
cat("============================================================\n")

results_list <- list()
result_id <- 0L

for (workers in worker_list) {
    cat("\n")
    cat("------------------------------------------------------------\n")
    cat(sprintf("Workers: %d\n", workers))
    cat("------------------------------------------------------------\n")

    # Warm-up runs are deliberately excluded from recorded timings.
    if (warmup_runs > 0L) {
        cat("Warm-up ...\n")

        for (warmup in seq_len(warmup_runs)) {
            invisible(
                run_one_setup_case(
                    workers = workers,
                    seed = seed_value
                )
            )
        }
    }

    cat("Measured runs ...\n")

    for (run in seq_len(n_runs)) {
        # Identical seed for every run: results should be identical.
        out <- run_one_setup_case(
            workers = workers,
            seed = seed_value
        )

        result_id <- result_id + 1L

        results_list[[result_id]] <- data.frame(
            version = benchmark_version,
            workers = workers,
            run = run,
            elapsed_sec = out$elapsed_sec,
            n_mc = out$n_mc,
            beta = out$beta,
            pf = out$pf
        )

        cat(sprintf(
            "  Run %2d/%2d: %.4f s | n_mc = %d | beta = %.10f\n",
            run,
            n_runs,
            out$elapsed_sec,
            out$n_mc,
            out$beta
        ))
    }
}

results <- do.call(rbind, results_list)

# ----------------------------------------------------------------------
# Consistency checks
# ----------------------------------------------------------------------

# All runs should use exactly one batch.
if (!all(results$n_mc == n_batch_setup)) {
    warning(
        "At least one run used a different number of samples than n_batch_setup. ",
        "Check the stopping behaviour before comparing timings."
    )
}

# For the fixed seed, results should be reproducible within each worker count.
for (workers in worker_list) {
    subset_results <- results[results$workers == workers, , drop = FALSE]

    if (length(unique(subset_results$n_mc)) != 1L) {
        warning(
            sprintf(
                "n_mc differs between repetitions for workers = %d.",
                workers
            )
        )
    }

    if (length(unique(subset_results$beta)) != 1L) {
        warning(
            sprintf(
                "beta differs between repetitions for workers = %d.",
                workers
            )
        )
    }
}

# ----------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------

summary_list <- lapply(worker_list, function(workers) {
    x <- results$elapsed_sec[results$workers == workers]

    data.frame(
        version = benchmark_version,
        workers = workers,
        n_runs = length(x),
        mean_time_sec = mean(x),
        median_time_sec = median(x),
        sd_time_sec = stats::sd(x),
        relative_sd = stats::sd(x) / mean(x),
        min_time_sec = min(x),
        max_time_sec = max(x),
        mean_n_mc = mean(results$n_mc[results$workers == workers]),
        mean_beta = mean(results$beta[results$workers == workers])
    )
})

summary_results <- do.call(rbind, summary_list)

# Speed-up relative to one worker.
baseline_mean <- summary_results$mean_time_sec[
    summary_results$workers == 1L
]

baseline_median <- summary_results$median_time_sec[
    summary_results$workers == 1L
]

summary_results$mean_speedup <- baseline_mean / summary_results$mean_time_sec
summary_results$median_speedup <- baseline_median / summary_results$median_time_sec

summary_results$mean_efficiency <-
    summary_results$mean_speedup / summary_results$workers

summary_results$median_efficiency <-
    summary_results$median_speedup / summary_results$workers

# ----------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("Raw benchmark results\n")
cat("============================================================\n")
print(results)

cat("\n")
cat("============================================================\n")
cat("Benchmark summary\n")
cat("============================================================\n")
print(summary_results)

utils::write.csv(
    results,
    file = output_csv,
    row.names = FALSE
)

utils::write.csv(
    summary_results,
    file = output_summary_csv,
    row.names = FALSE
)

cat("\n")
cat("Results written to:\n")
cat(sprintf("  %s\n", output_csv))
cat(sprintf("  %s\n", output_summary_csv))

invisible(
    list(
        raw = results,
        summary = summary_results
    )
)
