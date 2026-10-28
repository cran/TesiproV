# Benchmark MC_IS Future memory and runtime.
#
# Run from the TesiproV package root:
#   source("inst/extdata/benchmarks/benchmark_MC_IS_future_memory.R")
#
# Run once before and once after the chunking change.
# Change benchmark_version to keep both result files.

# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------

benchmark_version <- "after" # set to "after" for the second run

# Adjust these to a representative workload. This setting creates
# 32 MB for u_mat alone: 250,000 rows * 16 doubles * 8 bytes.
n_batch <- 250000L
n_max <- n_batch
n_vars <- 16L

n_repeats <- 3L
warmup_runs <- 1L
seed_value <- 12345L
sample_interval_sec <- 0.05

multisession_workers <- c(2L, 4L)
output_dir <- "inst/extdata/benchmarks/results"

# ---------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------

if (!file.exists("DESCRIPTION")) {
    stop("Run this script from the TesiproV package root.", call. = FALSE)
}

for (pkg in c("devtools", "callr", "ps", "future")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(
            sprintf("Package '%s' is required. Install it before running this benchmark.", pkg),
            call. = FALSE
        )
    }
}

package_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_file <- file.path(
    output_dir,
    sprintf("mc_is_future_memory_%s.csv", benchmark_version)
)

# ---------------------------------------------------------------------
# Process-tree RSS sampler
# ---------------------------------------------------------------------

# Sum the resident set size of the benchmark R process and its descendants.
# This includes multisession workers when they are descendants of the
# benchmark process. RSS is an approximate aggregate: shared pages can be
# counted in more than one process.
sample_process_tree_rss <- function(root_pid) {
    root <- tryCatch(ps::ps_handle(root_pid), error = function(e) NULL)
    if (is.null(root)) {
        return(NA_real_)
    }

    children <- tryCatch(
        ps::ps_children(root, recursive = TRUE),
        error = function(e) list()
    )

    handles <- c(list(root), children)

    rss <- vapply(handles, function(h) {
        tryCatch(
            as.numeric(ps::ps_memory_info(h)[["rss"]]),
            error = function(e) NA_real_
        )
    }, numeric(1))

    rss <- rss[is.finite(rss)]
    if (!length(rss)) NA_real_ else sum(rss)
}

# ---------------------------------------------------------------------
# One benchmark run in a fresh R process
# ---------------------------------------------------------------------

run_one <- function(plan_name, workers, run_id) {
    ready_file <- tempfile("mcis_ready_")
    result_file <- tempfile("mcis_result_", fileext = ".rds")

    on.exit(unlink(c(ready_file, result_file)), add = TRUE)

    child_fun <- function(package_dir, plan_name, workers,
                          n_batch, n_max, n_vars, seed_value,
                          ready_file, result_file) {
        suppressPackageStartupMessages(
            devtools::load_all(package_dir, quiet = TRUE)
        )

        if (identical(plan_name, "sequential")) {
            future::plan(future::sequential)
        } else {
            future::plan(future::multisession, workers = workers)
        }

        vars <- lapply(seq_len(n_vars), function(i) {
            PROB_BASEVAR(
                Id = i,
                Name = paste0("X", i),
                DistributionType = "norm",
                Mean = 0,
                Sd = 1
            )
        })

        lDistr <- lapply(vars, function(v) v$getlDistr()[[1]])

        # Vectorized LSF: failure occurs when the first variable is >= 4.
        lsf <- function(x) {
            if (is.matrix(x)) {
                4 - x[, 1L]
            } else {
                4 - x[1L]
            }
        }

        dps <- c(4, rep(0, n_vars - 1L))

        # Signal the parent that package loading and setup are complete.
        # The parent starts memory sampling from this point onward.
        writeLines("ready", ready_file)

        timing <- system.time({
            result <- MC_IS(
                lsf = lsf,
                lDistr = lDistr,
                backend = "future",
                # The active Future plan controls the worker count.
                use_threads = NULL,
                dps = dps,
                seed = seed_value,
                n_batch = n_batch,
                n_max = n_max,
                cov_user = 0,
                dataRecord = FALSE,
                debug.level = 0
            )
        })

        saveRDS(
            list(
                elapsed_sec = unname(timing[["elapsed"]]),
                result = result,
                actual_workers = future::nbrOfWorkers()
            ),
            result_file
        )

        invisible(NULL)
    }

    job <- callr::r_bg(
        func = child_fun,
        args = list(
            package_dir = package_dir,
            plan_name = plan_name,
            workers = workers,
            n_batch = n_batch,
            n_max = n_max,
            n_vars = n_vars,
            seed_value = seed_value,
            ready_file = ready_file,
            result_file = result_file
        ),
        stdout = "|",
        stderr = "|",
        supervise = TRUE
    )

    root_pid <- job$get_pid()
    peak_rss_bytes <- 0
    baseline_rss_bytes <- NA_real_
    started_sampling <- FALSE

    while (job$is_alive()) {
        if (file.exists(ready_file)) {
            rss_now <- sample_process_tree_rss(root_pid)

            if (is.finite(rss_now)) {
                if (!started_sampling) {
                    baseline_rss_bytes <- rss_now
                    started_sampling <- TRUE
                }
                peak_rss_bytes <- max(peak_rss_bytes, rss_now)
            }
        }

        job$poll_io(timeout = 0)
        Sys.sleep(sample_interval_sec)
    }

    job$wait()

    # Take one final sample in case the last peak occurred between polls.
    rss_final <- sample_process_tree_rss(root_pid)
    if (is.finite(rss_final)) {
        if (!started_sampling) {
            baseline_rss_bytes <- rss_final
            started_sampling <- TRUE
        }
        peak_rss_bytes <- max(peak_rss_bytes, rss_final)
    }

    if (!file.exists(result_file)) {
        output <- tryCatch(job$get_output(), error = function(e) "")
        error <- tryCatch(job$get_error(), error = function(e) NULL)

        stop(
            paste(
                "Benchmark child failed.",
                if (!is.null(error)) conditionMessage(error),
                output,
                sep = "\n"
            ),
            call. = FALSE
        )
    }

    child_result <- readRDS(result_file)

    list(
        plan = plan_name,
        workers = child_result$actual_workers,
        run = run_id,
        elapsed_sec = child_result$elapsed_sec,
        baseline_tree_rss_mb = baseline_rss_bytes / 1024^2,
        peak_tree_rss_mb = peak_rss_bytes / 1024^2,
        peak_incremental_tree_rss_mb =
            max(peak_rss_bytes - baseline_rss_bytes, 0) / 1024^2,
        n_mc = child_result$result$n_mc,
        beta = child_result$result$beta,
        pf = child_result$result$pf,
        log_pf = child_result$result$log_pf
    )
}

# ---------------------------------------------------------------------
# Scenarios and repetitions
# ---------------------------------------------------------------------

scenarios <- c(
    list(list(name = "sequential", workers = 1L)),
    lapply(multisession_workers, function(w) {
        list(name = paste0("multisession_", w), workers = w)
    })
)

rows <- list()
row_id <- 0L

cat("MC_IS Future memory benchmark\n")
cat("Version:", benchmark_version, "\n")
cat("Samples per run:", n_batch, "\n")
cat("Number of variables:", n_vars, "\n")
cat(
    "Approximate u_mat size (MB):",
    n_batch * n_vars * 8 / 1024^2, "\n"
)
cat("Warm-up runs per scenario:", warmup_runs, "\n")
cat("Measured runs per scenario:", n_repeats, "\n")

for (scenario in scenarios) {
    cat("\n------------------------------------------------------------\n")
    cat("Plan:", scenario$name, "| requested workers:", scenario$workers, "\n")
    cat("------------------------------------------------------------\n")

    for (i in seq_len(warmup_runs)) {
        cat("Warm-up", i, "...\n")
        invisible(run_one(scenario$name, scenario$workers, run_id = 0L))
    }

    for (i in seq_len(n_repeats)) {
        cat(sprintf("Measured run %d/%d ...\n", i, n_repeats))

        out <- run_one(scenario$name, scenario$workers, run_id = i)
        out$version <- benchmark_version
        out$n_vars <- n_vars
        out$n_batch <- n_batch

        row_id <- row_id + 1L
        rows[[row_id]] <- as.data.frame(out, stringsAsFactors = FALSE)

        cat(sprintf(
            paste0(
                "  %.3f s | workers = %d | peak process-tree RSS = %.1f MB",
                " | incremental peak = %.1f MB | beta = %.12f\n"
            ),
            out$elapsed_sec,
            out$workers,
            out$peak_tree_rss_mb,
            out$peak_incremental_tree_rss_mb,
            out$beta
        ))
    }
}

benchmark_results <- do.call(rbind, rows)

# ---------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------

summary_results <- do.call(
    rbind,
    lapply(split(benchmark_results, benchmark_results$plan), function(d) {
        data.frame(
            version = d$version[1],
            plan = d$plan[1],
            workers = d$workers[1],
            n_vars = d$n_vars[1],
            n_batch = d$n_batch[1],
            n_repeats = nrow(d),
            median_time_sec = median(d$elapsed_sec),
            median_baseline_tree_rss_mb = median(d$baseline_tree_rss_mb),
            median_peak_tree_rss_mb = median(d$peak_tree_rss_mb),
            median_peak_incremental_tree_rss_mb =
                median(d$peak_incremental_tree_rss_mb),
            beta = d$beta[1],
            pf = d$pf[1],
            max_beta_diff = max(abs(d$beta - d$beta[1])),
            max_pf_diff = max(abs(d$pf - d$pf[1]))
        )
    })
)

cat("\n================ Raw results ================\n")
print(benchmark_results, row.names = FALSE)

cat("\n================ Summary ====================\n")
print(summary_results, row.names = FALSE)

utils::write.csv(
    benchmark_results,
    output_file,
    row.names = FALSE
)

cat("\nResults written to:\n", normalizePath(output_file), "\n")
