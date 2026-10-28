# Benchmark MC_IS with future:
# sequential plan versus multisession with multiple workers.
#
# Run from the package root:
#   source("inst/extdata/benchmarks/benchmark_MC_IS_future_chunking.R")
#
# Run once before and once after the chunking change. Change
# benchmark_version between runs so the CSV files are not overwritten.

# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------

benchmark_version <- "after" # change to "after" for the second run

n_batch <- 100000L
n_max <- n_batch
n_repeats <- 5L
warmup_runs <- 1L
seed_value <- 12345L
output_dir <- "inst/extdata/benchmarks/results"

# Worker counts tested with future::multisession.
multisession_workers <- c(2L, 4L)

# ---------------------------------------------------------------------
# Setup and validation
# ---------------------------------------------------------------------

if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the TesiproV package root.", call. = FALSE)
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required.", call. = FALSE)
}

if (!requireNamespace("future", quietly = TRUE)) {
  stop("Package 'future' is required.", call. = FALSE)
}

devtools::load_all(".", quiet = TRUE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_file <- file.path(
  output_dir,
  sprintf("mc_is_future_chunking_%s.csv", benchmark_version)
)

# ---------------------------------------------------------------------
# Deterministic benchmark problem
# ---------------------------------------------------------------------

# A small, vectorizable LSF. Failure occurs when X >= 4 for X ~ N(0, 1).
# The supplied design point avoids running FORM inside the timed call.
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

# ---------------------------------------------------------------------
# Plan helper: set a plan for one scenario, then restore the user's plan
# ---------------------------------------------------------------------

with_benchmark_plan <- function(strategy, code, ...) {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(strategy, ...)
  force(code)
}

# ---------------------------------------------------------------------
# Run one fixed-size MC_IS calculation
# ---------------------------------------------------------------------

run_one <- function() {
  timing <- system.time({
    result <- MC_IS(
      lsf = lsf,
      lDistr = list(marginal),
      backend = "future",
      # Do not supply use_threads: the active Future plan sets the
      # number of Future workers.
      dps = 4,
      seed = seed_value,
      n_batch = n_batch,
      n_max = n_max,
      cov_user = 0,
      dataRecord = FALSE,
      debug.level = 0
    )
  })

  list(
    elapsed_sec = unname(timing[["elapsed"]]),
    result = result
  )
}

# ---------------------------------------------------------------------
# Scenario definitions
# ---------------------------------------------------------------------

scenarios <- c(
  list(list(name = "sequential", strategy = future::sequential, workers = NA_integer_)),
  lapply(multisession_workers, function(w) {
    list(name = paste0("multisession_", w), strategy = future::multisession, workers = w)
  })
)

# ---------------------------------------------------------------------
# Run benchmark
# ---------------------------------------------------------------------

results <- list()
row_id <- 0L

for (scenario in scenarios) {
  cat("\n------------------------------------------------------------\n")
  cat(sprintf(
    "Version: %s | Plan: %s | Workers: %s\n",
    benchmark_version,
    scenario$name,
    if (is.na(scenario$workers)) "1 (sequential)" else scenario$workers
  ))
  cat("------------------------------------------------------------\n")

  scenario_results <- with_benchmark_plan(
    scenario$strategy,
    {
      if (!is.na(scenario$workers)) {
        future::plan(scenario$strategy, workers = scenario$workers)
      }

      # Confirm the worker count actually selected by the plan.
      workers_actual <- future::nbrOfWorkers()
      cat("Active Future workers:", workers_actual, "\n")

      # Warm-up is excluded from the measured timings.
      for (i in seq_len(warmup_runs)) {
        invisible(run_one())
      }

      measured <- vector("list", n_repeats)
      for (i in seq_len(n_repeats)) {
        measured[[i]] <- run_one()
        r <- measured[[i]]$result

        cat(sprintf(
          "Run %d/%d: %.3f s | n_mc = %d | beta = %.12f | pf = %.12e\n",
          i,
          n_repeats,
          measured[[i]]$elapsed_sec,
          r$n_mc,
          r$beta,
          r$pf
        ))
      }

      measured
    }
  )

  # First measured result is the reference for within-scenario differences.
  ref <- scenario_results[[1]]$result

  for (i in seq_along(scenario_results)) {
    out <- scenario_results[[i]]
    r <- out$result
    row_id <- row_id + 1L

    results[[row_id]] <- data.frame(
      version = benchmark_version,
      plan = scenario$name,
      workers = if (is.na(scenario$workers)) 1L else scenario$workers,
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

benchmark_results <- do.call(rbind, results)

# ---------------------------------------------------------------------
# Summary and output
# ---------------------------------------------------------------------

summary_results <- do.call(
  rbind,
  lapply(split(benchmark_results, benchmark_results$plan), function(d) {
    data.frame(
      version = d$version[1],
      plan = d$plan[1],
      workers = d$workers[1],
      n_repeats = nrow(d),
      median_time_sec = median(d$elapsed_sec),
      mean_time_sec = mean(d$elapsed_sec),
      sd_time_sec = stats::sd(d$elapsed_sec),
      median_samples_per_sec = median(d$samples_per_sec),
      beta = d$beta[1],
      pf = d$pf[1],
      max_beta_diff = max(d$beta_diff_from_first_run),
      max_pf_diff = max(d$pf_diff_from_first_run)
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
