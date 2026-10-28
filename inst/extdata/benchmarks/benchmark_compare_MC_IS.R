# Benchmark: MC_IS before the merge with master (own branch) vs. after the merge
# (with the log-space implementation of MC_IS_single).
#
# Problem: Ex. 5.11 Nowak & Collins (2000) as in tests/testthat/long_test_empdist.R
# (student-t, empirical and Gumbel variable), MC_IS with the same options as there.
#
# Start from the package root (R console):
#   source("inst/extdata/benchmarks/benchmark_compare_MC_IS.R")
#
# The script
#   1. finds the last merge commit and its first parent (= state before the merge),
#   2. checks that state (and HEAD) out as git worktrees in the temp folder of R,
#   3. installs both versions into separate temporary libraries,
#   4. runs MC_IS for several seeds in a fresh Rscript session per run,
#   5. prints the results, writes a csv file and removes the worktree.
#
# Runtime: about 2 min per run, i.e. about 12 min for 2 versions x 3 seeds.

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
seeds    <- c(1234, 1, 2)       # first value as in long_test_empdist.R
cov_user <- 0.05
n_max    <- 1e7
n_batch  <- 10000
beta_ref <- 3.9585857628075     # reference value of the long test (beta_mcis)
out_csv  <- "benchmark_compare_MC_IS.csv"

# Leave old_commit = NULL to use the first parent of the last merge commit.
# Otherwise set a commit hash, e.g. old_commit <- "abc1234".
old_commit <- NULL

# ---------------------------------------------------------------------------
# Versions
# ---------------------------------------------------------------------------
dir_new <- normalizePath(".")
if (!file.exists(file.path(dir_new, "DESCRIPTION"))) {
  stop("Please start the script from the package root.")
}
git <- function(...) system2("git", c(...), stdout = TRUE, stderr = FALSE)

if (is.null(old_commit)) {
  merge_commit <- git("log", "--merges", "-1", "--format=%H")
  if (length(merge_commit) != 1L) stop("No merge commit found, set old_commit manually.")
  old_commit <- git("rev-parse", paste0(merge_commit, "^1"))
  cat("Merge commit:", merge_commit, "\n")
}
cat("State before the merge:", old_commit, "\n")
cat(paste(git("log", "--oneline", "-1", old_commit), collapse = "\n"), "\n\n")

# Both versions are checked out as git worktrees into the temporary folder of R.
# That path contains no blanks or umlauts, which R CMD INSTALL does not handle
# reliably on Windows.
dir_old <- file.path(tempdir(), "tesipro_old")
dir_new_wt <- file.path(tempdir(), "tesipro_new")
for (d in c(dir_old, dir_new_wt)) {
  system2("git", c("worktree", "remove", "--force", shQuote(d)), stdout = FALSE, stderr = FALSE)
}
if (system2("git", c("worktree", "add", "--detach", shQuote(dir_old), old_commit),
            stdout = FALSE, stderr = FALSE) != 0L) stop("git worktree add (old) failed")
if (system2("git", c("worktree", "add", "--detach", shQuote(dir_new_wt), "HEAD"),
            stdout = FALSE, stderr = FALSE) != 0L) stop("git worktree add (new) failed")
if (!file.exists(file.path(dir_old, "DESCRIPTION"))) {
  system2("git", c("worktree", "remove", "--force", shQuote(dir_old)), stdout = FALSE, stderr = FALSE)
  stop("No DESCRIPTION in ", dir_old)
}
# Uncommitted changes of tracked files are copied into the worktree of the new version,
# so that the current working state is measured and not only the last commit.
changed <- git("diff", "--name-only", "HEAD")
changed <- changed[file.exists(file.path(dir_new, changed))]
if (length(changed)) {
  cat("Uncommitted changes used for the 'new' version:\n")
  cat(paste0("  ", changed, collapse = "\n"), "\n\n")
  for (fn in changed) {
    file.copy(file.path(dir_new, fn), file.path(dir_new_wt, fn), overwrite = TRUE)
  }
}
# on.exit() does not work reliably in a sourced script, so the cleanup is called explicitly.
cleanup <- function() {
  for (d in c(dir_old, dir_new_wt)) {
    system2("git", c("worktree", "remove", "--force", shQuote(d)), stdout = FALSE, stderr = FALSE)
  }
  system2("git", c("worktree", "prune"), stdout = FALSE, stderr = FALSE)
}

# Both versions are installed into separate libraries. Installed packages are
# used on purpose, because future workers load the installed package and not
# the load_all state of the main session.
libs <- list()
rbin <- file.path(R.home("bin"), "R")
for (variant in c("old", "new")) {
  libs[[variant]] <- file.path(tempdir(), paste0("lib_", variant))
  dir.create(libs[[variant]], showWarnings = FALSE)
  src <- if (variant == "old") dir_old else dir_new_wt
  cat(sprintf("Installing %s version from %s ...\n", variant, src))
  out <- suppressWarnings(
    system2(rbin, c("CMD", "INSTALL", "--no-docs", "--no-multiarch",
                    paste0("--library=", shQuote(libs[[variant]])), shQuote(src)),
            stdout = TRUE, stderr = TRUE)
  )
  st <- attr(out, "status")
  if (!is.null(st) && st != 0L) {
    cat("\nInstallation output (last lines):\n")
    print(utils::tail(out, 25))
    cleanup()
    stop("Installation of the ", variant, " version failed")
  }
}

# ---------------------------------------------------------------------------
# Worker script (one fresh R session per run)
# ---------------------------------------------------------------------------
worker <- tempfile(fileext = ".R")
writeLines(r"---(
args <- commandArgs(trailingOnly = TRUE)
lib_dir <- args[1]
seed <- as.numeric(args[2])
cov_user <- as.numeric(args[3])
n_max <- as.numeric(args[4])
n_batch <- as.numeric(args[5])
rds_out <- args[6]

.libPaths(c(lib_dir, .libPaths()))
suppressMessages(library(TesiproV))

# Fy: empirical distribution from a log-normal sample (Mean = 40, Cov = 0.1)
set.seed(1)
sig_fy <- sqrt(log(1 + 0.1^2))
mue_fy <- log(40) - sig_fy^2 / 2
data_fy <- rlnorm(100000, mue_fy, sig_fy)

X1 <- PROB_BASEVAR(
  Id = 1, Name = "Z", DistributionType = "st",
  DistributionParameters = c(100, 4, 10000, 10000)
)
X2 <- PROB_BASEVAR(
  Id = 2, Name = "Fy", DistributionType = "emp",
  Package = "EnvStats", DistributionParameters = data_fy
)
X3 <- PROB_BASEVAR(
  Id = 3, Name = "M", DistributionType = "gumbel",
  Package = "evd", Mean = 2000, Cov = 0.1
)

lsf1 <- SYS_LSF(vars = list(X1, X2, X3), name = "Nowak & Collins Exp5.11_1")
lsf1$func <- function(M, Fy, Z) Z * Fy - M
lsf1$check()

mcis <- PROB_MACHINE(
  name = "MC IS", fCall = "MC_IS",
  options = list(
    cov_user = cov_user, n_max = n_max, n_batch = n_batch,
    use_threads = 1L, seed = seed, backend = "future"
  )
)
ps <- SYS_PROB(sys_input = list(lsf1), probMachines = list(mcis), debug.level = 0)

t0 <- proc.time()[["elapsed"]]
invisible(suppressMessages(ps$runMachines()))
el <- proc.time()[["elapsed"]] - t0

r <- ps$res_single[[1]][[1]]
saveRDS(r, rds_out)
pick <- function(nms) {
  for (nm in nms) if (!is.null(r[[nm]])) return(as.numeric(r[[nm]][1]))
  NA_real_
}
cat(sprintf("RESULT;%.6f;%.6e;%.5f;%.0f;%.2f\n",
            ps$beta_single[1], pick("pf"), pick(c("cov", "CoV")),
            pick(c("n_sim", "n_total", "N", "n_samples", "nsim")), el))
cat("FIELDS;", paste(names(r), collapse = ","), "\n", sep = "")
)---", worker)

# ---------------------------------------------------------------------------
# Runs
# ---------------------------------------------------------------------------
rows <- list()
raw <- list()
fields <- NULL
for (s in seeds) {
  for (variant in c("old", "new")) {
    cat(sprintf("[%s] seed %s ...\n", variant, format(s)))
    rds_file <- file.path(tempdir(), sprintf("r_%s_%s.rds", variant, format(s)))
    out <- system2("Rscript",
                   c(shQuote(worker), shQuote(libs[[variant]]), s, cov_user, n_max, n_batch,
                     shQuote(rds_file)),
                   stdout = TRUE, stderr = FALSE)
    line <- grep("^RESULT;", out, value = TRUE)
    if (length(line) != 1L) {
      cat("  run failed, last output:\n")
      print(utils::tail(out, 8))
      next
    }
    if (is.null(fields)) {
      fl <- grep("^FIELDS;", out, value = TRUE)
      if (length(fl)) fields <- sub("^FIELDS;", "", fl[1])
    }
    if (file.exists(rds_file)) raw[[paste(variant, format(s))]] <- readRDS(rds_file)
    p <- suppressWarnings(as.numeric(strsplit(sub("^RESULT;", "", line), ";")[[1]]))
    rows[[length(rows) + 1L]] <- data.frame(
      variant = variant, seed = s, beta = p[1], pf = p[2],
      cov = p[3], n_sim = p[4], seconds = p[5]
    )
  }
}

cleanup()
if (!length(rows)) stop("No run succeeded.")
res <- do.call(rbind, rows)
utils::write.csv(res, out_csv, row.names = FALSE)

cat("\nSingle runs:\n")
print(res, row.names = FALSE)

cat("\nMeans per version:\n")
num <- c("beta", "pf", "cov", "n_sim", "seconds")
agg <- do.call(rbind, lapply(split(res, res$variant), function(d) {
  data.frame(variant = d$variant[1], t(colMeans(d[num], na.rm = TRUE)))
}))
print(agg, row.names = FALSE)

cat("\nRelative spread over the seeds (sd / mean), the target of cov_user is", cov_user, ":\n")
for (v in unique(res$variant)) {
  d <- res[res$variant == v, ]
  cat(sprintf("  %-4s pf: %.4f   beta: %.5f   seconds: %.2f\n", v,
              sd(d$pf) / mean(d$pf), sd(d$beta) / mean(d$beta), sd(d$seconds) / mean(d$seconds)))
}

if (all(c("old", "new") %in% res$variant)) {
  t_old <- mean(res$seconds[res$variant == "old"])
  t_new <- mean(res$seconds[res$variant == "new"])
  cat(sprintf("\nRuntime new/old: %.2f\n", t_new / t_old))
  for (v in c("old", "new")) {
    b <- res$beta[res$variant == v]
    cat(sprintf("beta %s: mean %.4f, deviation from reference %.2f %%\n",
                v, mean(b), 100 * (mean(b) - beta_ref) / beta_ref))
  }
}
if (length(raw)) {
  saveRDS(raw, "benchmark_compare_MC_IS_raw.rds")
  cat("\nStructure of the first result object per version:\n")
  for (v in c("old", "new")) {
    k <- grep(paste0("^", v, " "), names(raw))[1]
    if (!is.na(k)) {
      cat("\n--", v, "--\n")
      utils::str(raw[[k]], max.level = 1, give.attr = FALSE, vec.len = 4)
    }
  }
  cat("\nAll result objects saved in benchmark_compare_MC_IS_raw.rds\n")
}
cat("\nResults written to", out_csv, "\n")
