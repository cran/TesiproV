# ======================================================================
# diagnose_MC_IS_log_weights.R
#
# Diagnostic script for numerical underflow in MC_IS_single().
#
# It demonstrates that direct exponentiation of very small log-weights
# loses all information, whereas log-sum-exp remains finite.
# ======================================================================

options(digits = 16)

# ----------------------------------------------------------------------
# Helper: stable log(sum(exp(x)))
# ----------------------------------------------------------------------

diagnostic_log_sum_exp <- function(x) {
    if (!is.numeric(x)) {
        stop("'x' must be numeric.", call. = FALSE)
    }

    if (length(x) == 0L) {
        return(-Inf)
    }

    m <- max(x)

    if (!is.finite(m)) {
        return(m)
    }

    m + log(sum(exp(x - m)))
}

# ----------------------------------------------------------------------
# Part 1: Pure numerical demonstration
# ----------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("Part 1: Direct exponentiation versus log-sum-exp\n")
cat("============================================================\n")

log_w <- c(-900, -901, -905)

naive_sum_w <- sum(exp(log_w))
stable_log_sum_w <- diagnostic_log_sum_exp(log_w)

cat("log_w:                 ", paste(log_w, collapse = ", "), "\n")
cat("sum(exp(log_w)):       ", naive_sum_w, "\n")
cat("log_sum_exp(log_w):    ", stable_log_sum_w, "\n")

stopifnot(naive_sum_w == 0)
stopifnot(is.finite(stable_log_sum_w))

cat("\nResult: Direct exponentiation underflows to zero,\n")
cat("while log-sum-exp retains the finite result.\n")

# ----------------------------------------------------------------------
# Part 2: MC_IS rare-event diagnostic
# ----------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("Part 2: MC_IS diagnostic for an ultra-rare event\n")
cat("============================================================\n")

# X ~ N(0, 1)
x_var <- PROB_BASEVAR(
    Id = 1,
    Name = "X",
    DistributionType = "norm",
    Mean = 0,
    Sd = 1
)

# Failure event:
#
#   g(x) = 40 - x <= 0
#
# thus:
#
#   X >= 40
#
# The exact reliability index is approximately beta = 40.
lsf_rare <- function(x) {
    40 - x[1]
}

# MC_IS receives a funlist, not the wrapper list(funlist, parameters).
l_distr <- list(x_var$getlDistr()[[1]])

res <- MC_IS(
    lsf = lsf_rare,
    lDistr = l_distr,
    dps = 40,
    cov_user = 1.0,
    n_batch = 2000L,
    n_max = 2000L,
    use_threads = 1L,
    backend = "parallel",
    seed = 12345L,
    dataRecord = FALSE,
    stability_mode = "robust",
    debug.level = 2
)

cat("\nMC_IS result:\n")
print(res[c("pf", "log_pf", "beta", "n_mc", "var", "cov_mc")])

cat("\nExpected behaviour after the log-weight fix:\n")
cat("  - n_mc should equal 2000\n")
cat("  - beta should be finite and approximately 40\n")
cat("  - pf may be too small to represent as an ordinary R double\n")
cat("    and should eventually be accompanied by a log-scale value.\n")

if (
    !is.null(res$log_pf) &&
        is.finite(res$log_pf) &&
        !is.null(res$beta) &&
        is.finite(res$beta)
) {
    cat("\nLOG-DOMAIN CALCULATION SUCCESSFUL:\n")
    cat("  log_pf remains finite although pf underflows to zero.\n")
    cat(sprintf("  Calculated beta: %.8f\n", res$beta))
} else {
    cat("\nDIAGNOSIS CONFIRMED:\n")
    cat("  The log-domain calculation did not retain a finite\n")
    cat("  failure probability or reliability index.\n")
}
