# library(devtools)
# devtools::load_all("~/projects/tesiprov")

options(digits = 14)

############################################################
## Common MC_IS settings
############################################################

cov_user <- 0.0
n_batch <- 250000L
n_max <- 2e6
seed_val <- 12345

thread_list <- c(1, 4, 8) # , 8)
n_repeat <- 3L

############################################################
## Base variables
############################################################

X1 <- PROB_BASEVAR(Id = 1, Name = "Z", DistributionType = "norm", Mean = 100, Cov = 0.04)
X2 <- PROB_BASEVAR(Id = 2, Name = "Fy", DistributionType = "lnorm", Mean = 40, Cov = 0.1)
X3 <- PROB_BASEVAR(Id = 3, Name = "M", DistributionType = "gumbel", Package = "evd", Mean = 2000, Cov = 0.1)

############################################################
## LSF Definitions
############################################################

lsf_simple <- function(x) {
    if (is.matrix(x)) {
        z <- x[, 1L]
        fy <- x[, 2L]
        m <- x[, 3L]
    } else {
        z <- x[1L]
        fy <- x[2L]
        m <- x[3L]
    }

    z * fy - m
}

lsf_medium <- function(x) {
    if (is.matrix(x)) {
        z <- x[, 1L]
        fy <- x[, 2L]
        m <- x[, 3L]
    } else {
        z <- x[1L]
        fy <- x[2L]
        m <- x[3L]
    }

    z * fy - m + 0.01 * z^2
}

lsf_heavy <- function(x) {
    if (is.matrix(x)) {
        z <- x[, 1L]
        fy <- x[, 2L]
        m <- x[, 3L]
    } else {
        z <- x[1L]
        fy <- x[2L]
        m <- x[3L]
    }

    val <- numeric(length(z))

    for (k in seq_len(200L)) {
        val <- val +
            sin(z * 0.001 * k) *
                cos(fy * 0.002 * k)
    }

    z * fy - m + val * 1e-6
}

############################################################
## Benchmark runner
############################################################

run_benchmark <- function(lsf_fun, label) {
    cat("\n==============================\n")
    cat("Testing:", label, "\n")
    cat("==============================\n")

    lsf <- SYS_LSF(vars = list(X1, X2, X3), name = label)
    lsf$func <- lsf_fun
    lsf$check()

    lDistr_form <- lapply(lsf$vars, function(v) v$getlDistr())
    form_result <- FORM(lsf$getLSF(), lDistr_form)
    dps <- form_result$x_points

    mean_time <- numeric(length(thread_list))
    median_time <- numeric(length(thread_list))
    sd_time <- numeric(length(thread_list))
    min_time <- numeric(length(thread_list))
    max_time <- numeric(length(thread_list))
    mean_n_mc <- numeric(length(thread_list))
    mean_beta <- numeric(length(thread_list))

    for (t in seq_along(thread_list)) {
        cores <- thread_list[t]
        cat("Threads:", cores, "\n")

        runtime_vec <- numeric(n_repeat)
        n_mc_vec <- numeric(n_repeat)
        beta_vec <- numeric(n_repeat)

        for (r in 1:n_repeat) {
            machine_mcis <- PROB_MACHINE(
                name = "MC IS",
                fCall = "MC_IS",
                options = list(
                    "n_max" = n_max,
                    "cov_user" = cov_user,
                    "use_threads" = cores,
                    "n_batch" = n_batch,
                    "backend" = "parallel",
                    "seed" = seed_val,
                    "dataRecord" = FALSE,
                    "densityType" = "norm",
                    dps = dps
                )
            )

            ps <- SYS_PROB(
                sys_input    = list(lsf),
                probMachines = list(machine_mcis),
                debug.level  = 0
            )

            tic <- proc.time()
            ps$runMachines()
            runtime_vec[r] <- (proc.time() - tic)[3]

            res <- ps$res_single[[1]][[1]]
            n_mc_vec[r] <- res$n_mc
            beta_vec[r] <- res$beta
        }

        # times[t] <- mean(runtime_vec)
        mean_time[t] <- mean(runtime_vec)
        median_time[t] <- median(runtime_vec)
        sd_time[t] <- sd(runtime_vec)
        min_time[t] <- min(runtime_vec)
        max_time[t] <- max(runtime_vec)
        mean_n_mc[t] <- mean(n_mc_vec)
        mean_beta[t] <- mean(beta_vec)
    }

    mean_speedup <- mean_time[1] / mean_time
    median_speedup <- median_time[1] / median_time

    relative_sd <- ifelse(mean_time > 0, sd_time / mean_time, NA_real_)

    result <- data.frame(
        threads = thread_list,
        mean_time_sec = mean_time,
        median_time_sec = median_time,
        sd_time_sec = sd_time,
        relative_sd = relative_sd,
        min_time_sec = min_time,
        max_time_sec = max_time,
        mean_speedup = mean_speedup,
        median_speedup = median_speedup,
        efficiency = mean_speedup / thread_list,
        samples = mean_n_mc,
        beta = mean_beta
    )

    print(result)
    return(result)
}

############################################################
## Execute benchmarks
############################################################

res_simple <- run_benchmark(lsf_simple, "Simple LSF")
res_medium <- run_benchmark(lsf_medium, "Medium LSF")
res_heavy <- run_benchmark(lsf_heavy, "Heavy LSF")
