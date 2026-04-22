#library(devtools)
#devtools::load_all("~/projects/tesiprov")

options(digits = 14)

############################################################
## Common MC_IS settings
############################################################

cov_user <- 0.05
n_batch <- 250000
n_max <- 2e6
seed_val <- 12345

thread_list <- c(1, 4, 8) # , 8)
n_repeat <- 2

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
    x[1] * x[2] - x[3]
}

lsf_medium <- function(x) {
    z <- x[1]
    fy <- x[2]
    m <- x[3]
    z * fy - m + 0.01 * z^2
}

lsf_heavy <- function(x) {
    z <- x[1]
    fy <- x[2]
    m <- x[3]

    val <- 0
    for (k in 1:200) {
        val <- val + sin(z * 0.001 * k) * cos(fy * 0.002 * k)
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

    times <- numeric(length(thread_list))

    for (t in seq_along(thread_list)) {
        cores <- thread_list[t]
        cat("Threads:", cores, "\n")

        runtime_vec <- numeric(n_repeat)

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
                    "densityType" = "norm"
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
        }

        times[t] <- mean(runtime_vec)
    }

    speedup <- times[1] / times

    result <- data.frame(
        threads = thread_list,
        time_sec = times,
        speedup = speedup,
        efficiency = speedup / thread_list
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
