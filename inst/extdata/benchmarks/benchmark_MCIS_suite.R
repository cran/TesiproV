# ============================================================
# TesiproV Benchmark Suite
# MC_IS_system vs MC_IS_gsys
# ============================================================

library(TesiproV)

set.seed(123)

cat("\n===========================================\n")
cat("TesiproV Extended Benchmark Suite\n")
cat("===========================================\n\n")

# ------------------------------------------------------------
# Global benchmark parameters
# ------------------------------------------------------------

n_batch <- 4000
n_max <- 8e4
threads <- 1

kappa_values <- c(10, 50, 200)

cat("Simulation parameters\n")
cat("n_batch:", n_batch, "\n")
cat("n_max:", n_max, "\n\n")

# ------------------------------------------------------------
# Helper: create variables
# ------------------------------------------------------------

create_vars <- function(dim) {
    lapply(seq_len(dim), function(i) {
        PROB_BASEVAR(
            Id = i,
            Name = paste0("X", i),
            DistributionType = "norm",
            Mean = 0,
            Sd = 1
        )
    })
}

# ------------------------------------------------------------
# Helper: build SYS_LSF objects
# ------------------------------------------------------------

build_lsf_objects <- function(lsf_list, vars) {
    lapply(seq_along(lsf_list), function(i) {
        obj <- SYS_LSF(
            vars = vars,
            name = paste0("LSF", i)
        )

        obj$func <- lsf_list[[i]]
        obj$check()

        obj
    })
}

# ------------------------------------------------------------
# Scenario definitions
# ------------------------------------------------------------

scenario_list <- list(
    list(
        name = "independent",
        dim = 2,
        lsf = list(
            function(x) 3 - x[1],
            function(x) 3 - x[2]
        )
    ),
    list(
        name = "correlated",
        dim = 2,
        lsf = list(
            function(x) 3 - x[1],
            function(x) 3 - (0.7 * x[1] + 0.7 * x[2])
        )
    ),
    list(
        name = "orthogonal",
        dim = 2,
        lsf = list(
            function(x) x[1] - 2,
            function(x) x[2] - 2
        )
    )
)

# ------------------------------------------------------------
# Benchmark container
# ------------------------------------------------------------

results <- list()

row_id <- 1

# ============================================================
# Benchmark loop
# ============================================================

for (scenario in scenario_list) {
    cat("\n-------------------------------------------\n")
    cat("Scenario:", scenario$name, "\n")
    cat("-------------------------------------------\n")

    vars <- create_vars(scenario$dim)

    sys_input <- build_lsf_objects(
        scenario$lsf,
        vars
    )

    # Prepare distribution list for MC_IS_system

    # Normalize distribution objects exactly as done in tests

    normalize_distr <- function(v) {
        d <- v$getlDistr()
        if (is.list(d) && length(d) >= 1 && is.list(d[[1]]) && !is.null(d[[1]]$p)) {
            d[[1]]
        } else {
            d
        }
    }

    d1 <- normalize_distr(vars[[1]])
    d2 <- normalize_distr(vars[[2]])

    lDistr <- list(
        list(d1),
        list(d1, d2)
    )

    # ------------------------------------------------------------
    # Run classical MC_IS_system
    # ------------------------------------------------------------

    cat("Running MC_IS_system...\n")

    time_system <- system.time({
        res_system <- MC_IS(
            lsf = as.list(scenario$lsf),
            lDistr = lDistr,
            cov_user = 0.05,
            n_batch = n_batch,
            n_max = n_max,
            use_threads = threads
        )
    })

    cat("Finished MC_IS_system\n")

    # ------------------------------------------------------------
    # Run smooth system for multiple kappa
    # ------------------------------------------------------------

    for (kappa in kappa_values) {
        cat("Running MC_IS_gsys (kappa =", kappa, ")...\n")

        time_gsys <- system.time({
            res_gsys <- MC_IS_gsys(
                sys_input = sys_input,
                sys_structure = "serial",
                smooth_kappa = kappa,
                cov_user = 0.05,
                n_batch = n_batch,
                n_max = n_max,
                use_threads = threads,
                seed = 123
            )
        })

        results[[row_id]] <- data.frame(
            scenario = scenario$name,
            kappa = kappa,
            pf_system = res_system$pf,
            pf_gsys = res_gsys$pf,
            beta_system = res_system$beta,
            beta_gsys = res_gsys$beta,
            samples_system = res_system$n_mc,
            samples_gsys = res_gsys$n_mc,
            runtime_system = time_system[3],
            runtime_gsys = time_gsys[3],
            speedup = time_system[3] / time_gsys[3],
            pf_diff = abs(res_system$pf - res_gsys$pf)
        )

        row_id <- row_id + 1
    }
}

# ============================================================
# Combine results
# ============================================================

benchmark_table <- do.call(rbind, results)

cat("\n===========================================\n")
cat("Benchmark Results\n")
cat("===========================================\n\n")

print(benchmark_table)

# ============================================================
# Optional plots
# ============================================================

if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)

    cat("\nGenerating plots...\n")

    p1 <- ggplot(
        benchmark_table,
        aes(
            x = factor(kappa),
            y = speedup,
            fill = scenario
        )
    ) +
        geom_bar(
            stat = "identity",
            position = "dodge"
        ) +
        labs(
            title = "Speedup: MC_IS_system / MC_IS_gsys",
            x = "kappa",
            y = "Speedup"
        ) +
        theme_bw()

    p2 <- ggplot(
        benchmark_table,
        aes(
            x = factor(kappa),
            y = pf_diff,
            fill = scenario
        )
    ) +
        geom_bar(
            stat = "identity",
            position = "dodge"
        ) +
        labs(
            title = "Pf difference",
            x = "kappa",
            y = "|Pf_system - Pf_gsys|"
        ) +
        theme_bw()

    print(p1)
    print(p2)
}

cat("\n===========================================\n")
cat("Benchmark Suite completed\n")
cat("===========================================\n")
