# TesiproV 1.0.0

Consolidated release of the package, based on more than two years of use.
New in this version: the influence of quality controls on the distribution of 
basic variables can be taken into account.

## New: quality control
- The influence of quality controls on the distribution of a basic variable
  can now be taken into account. A sampling plan is defined with the new
  class `PROB_QUALCON`. The conformity criterion of EN 206 for the compressive
  strength (`"fck_EN-206"`) and an attribute plan for dimensions and other
  characteristics (`"dim_k||n"`) are available. The acceptance probability of
  a plan (OC curve) can be computed, plotted or entered manually.
- `PROB_BASEVAR$updateBasevar()` returns the basic variable after acceptance of
  the lot (posterior distribution). The prior must be a Student-t (`"st"`) or
  log-Student-t (`"lt"`) distribution with hyperparameters. The message compares
  mean value and standard deviation of prior and posterior.
- The vignette contains a new section with examples.

## System reliability
- The importance sampling for systems (`"MCIS"`) was revised. It estimates the
  system failure probability more accurately and also handles parallel systems
  whose design points are far apart.
- Available methods of `calculateSystemProbability()`: `"simpleBounds"`,
  `"MCIS"`, `"MCC"` (crude Monte Carlo) and `"MCSUS"` (subset simulation).
  FORM and SORM provide single results, which enter the simple bounds. A system
  analysis by FORM or SORM is not available.
- Variables that are shared by several limit-state functions must have the
  same `Name`.
- The effective sample size reported by `MC_IS` now counts the samples in the
  failure domain.

## Other changes
- FORM with `optim_type = "auto"` (default) switches to a second solver if the
  first one does not converge. The result field `fallback` shows whether this
  happened.
- The console output of the Monte Carlo methods is reduced and controlled by
  `debug.level`.
- `PROB_BASEVAR$plot()` plots the probability density of a basic variable. For a variable created by `updateBasevar()`, the prior distribution can be added with `prior = TRUE` and the fitted curve and the empirical resample with `compare = TRUE`.
- Option descriptions and notes on parallel computing in the vignette were
  updated.


---

# TesiproV 0.9.6

## MC_IS estimator correction and reliability analysis improvements

### Monte Carlo Importance Sampling (MC_IS)
- Corrected the MC_IS estimator to remove bias caused by the previous self-normalized implementation.
- Refactored worker aggregation in `MC_IS_single` to ensure correct accumulation of 
  weighted samples and consistent global statistics.
- Workers now compute weights directly and return aggregated statistics instead of raw sample vectors.

### System reliability analysis
- Fixed incorrect propagation of `sys_type` from `calculateSystemProbability()` to `MC_IS_system()`, 
  ensuring that serial and parallel systems are evaluated correctly.
- Resolved cases of uninitialized variables in the system worker (`mode_ids`, `shift`, and references to `I_matrix`).
- Enabled multimodal sampling for both serial and parallel systems.

### Performance improvements
- Optimized the transformation from u-space to x-space by removing large intermediate 
  matrices and computing quantiles column-wise.
- Replaced large failure indicator matrices with incremental counters to significantly reduce memory usage.
- Reduced worker-to-master communication in parallel runs by returning only aggregated statistics.

### Validation and testing
- Added regression tests for `MC_IS_single`, including linear limit-state cases, 
  non-normal distributions, and rare-event scenarios.
- Added system reliability tests comparing results with `MC_CRUDE`.
- Verified consistency with FORM results within expected Monte Carlo error.


---

# TesiproV 0.9.5

## Major refactoring of object classes

### System and reliability objects
- Reworked all reference classes (`SYS_PROB`, `SYS_PARAM`, `SYS_LSF`, `PROB_BASEVAR`, …).
- Introduced explicit constructors with argument validation for better error handling.
- Improved consistency between mean, standard deviation and coefficient of variation.
- Added robust forward/backward transformations for all supported distributions.
- Fixed handling of zero mean values (`Cov <- Inf`) to reflect infinite relative scatter.

### Limit-state functions
- Field type for `func` changed from `"function"` to `"ANY"` with explicit default `NULL`.
  This prevents creation of dummy closures and allows proper missing-function detection.
- Enhanced `$check()` method:
  - Distinguishes between *missing*, *empty*, and *valid* limit-state functions.
  - Provides informative error messages used by automated tests.

### Parametric study objects
- Constructor of `SYS_PARAM` now correctly forwards arguments to parent class (`callSuper(...)`)
  ensuring that inherited fields like `sys_input` are preserved.
- Parallel execution restructured using the global future plan set in `.onLoad()`.
  The outer level is capped by environment variable `TesiproV.max_workers`
  (default: **4 cores**) while inner Monte-Carlo methods manage their own threads.

### Probabilistic base variables
- Complete rewrite of transformation logic in `$prepare()`:
  - Supports normal, lognormal, Gumbel, gamma, beta and Weibull distributions consistently.
  - **Corrected implementation of the Gumbel distribution**  
    (location parameter now properly computed as  
    `location = Mean + digamma(1) * scale`; scale derived via `(Sd * sqrt(6)) / π`).
  - Added numerical root search for Weibull shape parameter when empirical formula fails.
  - Added empirical, shifted lognormal, Student-t and logStudent-t distribution.
- Distribution caching via digest hash implemented in `$getlDistr()` for performance improvement.

### Package infrastructure
- Added internal `.onLoad()` function for package-specific default options that can later be used by TesiproV
functions to configure parallel execution in a controlled and CRAN-compliant way.

## Minor improvements
- Better formatted output in `$printResults()`.
- Consistent English comments throughout codebase for clarity and documentation generation via roxygen2.
- Updated examples in help pages to match new constructors and usage patterns.

## Compatibility notes
Existing scripts using older class definitions should continue to work after minor adjustments:
replace direct field assignments by constructor arguments where applicable,
and ensure that limit-state functions are defined before calling `$check()`.


---

# TesiproV 0.9.2.0

Better implementation of parallel computing in Unix and macOS environments.  
Changed RNG handling in parallel approach.


---

# TesiproV 0.9.1.0

* First release.


---

© 2021–2026 T. Lux, K. Nille-Hauf, T. Feiri, M. Ricker -- Hochschule Biberach (until 2022), TU Dortmund University – Chair of Structural Concrete (since 2023)
