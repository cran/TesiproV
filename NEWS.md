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

© 2021–2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022), TU Dortmund University – Chair of Structural Concrete (since 2023)
