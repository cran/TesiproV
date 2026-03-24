
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TesiproV

<!-- badges: start -->
<!-- badges: end -->

The goal of TesiproV is to provide a flexible framework for
probabilistic reliability analysis of structural systems. It implements
several state-of-the-art algorithms (FORM, SORM, MVFOSM and Monte-Carlo
methods) and supports parametric studies through reference classes such
as `SYS_PROB`, `SYS_PARAM`, and `SYS_LSF`.

## Installation

The stable release will be available on CRAN soon. Until then you can
install directly from our GitLab server:

``` r
# install.packages("remotes")
remotes::install_gitlab("ls-massivbau/tesiprov", host = "gitlab.tu-dortmund.de")
```

You can install the development version of TesiproV via the gitLab
Server of the TU Dortmund University:
<https://gitlab.tu-dortmund.de/ls-massivbau/tesiprov>

If you need any permissions, please contact the project authors (see
below).

## Parallel computing configuration

TesiproV automatically configures parallel execution using the parallel
and future packages. By default up to 4 cores are used for outer
parallelisation (ProbMachines), while inner Monte-Carlo simulations
manage their own threads. You can change this behaviour with environment
variables before loading the package:

``` r
Sys.setenv(TesiproV.max_workers = "8") # allow up to 8 workers
Sys.setenv(TesiproV.future = "multisession") # choose backend explicitly
library(TesiproV)
```

For details see the help page of `.onLoad()` or your system’s
documentation.

## Example

This is a basic example which shows how to define a limit-state function
and run a FORM analysis:

``` r
library("TesiproV")
Var1 <- PROB_BASEVAR(Id = 1, Name = "X1", DistributionType = "norm", Mean = 0.25, Sd = 1) # kN/m²
Var2 <- PROB_BASEVAR(Id = 2, Name = "X2", DistributionType = "norm", Mean = 0.25, Sd = 1) # m

lsf <- SYS_LSF(vars = list(Var1, Var2), name = "UQLab 2d_hat")
lsf$func <- function(X1, X2) {
  return(20 - (X1 - X2)^2 - 8 * (X1 + X2 - 4)^3)
}

form <- PROB_MACHINE(name = "FORM Rack.-Fieß.", fCall = "FORM")

ps <- SYS_PROB(
  sys_input = list(lsf),
  probMachines = list(form)
)
ps$runMachines()
```

``` r
ps$beta_single
#>                     UQLab 2d_hat
#> FORM Rack.-Fieß. 3.4345652859862
```

## Parametric study example

A simple sweep over different mean values can be performed with
`PARAM_BASEVAR` and `SYS_PARAM`:

``` r
pvar <- PARAM_BASEVAR(
  Name = "E_mod",
  DistributionType = "norm",
  ParamType = "Mean",
  ParamValues = c(30e3, 35e3, 40e3)
)

lsf_param <- SYS_LSF(vars = list(pvar), name = "ParamStudy")
lsf_param$func <- function(E_mod) {
  E_mod / 1000 - 30
}

machine_form <- PROB_MACHINE(name = "FORM", fCall = "FORM")

ps_param <- SYS_PARAM(
  sys_input = list(lsf_param),
  probMachines = list(machine_form)
)
ps_param$runMachines()
print(ps_param$beta_params)
```

## Recent changes

Version 0.9.5 introduces major refactoring:

- Improved consistency between probabilistic objects (`SYS_PROB`,
  `SYS_PARAM`, …).
- Corrected implementation of the Gumbel distribution in `$prepare()`.
- Enhanced error handling in `$check()` for limit-state functions.
- Automatic parallel plan configuration via `.onLoad()` with capped
  worker pool (default 4 cores).

For a complete changelog see **NEWS.md**.

## Further information

For more, check out the vignette! There are plenty more demonstrations
and detailed explanations for each object class.

You can open it directly in R:

``` r
vignette("TesiproV")
#> Warning: vignette 'TesiproV' not found
```

# Contact & Citation

If you use TesiproV for research or teaching, please cite it as:

> Nille-Hauf K., Lux T., Feiri T., Ricker M. (2026): *TesiproV —
> Probabilistic Reliability Analysis Framework*. TU Dortmund University.

For questions or bug reports please contact
<massivbau.ab@tu-dortmund.de>.

------------------------------------------------------------------------

## License

TesiproV is licensed under the [MIT
License](https://opensource.org/licenses/MIT).  
You are free to use, modify and redistribute the software provided that
the original copyright notice and license terms are retained.

© 2021-2026 Konstantin Nille-Hauf, Til Lux, Tania Feiri, Marcus Ricker —
Hochschule Biberach / TU Dortmund University.
