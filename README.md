
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TesiproV

<!-- badges: start -->
<!-- badges: end -->

The goal of TesiproV is to …

## Installation

Install the package through CRAN.

``` r
# install.packages("TesiproV")
```

You can install the development version of TesiproV via the gitLab
Server of the RWTH Aachen:
<https://git.rwth-aachen.de/jan.philip.schulze-ardey/tesiprov>

If you need any permissions, just let us know (email us).

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library("TesiproV")
Var1 <- PROB_BASEVAR(Id=1,Name="X1",DistributionType="norm",Mean=0.25,Sd=1) #kN/m²
Var2 <- PROB_BASEVAR(Id=2,Name="X2",DistributionType="norm",Mean=0.25,Sd=1) #m

lsf<-SYS_LSF(vars=list(Var1,Var2),name="UQLab 2d_hat")
lsf$func <- function(X1,X2){
  return(20-(X1-X2)^2-8*(X1+X2-4)^3)
}

form<-PROB_MACHINE(name="FORM Rack.-Fieß.",fCall="FORM")

ps <- SYS_PROB(
  sys_input=list(lsf),
  probMachines = list(form)
)
ps$runMachines()
```

``` r
ps$beta_single
#>                  UQLab 2d_hat
#> FORM Rack.-Fieß.     3.434565
```

For more, check out the vignette! There are plenty more examples and
information about every single object.
