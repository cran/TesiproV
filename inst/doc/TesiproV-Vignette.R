## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
#install.package("TesiproV")
library(TesiproV)

## ---- eval=FALSE--------------------------------------------------------------
#  PROB_BASEVAR(Id, Name, Description, Package, DistributionType, DistributionParameters, Mean, Sd, Cov, X0)

## -----------------------------------------------------------------------------
var.f_ck <- PROB_BASEVAR(Id=1,Name="f_ck",Description="Char. Betondruckfestigkeit",Package="TesiproV",                         DistributionType="lt",DistributionParameters=c(3.85,0.09,3,10),Mean=47.35,Sd=5.86)

## ----eval= FALSE--------------------------------------------------------------
#  PARAM_BASEVAR(Id, Name, Description, Package, DistributionType, Sd, ParamType, ParamValues)

## ---- eval=FALSE--------------------------------------------------------------
#  var.d <- PARAM_BASEVAR(Id=1,Name="d",Description="Stat. Nutzhöhe",DistributionType="norm",Sd=10,ParamType="Mean",ParamValues= c(seq(200,950,50),seq(1000,3000,100))+10)

## -----------------------------------------------------------------------------
var.gamma_c <- PROB_DETVAR(Id=1,Name="gamma_c",Description="TSB Beton", Value=1.5)

## -----------------------------------------------------------------------------
var.V_Ed <- PARAM_DETVAR(Id=1,Name="V_Ed",Description="Einwirkende Querkraft", ParamValues=c(1,2,3,4,seq(5,10,0.5)))

## -----------------------------------------------------------------------------
var1 <- PROB_BASEVAR(Name="E",Description="Effect",DistributionType="norm",Mean=10,Sd=1)
var2 <- PROB_BASEVAR(Name="R",Description="Resistance",DistributionType="norm",Mean=15,Sd=1.5)


lsf<-SYS_LSF(vars=list(var1,var2),name="LSF XY")
lsf$func <- function(E,R){R-E}
# you can run this check to see if all transformations of the variables worked, all needed data is available and the set of vars fits to the limit state function
lsf$check()


## -----------------------------------------------------------------------------
machine <- PROB_MACHINE(name="MVFOSM",fCall="MVFOSM",options=list("isExpression"= FALSE,"h"=0.1))

## -----------------------------------------------------------------------------
machine <- PROB_MACHINE(name="FORM Rack.-Fieß.",fCall="FORM",options=list("n_optim"=20, "loctol"=0.001, "optim_type"="rackfies"))

machine <- PROB_MACHINE(name="FORM Lagrange.",fCall="FORM",options=list("optim_type"="auglag"))

## -----------------------------------------------------------------------------
machine <- PROB_MACHINE(name="SORM",fCall="SORM")

## -----------------------------------------------------------------------------
machine <- PROB_MACHINE(name="MC CoV 0.05",fCall="MC_CRUDE",options=list("n_max"=1e6, "cov_user"=0.05))

## -----------------------------------------------------------------------------
machine <- PROB_MACHINE(name="MC IS",fCall="MC_IS",options=list("cov_user" = 0.05, "n_max"=300000))

## -----------------------------------------------------------------------------
machine <- PROB_MACHINE(name="MC Subset Sampling",fCall="MC_SubSam",options=list("variance"="uniform"))

## ---- results="hide"----------------------------------------------------------
var1 <- PROB_BASEVAR(Name="E",Description="Effect",DistributionType="norm",Mean=10,Sd=1)
var2 <- PROB_BASEVAR(Name="R",Description="Resistance",DistributionType="norm",Mean=15,Sd=1.5)


lsf<-SYS_LSF(vars=list(var1,var2),name="LSF XY")
lsf$func <- function(E,R){R-E}
# you can run this check to see if all transformations of the variables worked, all needed data is available and the set of vars fits to the limit state function
lsf$check()


machine <- PROB_MACHINE(name="FORM Rack.-Fieß.",fCall="FORM",options=list("n_optim"=20, "loctol"=0.001, "optim_type"="rackfies"))

ps <- SYS_PROB(
  sys_input=list(lsf),
  probMachines = list(machine),
  debug.level=0
)

ps$runMachines()


# Systemberechnungn (System aus zwei gleichen LSF´s)
ps2 <- SYS_PROB(
  sys_input=list(lsf,lsf),
  probMachines = list(machine),
  sys_type="serial"
)

ps2$runMachines()
ps2$calculateSystemProbability()

#ps2$calculateSystemProbability("MCSUS")

#params <- list("cov_user"=0.10)
#ps2$calculateSystemProbability("MCC",params)
#ps2$calculateSystemProbability("MCIS",params)

## ---- eval=FALSE--------------------------------------------------------------
#  var1 <- PARAM_BASEVAR(Name="E",Description="Effect",DistributionType="Norm",Sd=1,ParamType="Mean",Values=c(9,11,0.1))
#  
#  #[...]
#  
#  ps <- SYS_PARAM(
#    sys_input=list(lsf),
#    probMachines = list(form)
#  )
#  
#  ps$runMachines()

## ---- eval=FALSE--------------------------------------------------------------
#  ps <- SYS_PARAM(
#    sys_input=list(lsf),
#    probMachines = list(form)
#  )
#  
#  ps$runMachines()
#  
#  local_beta <- ps$beta
#  local_runtime <- ps$res_single[[1]][[1]]$runtime
#  
#  #define path with setwd("...")
#  ps$printResults("result_file") #Prints a report file in .txt format with all informations about the calculation
#  ps$saveProject(4,"project_file") #stores the project in a file
#  

## ---- results="hide"----------------------------------------------------------
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

## -----------------------------------------------------------------------------
ps$beta_single

## ---- results ='hide'---------------------------------------------------------
library("evd")
X1 <- PROB_BASEVAR(Id=1,Name="Z",DistributionType="norm",Mean=100,Cov=0.04) #kN/m²
X2 <- PROB_BASEVAR(Id=2,Name="Fy",DistributionType="lnorm",Mean=40,Cov=0.1) #m
X3 <- PROB_BASEVAR(Id=2,Name="M",DistributionType="gumbel",Package="evd",Mean=2000,Cov=0.1) #m

lsf<-SYS_LSF(vars=list(X1,X2,X3),name="Nowak & Collins Exp5.11_1")
lsf$func <- function(M,Fy,Z){
  K <- Z*Fy
  return(0.9*K-M)
}

form<-PROB_MACHINE(name="FORM Rack.-Fieß.",fCall="FORM")

ps <- SYS_PROB(
  sys_input=list(lsf),
  probMachines = list(form)
)
ps$runMachines()

## -----------------------------------------------------------------------------
ps$beta_single

## ---- results="hide"----------------------------------------------------------

Var1 <- PROB_BASEVAR(Id=1,Name="X1",DistributionType="norm",Mean=0.25,Sd=1) #kN/m²
Var2 <- PROB_BASEVAR(Id=2,Name="X2",DistributionType="norm",Mean=0.25,Sd=1) #m

lsf1<-SYS_LSF(vars=list(Var1,Var2),name="UQLab 2d_hat")
lsf1$func <- function(X1,X2){
  return(20-(X1-X2)^2-8*(X1+X2-4)^3)
}


X1 <- PROB_BASEVAR(Id=1,Name="Z",DistributionType="norm",Mean=100,Cov=0.04) #kN/m²
X2 <- PROB_BASEVAR(Id=2,Name="Fy",DistributionType="lnorm",Mean=40,Cov=0.1) #m
X3 <- PROB_BASEVAR(Id=2,Name="M",DistributionType="gumbel",Package="evd",Mean=2000,Cov=0.1) #m

lsf2<-SYS_LSF(vars=list(X1,X2,X3),name="Nowak & Collins Exp5.11_1")
lsf2$func <- function(M,Fy,Z){
  K <- Z*Fy
  return(0.75*K-M)
}

form<-PROB_MACHINE(name="FORM Rack.-Fieß.",fCall="FORM")

ps <- SYS_PROB(
  sys_input =list(lsf1,lsf2),
  sys_type = "serial",
  probMachines = list(form)
)
ps$runMachines()


## -----------------------------------------------------------------------------
ps$beta_single

ps$calculateSystemProbability()
ps$beta_sys

## ---- eval=FALSE--------------------------------------------------------------
#  params <- list("cov_user"=0.05,"n_max"=50000)
#  ps$calculateSystemProbability("MC",params)
#  
#  params <- list("cov_user"=0.05,"n_max"=50000)
#  ps$calculateSystemProbability("MCIS",params)
#  
#  ps$calculateSystemProbability("MCSUS",params)

## ---- results="hide", warning=FALSE-------------------------------------------
var.d <- PARAM_BASEVAR(Id=1,Name="d",Description="Stat. Nutzhöhe",DistributionType="norm",Sd=10,ParamType="Mean",ParamValues= c(seq(200,800,50),seq(1000,3000,100))+10)
pre.d <- c(seq(200,800,50),seq(1000,3000,100))

var.f_ck <- PROB_BASEVAR(Id=1,Name="f_ck",Description="Char.Betondruckfestigkeit",Package="TesiproV",DistributionType="lt",DistributionParameters=c(3.85,0.09,3,10),Mean=47.35,Sd=5.86)
pre.f_ck <- 35

var.f_ywk <- PROB_BASEVAR(Id=1,Name="f_ywk",Description="Zugfestigkeit Bügelbewehrung",DistributionType="norm",Mean=560,Cov=0.02)
pre.f_ywk <- 500

var.f_yk <- PROB_BASEVAR(Id=1,Name="f_yk",Description="Zugfestigkeit Längsbewehrung",DistributionType="norm",Mean=560,Sd=30)
pre.f_yk <- 500

var.rho_l <- PROB_BASEVAR(Id=1,Name="rho_l",Description="Biegebewehrungsgrad",DistributionType="norm",Mean=0.01,Cov=0.02)
pre.rho_l <- var.rho_l$Mean

var.c_D <- PROB_BASEVAR(Id=1,Name="c_D",Description="Stützendurchmesser",DistributionType="norm",Mean=500+0.003*500,Sd=4+0.006*500)
pre.c_D <- 500

var.theta_c <- PROB_BASEVAR(Id=1,Name="theta_c",Description="Modellunsicherheit c nach DIBT",DistributionType="lnorm", Mean=1.0644,Cov=0.1679)

var.gamma_c <- PROB_DETVAR(Id=1,Name="gamma_c",Description="TSB Beton", Value=1.5)
var.gamma_s <- PROB_DETVAR(Id=1,Name="gamma_s",Description="TSB Stahl", Value=1.15)
var.alpha_cc <- PROB_DETVAR(Id=1,Name="alpha_cc",Description="Langzeitfaktor AlphaCC", Value=0.85)

V_Ed <- c(0.6412447,0.8760515,1.1424185,1.4399554,1.7554079,2.0188587,2.3003489,2.5997038,2.9167825,3.2514671,3.6036567,3.9732630,
          4.3602076,6.0800504,7.0424326,8.0726216,9.1702994,10.3351829,11.5670173,12.8655715,14.2306347,15.6620132,17.1595283,18.7230145,
          20.3523175,22.0472935,23.8078076,25.6337332,27.5249511,29.4813487,31.5028193,33.5892621,35.7405811,37.9566849)
var.V_Ed <- PARAM_DETVAR(Id=1,Name="V_Ed",Description="Einwirkende Querkraft", ParamValues=unlist(V_Ed))

lsf1 <-SYS_LSF(vars=list(var.d,var.f_ck,var.f_yk,var.f_ywk,var.rho_l,var.c_D,
                         var.theta_c,var.V_Ed,var.gamma_c,var.gamma_s,var.alpha_cc),
               name="Durchstanzwiderstand ohne Bewehrung")
lsf1$func <- function(d, f_ck,f_yk,f_ywk, rho_l, c_D, theta_c, V_Ed, gamma_c, gamma_s, alpha_cc){

  f_cd <- f_ck*alpha_cc/gamma_c
  f_yd <- f_yk/gamma_s
  f_ywd <- f_ywk/gamma_s

  u_0 <- c_D*pi
  u_1 <- (c_D/2+2*d)*2*pi
  u_0.5 <- (c_D/2+0.5*d)*2*pi

  k <- min(1+sqrt(200/d),2)

  if((u_0/d)<4){
    C_Rd_c <- 0.18 * (0.1*u_0/d + 0.6)
  }else{
    C_Rd_c <- 0.18
  }
  rho_l <- min(rho_l,0.5*f_cd/f_yd,0.02)
  C_Rd_c1 <- C_Rd_c * k * (100*rho_l*f_ck)^(1/3)

  v_Rd_c <- C_Rd_c1
  V_Rd_c <- v_Rd_c*d*u_1*1e-6*theta_c
  return(V_Rd_c-V_Ed)
}

form <- PROB_MACHINE(name="FORM RF.",fCall="FORM",options=list("n_optim"=20, "loctol"=0.0001, "optim_type"="rackfies"))

ps <- SYS_PARAM(
  sys_input=list(lsf1),
  probMachines = list(form)
)
ps$runMachines()

## -----------------------------------------------------------------------------
d <- pre.d
betas <- unlist(unname(ps$beta_params))
plot(x=d,y=betas,type="l")

