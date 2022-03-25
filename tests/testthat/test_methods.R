test_that("FORM and SORM", {

  X1 <- PROB_BASEVAR(Id=1,Name="Z",DistributionType="norm",Mean=100,Cov=0.04) #kN/m²
  X2 <- PROB_BASEVAR(Id=2,Name="Fy",DistributionType="lnorm",Mean=40,Cov=0.1) #m
  X3 <- PROB_BASEVAR(Id=2,Name="M",DistributionType="gumbel",Package="evd",Mean=2000,Cov=0.1) #m

  lsf1<-SYS_LSF(vars=list(X1,X2,X3),name="Nowak & Collins Exp5.11_1")
  lsf1$check()
  lsf1$func <- function(M,Fy,Z){
    K <- Z*Fy
    return(0.8*K-M)
  }

  form_rf<-PROB_MACHINE(name="FORM Rack.-Fieß.",fCall="FORM",options=list("n_optim"=20, "loctol"=0.001, "optim_type"="rackfies"))
  sorm <- PROB_MACHINE(name="SORM",fCall="SORM")
  mcis<-PROB_MACHINE(name="MC IS",fCall="MC_IS",options=list("cov_user" = 0.05, "n_max"=300000))
  mcsus<-PROB_MACHINE(name="MC SuS",fCall="MC_SubSam")


  ps <- SYS_PROB(
    sys_input=list(lsf1),
    probMachines = list(form_rf,sorm,mcis)
  )
  ps$runMachines()
  ps$beta_single

  # Function is linear, FORM and SORM should have the same outcome...
  expect_equal(ps$beta_single[1],ps$beta_single[2])

})


test_that("MCIS and high no of inputs", {
# Example 4
u1 <- PROB_BASEVAR(Id=1,Name="u1",DistributionType="norm",Mean=0.5,Sd=0.1)
u2 <- PROB_BASEVAR(Id=2,Name="u2",DistributionType="norm",Mean=0.2,Sd=0.1)
u3 <- PROB_BASEVAR(Id=3,Name="u3",DistributionType="norm",Mean=0.2,Sd=0.1)
u4 <- PROB_BASEVAR(Id=4,Name="u4",DistributionType="norm",Mean=0.2,Sd=0.1)
u5 <- PROB_BASEVAR(Id=5,Name="u5",DistributionType="norm",Mean=0.2,Sd=0.1)
u6 <- PROB_BASEVAR(Id=6,Name="u6",DistributionType="norm",Mean=0.2,Sd=0.1)
u7 <- PROB_BASEVAR(Id=7,Name="u7",DistributionType="norm",Mean=0.2,Sd=0.1)
u8 <- PROB_BASEVAR(Id=8,Name="u8",DistributionType="norm",Mean=0.2,Sd=0.1)
u9 <- PROB_BASEVAR(Id=9,Name="u9",DistributionType="norm",Mean=0.2,Sd=0.1)
u10 <- PROB_BASEVAR(Id=10,Name="u10",DistributionType="norm",Mean=0.2,Sd=0.1)
u11 <- PROB_BASEVAR(Id=11,Name="u11",DistributionType="norm",Mean=0.2,Sd=0.1)
u12 <- PROB_BASEVAR(Id=12,Name="u12",DistributionType="norm",Mean=0.2,Sd=0.1)
u13 <- PROB_BASEVAR(Id=13,Name="u13",DistributionType="norm",Mean=0.2,Sd=0.1)
u14 <- PROB_BASEVAR(Id=14,Name="u14",DistributionType="norm",Mean=0.2,Sd=0.1)
u15 <- PROB_BASEVAR(Id=15,Name="u15",DistributionType="norm",Mean=0.2,Sd=0.1)
u16 <- PROB_BASEVAR(Id=16,Name="u16",DistributionType="norm",Mean=0.2,Sd=0.1)
u17 <- PROB_BASEVAR(Id=17,Name="u17",DistributionType="norm",Mean=0.2,Sd=0.1)
u18 <- PROB_BASEVAR(Id=18,Name="u18",DistributionType="norm",Mean=0.2,Sd=0.1)
u19 <- PROB_BASEVAR(Id=19,Name="u19",DistributionType="norm",Mean=0.2,Sd=0.1)
u20 <- PROB_BASEVAR(Id=20,Name="u20",DistributionType="norm",Mean=0.2,Sd=0.1)
u21 <- PROB_BASEVAR(Id=21,Name="u21",DistributionType="norm",Mean=0.2,Sd=0.1)
u22 <- PROB_BASEVAR(Id=22,Name="u22",DistributionType="norm",Mean=0.2,Sd=0.1)
u23 <- PROB_BASEVAR(Id=23,Name="u23",DistributionType="norm",Mean=0.2,Sd=0.1)
u24 <- PROB_BASEVAR(Id=24,Name="u24",DistributionType="norm",Mean=0.2,Sd=0.1)
u25 <- PROB_BASEVAR(Id=25,Name="u25",DistributionType="norm",Mean=0.2,Sd=0.1)
u26 <- PROB_BASEVAR(Id=26,Name="u26",DistributionType="norm",Mean=0.2,Sd=0.1)
lsf4 <- SYS_LSF(vars=list(u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26), name="RackwitzProb.")
lsf4$check()
lsf4$func <- function(u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26){

  r <- u1-(((u2^2)/1)+((u3^2)/2)+((u4^2)/3)+((u5^2)/4)+((u6^2)/5)+((u7^2)/6)+((u8^2)/7)+((u9^2)/8)+((u10^2)/9)+((u11^2)/10)+((u12^2)/11)+((u13^2)/12)+((u14^2)/13)+((u15^2)/14)+((u16^2)/15)+((u17^2)/16)+((u18^2)/17)+((u19^2)/18)+((u20^2)/19)+((u21^2)/20)+((u22^2)/21)+((u23^2)/22)+((u24^2)/23)+((u25^2)/24)+((u26^2)/25))

  return(r)
}

mcis<-PROB_MACHINE(name="MC IS",fCall="MC_IS",options=list("cov_user" = 0.05, "n_max"=300000))

ps <- SYS_PROB(
  sys_input=list(lsf4),
  probMachines = list(mcis)
)
ps$runMachines()

expect_true((ps$beta_single-2.7)<0)

})
