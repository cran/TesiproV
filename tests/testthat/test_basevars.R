test_that("Basis var, transformation", {

  a <- PROB_BASEVAR(Id=1,Name="a",DistributionType="norm", Mean=1, Sd=0.1)
  a$prepare()
  mean_ <- a$DistributionParameters[1]
  sd_ <- a$DistributionParameters[2]
  expect_equal(mean_, 1)
  expect_equal(sd_,0.1)

  b <- PROB_BASEVAR(Id=1,Name="b",DistributionType="lnorm", Mean=1, Sd=0.1)
  b$prepare()
  mean_ <- b$DistributionParameters[1]
  sd_ <- b$DistributionParameters[2]
  expect_equal(mean_,-0.004975165)
  expect_equal(sd_,0.099751345)

  b <- PROB_BASEVAR(Id=1,Name="c",DistributionType="lnorm",DistributionParameters=c(mean_, sd_))
  b$prepare()
  expect_equal(b$Mean,1)
  expect_equal(b$Sd,0.1)
})



test_that("Basis var, parameters", {

  params <- c(2,3,4,5)
  means <- round(c(0.69189874,1.098057,1.385982,1.609238,0.6918987),2)
  sds <- round(c(0.04996879,0.03332408,0.0249961,0.019998,0.04996879),4)
  a <- PARAM_BASEVAR(Id=1,Name="a",DistributionType="lnorm", ParamType="Mean", Sd=0.1, ParamValues=params)

  for(i in 1:5){

  mean_ <- round(a$DistributionParameters[1],2)
  sd_ <- round(a$DistributionParameters[2],4)
  expect_equal(mean_, means[i])
  expect_equal(sd_,sds[i])
  a$nextParam()
  }

})



