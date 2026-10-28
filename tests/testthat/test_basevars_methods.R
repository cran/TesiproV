test_that("fromJCSS methods generates Basevar correctly (fc, C30)", {
  a <- PROB_BASEVAR(Id = 1, Name = "a", DistributionType = "lnorm", Mean = 10, Sd = 0.1)
  a <- a$fromJCSS("fc", "C30", "ready-mixed")

  a_2 <- PROB_BASEVAR(Id = 1, Name = "a", DistributionType = "lt", DistributionParameters = c(3.75, 0.105, 3.0, 10))

  # Compare the fields, expect_equal() on the objects also compares the
  # environment including already used methods (edition 3)
  expect_equal(a$Id, a_2$Id)
  expect_equal(a$Name, a_2$Name)
  expect_equal(a$DistributionType, a_2$DistributionType)
  expect_equal(a$DistributionParameters, a_2$DistributionParameters)
  expect_equal(a$Package, a_2$Package)
})

test_that("fromJCSS methods generates Basevar correctly (fy, B500)", {
  a <- PROB_BASEVAR(Id = 1, Name = "a", DistributionType = "lnorm", Mean = 10, Sd = 0.1)
  a <- a$fromJCSS("fy", "B500", "reinforcing steel")

  a_2 <- PROB_BASEVAR(Id = 1, Name = "a", DistributionType = "norm", Mean = 560, Sd = 30)

  expect_equal(a$Id, a_2$Id)
  expect_equal(a$Name, a_2$Name)
  expect_equal(a$DistributionType, a_2$DistributionType)
  expect_equal(a$DistributionParameters, a_2$DistributionParameters)
  expect_equal(a$Package, a_2$Package)
  expect_equal(a$Mean, a_2$Mean)
  expect_equal(a$Sd, a_2$Sd)
})

test_that("genHyperParam methods generates Basevar correctly", {
  a <- PROB_BASEVAR(Id = 1, Name = "a", DistributionType = "gumbel", Mean = 10, Sd = 1)

  expect_error(a$genHyperParam(), "only implemented")

  b <- PROB_BASEVAR(Id = 2, Name = "b", DistributionType = "norm", Mean = 10, Sd = 1)
  b <- b$genHyperParam()
  expect_equal(length(b$DistributionParameters), 4)
  expect_equal(b$Package, "TesiproV")
})
