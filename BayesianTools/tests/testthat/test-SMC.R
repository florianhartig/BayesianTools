skip_on_cran()

context("Test SMC and utility functions")





test_that("resampling works", {
  
  weights = log(c(1,2,3)/6) 
  size = 10000
  
  x1 = replicate(size, resample(weights))
  dif = table(x1) - exp(weights) * size * 3
  expect_lt(sum(abs(dif)), 300)
  
  x2 = replicate(size, resample(weights, method = "residual"))
  dif = table(x2) - exp(weights) * size * 3
  expect_lt(sum(abs(dif)), 300)
  
  x3 = replicate(size, resample(weights, method = "systematic"))
  dif = table(x3) - exp(weights) * size  * 3
  expect_lt(sum(abs(dif)), 300)
  
}
