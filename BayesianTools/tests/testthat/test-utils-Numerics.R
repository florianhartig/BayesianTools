context("Utils Numerics work")

skip_on_cran()

test_that("logSumExp", {
  
  x = rep(log(1),10)
  expect_equal(exp(BayesianTools:::logSumExp(x)), 10)
  expect_equal(exp(BayesianTools:::logSumExp(x,mean = T)), 1)  
  
  x = c(rep(log(1),9), -Inf)  
  expect_equal(exp(BayesianTools:::logSumExp(x)), 9)
  expect_equal(exp(BayesianTools:::logSumExp(x, mean = T)), 9/10)

  x = rep(-100, 100) 
  expect_equal(log(sum(exp(x))),   BayesianTools:::logSumExp(x))

})

