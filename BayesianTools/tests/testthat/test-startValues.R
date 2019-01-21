context("Test start values/Z matrix inside prior range")

skip_on_cran()

test_that("test z matrix check", {
  ll = testDensityNormal
  bayesianSetup = createBayesianSetup(likelihood = ll, lower = 2, upper = 5)
  expect_warning(runMCMC(bayesianSetup, sampler = "DEzs",settings = list(iterations = 100, Z = matrix(6, nrow = 3, ncol =1), message = F)))
  expect_warning(runMCMC(bayesianSetup, sampler = "DREAMzs",settings = list(iterations = 100, Z = matrix(6, nrow = 3, ncol =1), message = F)))
  
})

test_that("test start value check", {
  ll = testDensityNormal
  bayesianSetup = createBayesianSetup(likelihood = ll, lower = 2, upper = 5)
  expect_warning(runMCMC(bayesianSetup, sampler = "DREAMzs",settings = list(iterations = 100, startValue = matrix(6, nrow = 3, ncol =1), message = F)))
  expect_warning(runMCMC(bayesianSetup, sampler = "DREAM",settings = list(iterations = 100, startValue = matrix(1:4, nrow = 4, ncol =1), message = F, DEpairs = 2)))
  expect_warning(runMCMC(bayesianSetup, sampler = "DEzs",settings = list(iterations = 100, startValue = matrix(6, nrow = 3, ncol =1), message = F)))
  expect_warning(runMCMC(bayesianSetup, sampler = "DE",settings = list(iterations = 100, startValue = matrix(6, nrow = 3, ncol =1), message = F)))
  expect_warning(runMCMC(bayesianSetup, sampler = "Metropolis",settings = list(iterations = 100, startValue = matrix(6, nrow = 1, ncol =1), message = F)))
  
  density = function(par) return(dunif(par[1],0,5, log = T))
  sampler = function(n=1) return(runif(n,-1,4))
  prior = createPrior(sampler = sampler, density = density, lower = 0,upper = 5)
  bayesianSetup = createBayesianSetup(likelihood = ll, prior = prior)
  expect_warning(runMCMC(bayesianSetup, sampler = "DE",settings = list(iterations = 100, message = F)))
})
