context("Test SMC and utility functions") # Must be at the top of the file (https://stackoverflow.com/questions/50083521/error-in-xmethod-attempt-to-apply-non-function-in-testthat-test-when)

library(testthat)
library(BayesianTools)

test <- "exact"

#skip_on_cran()



set.seed(3)



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
  
})

test_that("adaptive scheme works", {
  weights <- log(rep(1/100, 100))
  ess <- 1 / sum(exp(2 * weights))
  posteriorValues <- log(rnorm(100,10,1))
  posteriorValues <- posteriorValues - BayesianTools:::logSumExp(posteriorValues)
  importanceValues <- log(runif(100,0,1))
  
  ess.factor <- 0.95
  
  # Initial iteration: start with uniform importance distribution
  curExp <- 0
  interDist <- importanceValues
  
  inter.out <- beta.search(ess=ess, target.ess=ess.factor * ess, posteriorValues = posteriorValues, importanceValues = importanceValues, oldInter = interDist, curWeights = weights, curExp = curExp, tol=1)
  
  curExp <- inter.out$newExp
  weights <- inter.out$weights
  interDist <- inter.out$interDist
  new.ess <- inter.out$ess
  
  # First iteration: tests whether new ESS is less than original one, and within the tolerance limit specified as argument to beta.search (tol = 1)
  
  expect_lt(new.ess,ess)
  expect_gt(new.ess,(ess*ess.factor - 1))
  
  inter.out <- beta.search(ess=new.ess, target.ess=ess.factor * ess, posteriorValues = posteriorValues, importanceValues = importanceValues, oldInter = interDist, curWeights = weights, curExp = curExp, tol=1)
  new.new.ess <- inter.out$ess
  
  # Second iteration: tests whether new ESS is less than previous one, and within the tolerance limit specified as argument to beta.search (tol = 1)
  
  expect_lt(new.new.ess,new.ess)
  expect_gt(new.new.ess,(ess*ess.factor^2 - 1))
})

test_that("mutation works", {
  
  # Create test BayesianSetup
  ll = function(x) sum(dnorm(x, log = T))
  setup = createBayesianSetup(ll, lower = c(-10), upper = c(10))
  
  particles <- runif(n=5000, min=-10, max=10)
  importanceDensity <- function(x){return(dunif(x, min=-10, max=10, log = TRUE))}
  
  posteriorValues <- setup$posterior$density(particles)
  
  # Initialize proposal Generator
  proposalGenerator <- createProposalGenerator()
  
})