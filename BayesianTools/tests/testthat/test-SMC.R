context("Test SMC and utility functions") # Must be at the top of the file (https://stackoverflow.com/questions/50083521/error-in-xmethod-attempt-to-apply-non-function-in-testthat-test-when)

library(testthat)
library(BayesianTools)
library(Matching)

test <- "exact"

#skip_on_cran()



set.seed(123)



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
  particles <- matrix(particles, ncol=1)
  importanceDensity <- function(x){return(dunif(x, min=-10, max=10, log = TRUE))}
  
  ### Testing *Metropolis* implementation
  
  posteriorValues <- setup$posterior$density(particles)
  # Initialize proposal Generator
  proposalGenerator <- createProposalGenerator(covariance=var(particles))
  mutate.test <- mutate(setup = setup, particles = particles, proposalGenerator = proposalGenerator, posteriorValues = posteriorValues, importanceDensity = importanceDensity, method="Metropolis", steps=50, proposalScale = 0.5)
  
  # Outputs should be normally distributed after a couple of steps
  ks.pval <- ks.test(mutate.test$particles, rnorm(n=length(test.smc$particles)))$p.value
  expect_gt(ks.pval,0.05)
  
  ### Testing *Differential Evolution* implementation
  
  mutate.test <- mutate(setup = setup, particles = particles, proposalGenerator = proposalGenerator, posteriorValues = posteriorValues, importanceDensity = importanceDensity, method="DE", steps=50, proposalScale = 0.5)
  
  # Outputs should be normally distributed after a couple of steps
  ks.pval <- ks.test(mutate.test$particles, rnorm(n=length(test.smc$particles)))$p.value
  expect_gt(ks.pval,0.05)
})

test_that("univariate SMC sampler works", {
  ll = function(x) sum(dnorm(x, log = T))
  setup = createBayesianSetup(ll, lower = c(-10), upper = c(10))
  
  initialParticles <- list()
  initialParticles$particles <- runif(n=5000, min=-10, max=10)
  initialParticles$density <- function(x){return(dunif(x, min=-10, max=10, log = TRUE))}
  
  # Test with *Metropolis* mutation
  test.smc <- smcSampler(setup, initialParticles = initialParticles, mutate.method ="Metropolis",sampling="systematic",adaptive=TRUE, ess.factor=0.99, resamplingSteps = 1)
  ks.pval <- ks.boot(rnorm(n=length(test.smc$particles)), test.smc$particles)$ks.boot.pvalue
  expect_gt(ks.pval,0.05)
  
  # Test with *Differential Evolution* mutation
  test.smc <- smcSampler(setup, initialParticles = initialParticles, mutate.method ="DE",sampling="systematic",adaptive=TRUE, ess.factor=0.99, resamplingSteps = 1)
  ks.pval <- ks.boot(rnorm(n=length(test.smc$particles)), test.smc$particles)$ks.boot.pvalue
  expect_gt(ks.pval,0.05)
})

test_that("2D sampler works", {
  ll = function(x) sum(dnorm(x, log = T))
  setup = createBayesianSetup(ll, lower = c(-10,-10), upper = c(10,10))
  
  initialParticles <- list()
  part1 <- runif(n=5000, min=-10, max=10)
  part2 <- runif(n=5000, min=-10, max=10)
  initialParticles$particles <- data.frame(part1,part2)
  initialParticles$density <- function(x){return(rep(dunif(1, min=-10, max=10, log = TRUE), nrow(x)))}
  
  ### Test with *Metropolis* mutation
  test.smc <- smcSampler(setup, initialParticles = initialParticles, mutate.method ="Metropolis",sampling="systematic",adaptive=TRUE, ess.factor=0.99, resamplingSteps = 1)
  
  #ks.pval1 <- ks.test(test.smc$particles[,1], rnorm(n=length(test.smc$particles)))$p.value
  ks.pval1 <- ks.boot(rnorm(n=length(test.smc$particles[,1])), test.smc$particles[,1])$ks.boot.pvalue
  expect_gt(ks.pval1,0.05)

  ks.pval2 <- ks.boot(rnorm(n=length(test.smc$particles[,2])), test.smc$particles[,2])$ks.boot.pvalue
  expect_gt(ks.pval2,0.05)
  
  ### Test with *Differential Evolution* mutation
  test.smc <- smcSampler(setup, initialParticles = initialParticles, mutate.method ="DE",sampling="systematic",adaptive=TRUE, ess.factor=0.99, resamplingSteps = 1)
  
  ks.pval1 <- ks.boot(rnorm(n=length(test.smc$particles[,1])), test.smc$particles[,1])$ks.boot.pvalue
  expect_gt(ks.pval1,0.05)
  
  ks.pval2 <- ks.boot(rnorm(n=length(test.smc$particles[,2])), test.smc$particles[,2])$ks.boot.pvalue
  expect_gt(ks.pval2,0.05)
})
