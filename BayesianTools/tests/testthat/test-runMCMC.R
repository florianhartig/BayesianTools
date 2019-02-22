context("Test runMCMC interface")

set.seed(1)
library(BayesianTools)
library(testthat)


ll <- function(x) sum(dnorm(x, log = T))
bayesianSetup <- createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))

settings = list(iterations = 200)
out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)
expect_error(runMCMC(bayesianSetup = ll, settings = settings))



