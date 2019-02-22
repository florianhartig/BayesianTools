context("Test runMCMC interface")

set.seed(1)
library(BayesianTools)
library(testthat)


ll <- function(x) sum(dnorm(x, log = T))
bayesianSetup <- createBayesianSetup(likelihood = ll, lower = rep(-10, 4), upper = rep(10, 4))

settings = list(iterations = 5)
out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)

x = gelmanDiagnostics(out)
summary(out)


settings = list(iterations = 10, nrChains = 3)
out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)
gelmanDiagnostics(out)
summary(out)


expect_error(runMCMC(bayesianSetup = ll, settings = settings))



