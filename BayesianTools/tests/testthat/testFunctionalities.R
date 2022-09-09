context("Test basic functionalities")

set.seed(1)
library(BayesianTools)
library(testthat)


settings = list(iterations = 100)
ll <- function(x) sum(dnorm(x, log = T))

testFunctions <-function(x){
  print(x)
  summary(x)
  plot(x)
  marginalPlot(x)
  getSample(x)
  DIC(x)
  MAP(x)
  getVolume(x)
  marginalLikelihood(x)  
}


# 1d par, 1d / 3d sampler
bayesianSetup <- createBayesianSetup(likelihood = ll, lower = rep(-10, 1), upper = rep(10, 1))

out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
testFunctions(out)

out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
testFunctions(out)


# 3d par, 1d / 3d sampler
bayesianSetup <- createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))

out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
testFunctions(out)

out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
testFunctions(out)


