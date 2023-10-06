

# likelihood function can be created with or without createBayesianSetup. ??createBayesianSetup provides additional properties.


data = rnorm(20)

# create standardized likelihood function for Gaussian likelihood function
likelihood <- function(x) likelihoodIidNormal(rnorm(20, mean = x), data, 2)
likelihoodClass <- createLikelihood(likelihood)
bayesianSetup <- createBayesianSetup(likelihood = likelihoodClass, lower = -10, upper = 10)
settings = list(iterations = 100)
plotSensitivity(bayesianSetup)
out <- runMCMC(bayesianSetup, sampler = "DEzs", settings)


# create likelihood function for AR1 likelihood function
likelihood <- function(x) likelihoodAR1(rnorm(20, mean = x), data, 2, 0.6)
bayesianSetup <- createBayesianSetup(likelihood = likelihood, lower = -10, upper = 10)
settings = list(iterations = 100)
out <- runMCMC(bayesianSetup, sampler = "DEzs", settings)


# create likelihood function for banana likelihood function
likelihood <- testDensityBanana
bayesianSetup <- createBayesianSetup(likelihood = likelihood, lower = rep(-10,2), upper = rep(10,2))
settings = list(iterations = 100)
out <- runMCMC(bayesianSetup, sampler = "DEzs", settings)

