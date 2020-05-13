library(BayesianTools)

ll <- generateTestDensityMultiNormal(sigma = "no correlation")
bayesianSetup <- createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))

settings = list(iterations = 200)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
summary(out)

# all DE family samplers have a number of internal chains that can be
# parallelized. Here examples how to change the internal chains
settings = list(startValue = 4, iterations = 200)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
summary(out)

# Start with explicit start values
settings = list(startValue = matrix(rnorm(12), nrow = 4, ncol = 3), 
                iterations = 200)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
summary(out)

# Start with explicit start values for the Z matrix - this is often useful
# to improve sampler convergence, 
# see also https://github.com/florianhartig/BayesianTools/issues/79
settings = list(startValue = matrix(rnorm(12), nrow = 4, ncol = 3),
                Z = matrix(rnorm(300), nrow = 100, ncol = 3),
                iterations = 200)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
summary(out)


