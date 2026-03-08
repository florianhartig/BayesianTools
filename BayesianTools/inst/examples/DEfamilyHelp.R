library(BayesianTools)

ll <- generateTestDensityMultiNormal(sigma = "no correlation")
bayesianSetup <- createBayesianSetup(likelihood = ll, 
                                     lower = rep(-10, 3), 
                                     upper = rep(10, 3))

settings = list(iterations = 200)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
summary(out)

# DE family samplers are population MCMCs that run a number of internal chains
# in parallel. Here examples how to change the internal chains
# note that internal chains can be executedi n parallel
settings = list(startValue = 4, iterations = 200)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
summary(out)

# Modify the start values of the internal chains (note that this is a matrix
# of dim nChain * nPar)
settings = list(startValue = matrix(rnorm(12), nrow = 4, ncol = 3), 
                iterations = 200)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
summary(out)

# In the DE sampler family with Z matrix, the previous chains are written in 
# a common matrix, from which proposals are generated. Per default this matrix
# is started with samples from the prior, but we can change this. Often useful
# to improve sampler convergence, 
# see  https://github.com/florianhartig/BayesianTools/issues/79
settings = list(startValue = matrix(rnorm(12), nrow = 4, ncol = 3),
                Z = matrix(rnorm(300), nrow = 100, ncol = 3),
                iterations = 200)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
summary(out)


