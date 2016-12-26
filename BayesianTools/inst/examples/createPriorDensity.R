# Create a BayesianSetup
ll <- generateTestDensityMultiNormal(sigma = "no correlation")
bayesianSetup = createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))

settings = list(iterations = 1000)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)

# now use the output as a new prior TODO
#bayesianSetup = createBayesianSetup(likelihood = ll, prior = out)
#out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)

