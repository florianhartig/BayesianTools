bayesianSetup = createBayesianSetup(likelihood = generateTestDensityMultiNormal(sigma = "no correlation"), lower = rep(-10, 3), upper = rep(10, 3))

out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = list(iterations = 2000, message = FALSE))

getVolume(out, prior = T)

bayesianSetup = createBayesianSetup(likelihood = generateTestDensityMultiNormal(sigma = "strongcorrelation"), lower = rep(-10, 3), upper = rep(10, 3))

out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = list(iterations = 2000, message = FALSE))

getVolume(out, prior = T)

