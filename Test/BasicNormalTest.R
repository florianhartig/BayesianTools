library(BayesianTools)

ll <- generateTestDensityMultiNormal(sigma = "no correlation")

bayesianSetup = createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))



settings = list(iterations = 5000, nrChains= 3, message = FALSE, startValue = list(rep(-9, 3), rep(0, 3), rep(9, 3)), optimize = F)

settings <- list(iterations = iter, adapt = T, DRlevels = 2, gibbsProbabilities = NULL, temperingFunction = NULL, optimize = T,nrChains= 1,burnin=7000,adaptationInterval=10)

settings <- list(iterations = 10000, adapt = T, DRlevels = 2 , gibbsProbabilities = NULL, temperingFunction = NULL, optimize = T, nrChains= 1, burnin=2000, adaptationInterval=100)

out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)


plot(out)


getSeett



applySettingsDefault(sampler = "Metropolis")
