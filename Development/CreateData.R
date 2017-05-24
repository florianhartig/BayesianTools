# ensure WD is BayesianTools package

library(BayesianTools)
ll <- generateTestDensityMultiNormal(sigma = "no correlation")
bayesianSetup = createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))
settings = list(iterations = 3000)
mcmcOutput <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
devtools::use_data(mcmcOutput, overwrite = TRUE)


