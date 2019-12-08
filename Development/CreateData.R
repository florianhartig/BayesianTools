# ensure WD is BayesianTools package

library(BayesianTools)
library(devtools)
library(usethis)

use_data_raw()

ll <- generateTestDensityMultiNormal(sigma = "no correlation")
bayesianSetup = createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))
settings = list(iterations = 10000, nrChains = 3)
mcmcOutput <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
use_data(mcmcOutput, overwrite = TRUE, name = "mcmcRun", internal = T)
