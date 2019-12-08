## code to prepare `DATASET` dataset goes here

# note limitations of use_data https://github.com/r-lib/usethis/issues/789

# usethis::proj_set("./BayesianTools/") does not work

# workaround is to
# create a new project in the package directory (I had my project one level above the package directory)
# use the function outside an R project

library(BayesianTools)

ll <- generateTestDensityMultiNormal(sigma = "no correlation")
bayesianSetup = createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))
settings = list(iterations = 10000, nrChains = 3)
mcmcOutput <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)

#use_data(mcmcOutput, overwrite = TRUE, name = "mcmcRun", internal = T)

usethis::use_data("exampleMCMC", overwrite = TRUE)


