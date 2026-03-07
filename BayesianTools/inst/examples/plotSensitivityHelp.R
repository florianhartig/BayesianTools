
ll <- testDensityBanana
bayesianSetup <- createBayesianSetup(likelihood = ll, lower = rep(-10, 2), upper = rep(10, 2))

plotSensitivity(bayesianSetup)

# can select parameters
plotSensitivity(bayesianSetup, selection = c(2,1))

# scale each parameter independently
plotSensitivity(bayesianSetup, equalScale = FALSE)
