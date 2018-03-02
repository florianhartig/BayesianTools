# Create bayesian setup with  
bayesianSetup <- createBayesianSetup(likelihood = testDensityNormal, 
                                     prior = createUniformPrior(lower = rep(-10,2),
                                                                upper = rep(10,2)))

# running MCMC

out = runMCMC(bayesianSetup = bayesianSetup)

# diagnostic plots

plotDiagnostic(out)

