# Create bayesian setup with  
bayesianSetup <- createBayesianSetup(likelihood = testDensityNormal, 
                                     prior = createUniformPrior(lower = -10,
                                                                upper = 10))

# running MCMC

out = runMCMC(bayesianSetup = bayesianSetup)

# diagnostic plots
\dontrun{
plotDiagnostic(out)
}
