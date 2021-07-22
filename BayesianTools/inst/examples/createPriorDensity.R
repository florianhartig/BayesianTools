# Create a BayesianSetup
ll <- generateTestDensityMultiNormal(sigma = "no correlation")
bayesianSetup = createBayesianSetup(likelihood = ll, 
                                    lower = rep(-10, 3), 
                                    upper = rep(10, 3))

settings = list(iterations = 1000)
out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)


newPrior = createPriorDensity(out, method = "multivariate",
                              eps = 1e-10, lower = rep(-10, 3),
                              upper =  rep(10, 3), best = NULL)

bayesianSetup <- createBayesianSetup(likelihood = ll, prior = newPrior)

\dontrun{
  # first prior / posterior
  marginalPlot(out)
  
  # second prior / posterior
  settings = list(iterations = 1000)
  out2 <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)
  marginalPlot(out2)
  
  # change scaling of the prior with the scaling option
  newPrior = createPriorDensity(out, method = "multivariate",
                                eps = 1e-10, lower = rep(-10, 3),
                                upper =  rep(10, 3), best = NULL, scaling = 0.5)
  
  apply(newPrior$sampler(1000),2,sd)
}






