## Generate a test likelihood function. 
ll <- generateTestDensityMultiNormal(sigma = "no correlation")

## Create a BayesianSetup object from the likelihood 
## is the recommended way of using the runMCMC() function.
bayesianSetup <- createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))

## Finally we can run the sampler and have a look
settings = list(iterations = 1000, adapt = FALSE)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)

## We can plot the marginals in several ways:
## violin plots
marginalPlot(out, densityOnly = F) 
marginalPlot(out, densityOnly = F, priorTop = T, best = T) #flip plot
marginalPlot(out, densityOnly = F, plotPrior = FALSE) #plot only posterior

## histogram
marginalPlot(out, histogram = TRUE)
marginalPlot(out, histogram = TRUE, dens = F, priorTop = T)

## distributions
marginalPlot(out, priorTop = F) 
marginalPlot(out, plotPrior = FALSE)

##Further options
# We can pass arguments to getSample (check ?getSample)
marginalPlot(out, singlePanel = TRUE, scale=TRUE, 
             col = c("red", "blue"), thin = 100)
marginalPlot(out, singlePanel = TRUE, scale=TRUE, 
             col = c("red", "blue"), thin = 100, start = 500)