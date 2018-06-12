## Generate a test likelihood function. 
ll <- generateTestDensityMultiNormal(sigma = "no correlation")

## Create a BayesianSetup object from the likelihood 
## is the recommended way of using the runMCMC() function.
bayesianSetup <- createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))

## Finally we can run the sampler and have a look
settings = list(iterations = 1000)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)

## Correlation density plots:
correlationPlot(out)
correlationPlot(out, scaleCorText = FALSE) # text to display correlation will be not scaled to the strength of the correlation
## We can also switch the method for calculating correllations
correlationPlot(out, scaleCorText = FALSE, method = "spearman")
correlationPlot(out, scaleCorText = FALSE, method = "kendall")

## Select parameters that should be plotted (by indices)
correlationPlot(out, scaleCorText = FALSE, method = "kendall", whichParameters = c(1,2))

## additional parameters can be passed to getSample (see ?getSample for further information)
correlationPlot(out, scaleCorText = FALSE, thin = 100, start = 200)