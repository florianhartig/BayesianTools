library(BayesianTools)


## Generate a test likelihood function. 
ll <- dnorm
  
## Create a BayesianSetup object from the likelihood 
## is the recommended way of using the runMCMC() function.
bayesianSetup <- createBayesianSetup(likelihood = ll, lower = -10, upper = 10)

## Finally we can run the sampler and have a look
settings = list(iterations = 1000, adapt = FALSE, nrChains = 2)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)

summary(out)


## out is of class bayesianOutput. There are various standard functions 
# implemented for this output

plot(out)
correlationPlot(out)
marginalPlot(out)

gelmanDiagnostics(out)

summary(out)


x = getSample(out, coda = T)
plot(x)

re

coda::gelman.diag(x)
coda::rejectionRate(x[[2]])




ll <- generateTestDensityMultiNormal(sigma = "no correlation")

bayesianSetup <- createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))

settings = list(iterations = 1000, adapt = FALSE, nrChains = 2)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)

getSample(out, start = 1, end = 5, which = 1)



summary(out)

## out is of class bayesianOutput. There are various standard functions 
# implemented for this output

plot(out)
correlationPlot(out)
marginalPlot(out)


## additionally, you can return the sample as a coda object, and make use of the coda functions
# for plotting and analysis

codaObject = getSample(out, start = 500, coda = TRUE)
coda::rejectionRate(codaObject)

