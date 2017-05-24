ll = function(x) sum(dnorm(x, log = T))

setup = createBayesianSetup(ll, lower = c(-10,-10), upper = c(10,10))

settings = list(nrChains = 2, iterations = 1000)
out <- runMCMC(bayesianSetup = setup, sampler = "DEzs", settings = settings)


# population MCMCs divide the interations by the number of internal chains, so the end of the 3 chains is 1000/3 = 334
sample <- getSample(out, start = 100, end = 334, thin = 10) 

# 
sample <- getSample(out, start = 100, numSamples = 60, coda = T)
plot(sample)
