# This example show how to post-hoc impose a prior on a MCMC sample from a likelihood with a flat prior

library(BayesianTools)

# likelihood sampled with flat prior 

ll = function(x) sum(dnorm(x, log = TRUE))
setup = createBayesianSetup(ll, lower = c(-100,-100), upper = c(100,100))
out <- runMCMC(bayesianSetup = setup, sampler = "DEzs")
summary(out)


# post-hoc impose prior by rejection / sampling

prior = function(x) sum(dnorm(x, log = TRUE, mean = 2))
setup2 = createBayesianSetup(prior, lower = c(-100,-100), upper = c(100,100))

out2 = smcSampler(setup2, initialParticles = getSample(out, start = 1000, numSamples = 2000), iterations = 30)

plot(out2, plotPrior = F)


# Alternative with MCMC and createPrior, summarizes old MCMC via multivariateNormal

oldPosterior = createPriorDensity(getSample(out, start = 1000, numSamples = 2000))
prior = function(x) sum(dnorm(x, log = TRUE, mean = 5))


setup3 = createBayesianSetup(likelihood = prior, prior = oldPosterior)

out3 = runMCMC(bayesianSetup = setup3)
summary(out3)


# should be identical to 

ll = function(x) sum(dnorm(x, log = TRUE))
prior = createTruncatedNormalPrior(mean = c(2,2), sd = c(1,1), lower = c(-100,-100), upper = c(100,100) )
setup = createBayesianSetup(ll, prior)
out <- runMCMC(bayesianSetup = setup, sampler = "DEzs")
summary(out)





