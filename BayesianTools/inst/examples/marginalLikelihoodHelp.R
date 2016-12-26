# Harmonic mean works OK for a low-dim case with 

likelihood <- function(x) sum(dnorm(x, log = TRUE))
prior = createUniformPrior(lower = rep(-1,2), upper = rep(1,2))
bayesianSetup <- createBayesianSetup(likelihood = likelihood, prior = prior)
out = runMCMC(bayesianSetup = bayesianSetup, settings = list(iterations = 5000))

plot(out)


marginalLikelihood(out, numSamples = 500)[[1]]
marginalLikelihood(out, method = "HM", numSamples = 500)[[1]]
marginalLikelihood(out, method = "Prior", numSamples =  500)[[1]]

# True marginal likelihood (brute force approximation)

marginalLikelihood(out, method = "Prior", numSamples =  10000)[[1]]


# Harmonic mean goes totally wrong for higher dimendsions - wide prior. 
# Could also be a problem of numeric stability of the implementation

likelihood <- function(x) sum(dnorm(x, log = TRUE))
prior = createUniformPrior(lower = rep(-10,3), upper = rep(10,3))
bayesianSetup <- createBayesianSetup(likelihood = likelihood, prior = prior)
out = runMCMC(bayesianSetup = bayesianSetup, settings = list(iterations = 5000))

plot(out)

marginalLikelihood(out, numSamples = 500)[[1]]
marginalLikelihood(out, method = "HM", numSamples = 500)[[1]]
marginalLikelihood(out, method = "Prior", numSamples =  500)[[1]]

# True marginal likelihood (brute force approximation)

marginalLikelihood(out, method = "Prior", numSamples =  10000)[[1]]


