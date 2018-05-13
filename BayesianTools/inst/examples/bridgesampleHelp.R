means <- c(0, 1, 2)
sds <- c(1, 0.6, 3)

# log-likelihood
ll <- function (x) {
  return(sum(dnorm(x, mean = means, sd = sds, log = TRUE)))
}

# lower and upper bounds for prior
lb <- rep(-10, 3)
ub <- rep(10, 3)

# create setup and run MCMC
setup <- createBayesianSetup(likelihood = ll,
                             lower = lb,
                             upper = ub)

out <- runMCMC(bayesianSetup = setup,
               settings = list(iterations = 1000),
               sampler = "DEzs")

# sample from MCMC output with "burn-in" of 25%
sample <- getSample(out$chain, start = 250, numSamples = 500)

# use bridge sampling to get marginal likelihood
bs_result <- bridgesample(chain = sample,
                          nParams = out$setup$numPars,
                          lower = lb,
                          upper = ub,
                          posterior = out$setup$posterior$density)
bs_result
