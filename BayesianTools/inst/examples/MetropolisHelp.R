
# Running the metropolis via the runMCMC with a proposal covariance generated from the prior 
# (can be useful for complicated priors)

ll = function(x) sum(dnorm(x, log = TRUE))
setup = createBayesianSetup(ll, lower = c(-10,-10), upper = c(10,10))

samples = setup$prior$sampler(1000)

generator = createProposalGenerator(diag(1, setup$numPars))
generator = updateProposalGenerator(generator, samples, manualScaleAdjustment = 1, message = TRUE)

settings =  list(proposalGenerator = generator, optimize = FALSE, iterations = 500)  

out = runMCMC(bayesianSetup = setup, sampler = "Metropolis", settings = settings)
