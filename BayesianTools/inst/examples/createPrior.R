# The BT package includes a number of convenience functions to specify
# prior distributions, including createUniformPrior, createTruncatedNormalPrior
# etc. If you want to specify a prior that corresponds to one of these
# distributions, you should use these functions, e.g.:
#
# The prior we choose depends on the prior information we have. For example, if 
# we have no prior information, we can choose a uniform prior. The normal 
# distribution is often used to model a wide range of phenomena in statistics, 
# such as the distribution of heights or weights in a population. Beta 
# distribution, on the other hand, is defined on the interval [0,1]. 
# It is often used to model random variables that represent proportions, 
# probabilities or other values that are constrained to lie within this interval.

# | createPrior | createBetaPrior | createUniformPrior | createTruncatedNormalPrior |
# | :----------:|:---------------:|:------------------:|:--------------------------:|
# |Any distribution provided by the user|Beta distribution|Uniform distribution|Normal distribution|

  
  

prior <- createUniformPrior(lower = c(0,0), upper = c(0.4,5))

prior$density(c(2, 3)) # outside of limits -> -Inf
prior$density(c(0.2, 2)) # within limits, -0.6931472

# All default priors include a sampling function, i.e. you can create
# samples from the prior via
prior$sampler()
# [1] 0.2291413 4.5410389

# if you want to specify a prior that does not have a default function, 
# you should use the createPrior function, which expects a density and 
# optionally a sampler function:

density = function(par){
  d1 = dunif(par[1], -2,6, log =TRUE)
  d2 = dnorm(par[2], mean= 2, sd = 3, log =TRUE)
  return(d1 + d2)
}

sampler = function(n=1){
  d1 = runif(n, -2,6)
  d2 = rnorm(n, mean= 2, sd = 3)
  return(cbind(d1,d2))
}

prior <- createPrior(density = density, sampler = sampler, 
                     lower = c(-10,-20), upper = c(10,20), best = NULL)

# note that the createPrior supports additional truncation


# To use a prior in an MCMC, include it in a BayesianSetup 

set.seed(123)
ll <- function(x) sum(dnorm(x, log = TRUE)) # multivariate normal ll
bayesianSetup <- createBayesianSetup(likelihood = ll, prior = prior)

settings = list(iterations = 100)
out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)

# use createPriorDensity to create a new (estimated) prior from MCMC output

newPrior = createPriorDensity(out, method = "multivariate",
                              eps = 1e-10, lower = c(-10,-20), 
                              upper = c(10,20), best = NULL, scaling = 0.5)

