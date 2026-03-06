# data comes from normal distribution 

data = rnorm(100, mean = 1.2, sd = 0.7)

# define a density function with unknown mean and sd
# needs to return LOG! likelihood

likelihood <- function(x) {
  LL = sum(dnorm(data, mean = x[1], sd = x[2], log = T))
  return(LL)
}

likelihood(x = c(1,1))

# To be used in BayesianTools, a standardized likelihood class should be 
# created via createLikelihood. 

likelihoodClass <- createLikelihood(likelihood)

# This class parallelizes the simple density function and to add
# some convenience functions, for example 
# to catch exceptions that are raised by the call to likelihood
# for example, we can now provide a matrix which will result in
# several likelihood values

likelihoodClass$density(matrix(c(1,2,1,1), nrow = 2))

# createLikelihood is automatically called when you provide the likelihood
# to a BayesianSetup, so it is usually sufficient to just provide 
# the likelihood to createBayesianSetup

bayesianSetup <- createBayesianSetup(likelihood = likelihoodClass, 
                                     lower = c(-2,0.001), 
                                     upper = c(5, 5))

settings = list(iterations = 500)
out <- runMCMC(bayesianSetup, sampler = "DEzs", settings)
summary(out, start = 50)
plot(out) 


