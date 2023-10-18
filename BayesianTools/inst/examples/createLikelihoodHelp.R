

# likelihood function can be created with or without createBayesianSetup. ??createBayesianSetup provides additional properties.

data = rnorm(20)

# create standardized likelihood density for Gaussian likelihood function
# fitting parameters mean and standard deviation

likelihood <- function(x) {
  predicted <- rep(x[1], length(data))
  LL = likelihoodIidNormal(predicted, data, x[2])
  return(LL)
}

likelihood(x = c(1,1))

# The next step is optional, you can provide the likelihood density directly
# to BayesianSetup and createLikelihood will be called automatically.
# The purpose of the createLikelihood function is to parallelize the 
# density function and to add some convenience functions, 
# for example to catch exceptions that are raised by the call to likelihood
likelihoodClass <- createLikelihood(likelihood)

bayesianSetup <- createBayesianSetup(likelihood = likelihoodClass, lower = c(-10,0.001), upper = c(10, 5))

settings = list(iterations = 1000)
out <- runMCMC(bayesianSetup, sampler = "DEzs", settings)
summary(out)



# simulate temporally autocorrelated data
AR1sim<-function(n, a){
  x = rep(NA, n)
  x[1] = 0
  for(i in 2:n){
    x[i] = a * x[i-1] + (1-a) * rnorm(1)
  }
  return(x)
}

set.seed(123)
data = AR1sim(1000, 0.5)
plot(data, type = "l") 


# create likelihood function for AR1 likelihood function
likelihood <- function(x) {
  predicted <- rep(x[1], length(data))
  LL = likelihoodAR1(predicted, data, x[2], x[3])
  return(LL)
}

bayesianSetup <- createBayesianSetup(likelihood = likelihood, 
                                lower = c(-10,0.001, 0), upper = c(10, 5, 0.99))

settings = list(iterations = 1000)
out <- runMCMC(bayesianSetup, sampler = "DEzs", settings)
summary(out, start = 200)
plot(out, start = 200)

