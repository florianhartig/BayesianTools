
# Comparison of ML for two regression models

sampleSize = 30
x <- (-(sampleSize-1)/2):((sampleSize-1)/2)
y <-  1 * x + 1*x^2 + rnorm(n=sampleSize,mean=0,sd=10)
#plot(x,y, main="Test Data")


# linear and quadratic effect
likelihood1 <- function(param){
  pred = param[1] + param[2]*x + param[3] * x^2
  singlelikelihoods = dnorm(y, mean = pred, sd = 1/(param[4]^2), log = T)
  return(sum(singlelikelihoods))  
}

# linear  effect
likelihood2 <- function(param){
  pred = param[1] + param[2]*x 
  singlelikelihoods = dnorm(y, mean = pred, sd = 1/(param[3]^2), log = T)
  return(sum(singlelikelihoods))  
}

setUp1 <- createBayesianSetup(likelihood1, lower = c(-5,-5,-5,0.01), upper = c(5,5,5,30))

setUp2 <- createBayesianSetup(likelihood2, lower = c(-5,-5,0.01), upper = c(5,5,30))

out1 <- runMCMC(bayesianSetup = setUp1)
M1 = marginalLikelihood(out1, start = 1000)

out2 <- runMCMC(bayesianSetup = setUp2)
M2 = marginalLikelihood(out2, start = 1000)


### Bayes factor

exp(M1$ln.ML - M2$ln.ML)

# BF > 1 means the evidence is in favor of M1. See Kass, R. E. & Raftery, A. E. 
# (1995) Bayes Factors. J. Am. Stat. Assoc., Amer Statist Assn, 90, 773-795.

### Posterior weight

exp(M1$ln.ML) / ( exp(M1$ln.ML) + exp(M2$ln.ML))

# If models have different model priors, multiply with the prior probabilities of each model. 

############################################################

### Performance comparison ### 


# Low dimensional case with narrow priors - all methods have low error

# we use a truncated normal for the likelihood to make sure that the density 
# integrates to 1 - makes it easier to calcuate the theoretical ML
likelihood <- function(x) sum(msm::dtnorm(x, log = TRUE, lower = -1, upper = 1))
prior = createUniformPrior(lower = rep(-1,2), upper = rep(1,2))
bayesianSetup <- createBayesianSetup(likelihood = likelihood, prior = prior)
out = runMCMC(bayesianSetup = bayesianSetup, settings = list(iterations = 5000))

# plot(out)

# theoretical value
theory = log(1/(2^2))

marginalLikelihood(out)$ln.ML - theory
marginalLikelihood(out, method = "Prior", numSamples =  500)$ln.ML - theory
marginalLikelihood(out, method = "HM", numSamples =  500)$ln.ML - theory
marginalLikelihood(out, method = "Bridge", numSamples =  500)$ln.ML - theory


# higher dimensions - wide prior - HM and bridge don't work.

likelihood <- function(x) sum(msm::dtnorm(x, log = TRUE, lower = -10, upper = 10))
prior = createUniformPrior(lower = rep(-10,3), upper = rep(10,3))
bayesianSetup <- createBayesianSetup(likelihood = likelihood, prior = prior)
out = runMCMC(bayesianSetup = bayesianSetup, settings = list(iterations = 5000))

# plot(out)

# theoretical value
theory = log(1/(20^3))

marginalLikelihood(out)$ln.ML - theory
marginalLikelihood(out, method = "Prior", numSamples =  500)$ln.ML - theory
marginalLikelihood(out, method = "HM", numSamples =  500)$ln.ML - theory
marginalLikelihood(out, method = "Bridge", numSamples =  500)$ln.ML - theory



