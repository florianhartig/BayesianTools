## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(fig.width=5, fig.height=5, warning=FALSE, cache = F)

## ----echo = F, message = F----------------------------------------------------
set.seed(123)

## ----eval = F-----------------------------------------------------------------
#  install.packages("BayesianTools")

## -----------------------------------------------------------------------------
library(BayesianTools)
citation("BayesianTools")

## -----------------------------------------------------------------------------
set.seed(123)

## ----eval = F-----------------------------------------------------------------
#  sessionInfo()

## -----------------------------------------------------------------------------
ll <- generateTestDensityMultiNormal(sigma = "no correlation")
bayesianSetup = createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))

## -----------------------------------------------------------------------------
iter = 10000
settings = list(iterations = iter, message = FALSE)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)

## -----------------------------------------------------------------------------
print(out)
summary(out)

## -----------------------------------------------------------------------------
plot(out) # plot internally calls tracePlot(out)
correlationPlot(out)
marginalPlot(out, prior = TRUE)

## -----------------------------------------------------------------------------
marginalLikelihood(out)
DIC(out)
MAP(out)

## ----eval = F-----------------------------------------------------------------
#  getSample(out, start = 100, end = NULL, thin = 5, whichParameters = 1:2)

## ----echo = T-----------------------------------------------------------------
iter = 1000
settings = list(iterations = iter, nrChains = 3, message = FALSE)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)


## -----------------------------------------------------------------------------
print(out)
summary(out)

## -----------------------------------------------------------------------------
plot(out)

## -----------------------------------------------------------------------------
#getSample(out, coda = F)
gelmanDiagnostics(out, plot = T)

## ----eval = F-----------------------------------------------------------------
#  ll = logDensity(x)

## -----------------------------------------------------------------------------
ll <- generateTestDensityMultiNormal(sigma = "no correlation")
bayesianSetup = createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))

## ----eval = FALSE-------------------------------------------------------------
#  ## Definition of likelihood function
#  likelihood <- function(matrix){
#  	# Calculate likelihood in parallel
#  	# Return vector of likelihood valus
#  }
#  
#  ## Create Bayesian Setup
#  BS <- createBayesianSetup(likelihood, parallel = "external", ...)
#  
#  ## Run MCMC
#  runMCMC(BS, sampler = "SMC", ...)

## ----eval = FALSE-------------------------------------------------------------
#  ## n = Number of cores
#  n=2
#  x <- c(1:10)
#  likelihood <- function(param) return(sum(dnorm(x, mean = param, log = T)))
#  bayesianSetup <- createBayesianSetup(likelihood, parallel = n, lower = -5, upper = 5)
#  
#  ## give runMCMC a matrix with n rows of proposals as startValues or sample n times from the previous created sampler
#  out <- runMCMC(bayesianSetup, settings = list(iterations = 1000))

## ----eval = FALSE-------------------------------------------------------------
#  ### Create cluster with n cores
#  cl <- parallel::makeCluster(n)
#  
#  ## Definition of the likelihood
#  likelihood  <- function(X) sum(dnorm(c(1:10), mean = X, log = T))
#  
#  ## Definition of the likelihood which will be calculated in parallel. Instead of the parApply function, we could also define a costly parallelized likelihood
#  pLikelihood <- function(param) parallel::parApply(cl = cl, X = param, MARGIN = 1, FUN = likelihood)
#  
#  ## export functions, dlls, libraries
#  # parallel::clusterEvalQ(cl, library(BayesianTools))
#  parallel::clusterExport(cl, varlist = list(likelihood))
#  
#  ## create BayesianSetup
#  bayesianSetup <- createBayesianSetup(pLikelihood, lower = -10, upper = 10, parallel = 'external')
#  
#  ## For this case we want to parallelize the internal chains, therefore we create a n row matrix with startValues, if you parallelize a model in the likelihood, do not set a n*row Matrix for startValue
#  settings = list(iterations = 100, nrChains = 1, startValue = bayesianSetup$prior$sampler(n))
#  
#  ## runMCMC
#  out <- runMCMC(bayesianSetup, settings, sampler = "DEzs")

## ----eval = FALSE-------------------------------------------------------------
#  ### Create cluster with n cores
#  cl <- parallel::makeCluster(n)
#  
#  ## export your model
#  # parallel::clusterExport(cl, varlist = list(complexModel))
#  
#  ## Definition of the likelihood
#  likelihood  <- function(param) {
#    # ll <- complexModel(param)
#    # return(ll)
#  }
#  
#  ## create BayesianSetup and settings
#  bayesianSetup <- createBayesianSetup(likelihood, lower = -10, upper = 10, parallel = 'external')
#  settings = list(iterations = 100, nrChains = 1)
#  
#  ## runMCMC
#  out <- runMCMC(bayesianSetup, settings)
#  

## ----eval = FALSE-------------------------------------------------------------
#  ### Definition of likelihood function
#  x <- c(1:10)
#  likelihood <- function(param) return(sum(dnorm(x, mean = param, log = T)))
#  
#  ## Create BayesianSetup and settings
#  bayesianSetup <- createBayesianSetup(likelihood, lower = -10, upper = 10, parallel = F)
#  settings = list(iterations = 100000)
#  
#  ## Start cluster with n cores for n chains and export BayesianTools library
#  cl <- parallel::makeCluster(n)
#  parallel::clusterEvalQ(cl, library(BayesianTools))
#  
#  ## calculate parallel n chains, for each chain the likelihood will be calculated on one core
#  out <- parallel::parLapply(cl, 1:n, fun = function(X, bayesianSetup, settings) runMCMC(bayesianSetup, settings, sampler = "DEzs"), bayesianSetup, settings)
#  
#  ## Combine the chains
#  out <- createMcmcSamplerList(out)

## -----------------------------------------------------------------------------
# Create a BayesianSetup
ll <- generateTestDensityMultiNormal(sigma = "no correlation")
bayesianSetup = createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))

settings = list(iterations = 2500,  message = FALSE)
out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)


newPrior = createPriorDensity(out, method = "multivariate", eps = 1e-10, lower = rep(-10, 3), upper =  rep(10, 3), best = NULL)
bayesianSetup <- createBayesianSetup(likelihood = ll, prior = newPrior)

settings = list(iterations = 1000,  message = FALSE)
out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)

## ----message = F--------------------------------------------------------------

ll <- generateTestDensityMultiNormal(sigma = "no correlation")

bayesianSetup = createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))

settings = list(iterations = 10000, nrChains= 3, message = FALSE)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)

plot(out)
marginalPlot(out, prior = T)
correlationPlot(out)
gelmanDiagnostics(out, plot=T)

# option to restart the sampler

settings = list(iterations = 1000, nrChains= 1, message = FALSE)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)

out2 <- runMCMC(bayesianSetup = out)

out3 <- runMCMC(bayesianSetup = out2)

#plot(out)
#plot(out3)

# create new prior from posterior sample 

newPriorFromPosterior <- createPriorDensity(out2)


## -----------------------------------------------------------------------------
iter = 10000

## -----------------------------------------------------------------------------
applySettingsDefault(sampler = "Metropolis")

## ----results = 'hide', eval = F-----------------------------------------------
#  settings <- list(iterations = iter, adapt = F, DRlevels = 1, gibbsProbabilities = NULL, temperingFunction = NULL, optimize = F, message = FALSE)
#  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
#  plot(out)

## ----results = 'hide', eval = F-----------------------------------------------
#  settings <- list(iterations = iter, adapt = F, DRlevels = 1, gibbsProbabilities = NULL, temperingFunction = NULL, optimize = T, message = FALSE)
#  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
#  plot(out)

## ----results = 'hide', eval = F-----------------------------------------------
#  settings <- list(iterations = iter, adapt = T, DRlevels = 1, gibbsProbabilities = NULL, temperingFunction = NULL, optimize = T,  message = FALSE)
#  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
#  plot(out)

## ----results = 'hide', eval = F-----------------------------------------------
#  settings <- list(iterations = iter, adapt = F, DRlevels = 2, gibbsProbabilities = NULL, temperingFunction = NULL, optimize = T,  message = FALSE)
#  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
#  plot(out)

## ----results = 'hide', eval = F-----------------------------------------------
#  settings <- list(iterations = iter, adapt = T, DRlevels = 2, gibbsProbabilities = NULL, temperingFunction = NULL, optimize = T,  message = FALSE)
#  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
#  plot(out)

## ----results = 'hide', eval = F-----------------------------------------------
#  settings <- list(iterations = iter, adapt = T, DRlevels = 1, gibbsProbabilities = c(1,0.5,0), temperingFunction = NULL, optimize = T,  message = FALSE)
#  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
#  plot(out)

## ----results = 'hide', eval = F-----------------------------------------------
#  temperingFunction <- function(x) 5 * exp(-0.01*x) + 1
#  settings <- list(iterations = iter, adapt = F, DRlevels = 1, gibbsProbabilities = c(1,1,0), temperingFunction = temperingFunction, optimize = T,  message = FALSE)
#  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
#  plot(out)

## ----results = 'hide', eval = F-----------------------------------------------
#  settings <- list(iterations = iter,  message = FALSE)
#  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DE", settings = settings)
#  plot(out)

## ----results = 'hide', eval = F-----------------------------------------------
#  settings <- list(iterations = iter,  message = FALSE)
#  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
#  plot(out)

## ----results = 'hide', eval = F-----------------------------------------------
#  settings <- list(iterations = iter,  message = FALSE)
#  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DREAM", settings = settings)
#  plot(out)

## ----results = 'hide', eval = FALSE-------------------------------------------
#  settings <- list(iterations = iter,  message = FALSE)
#  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DREAMzs", settings = settings)
#  plot(out)

## ----eval = F-----------------------------------------------------------------
#  settings = list(iterations = iter,  message = FALSE)
#  
#  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Twalk", settings = settings)

## ----eval = T-----------------------------------------------------------------
settings <- list(iterations = iter, nrChains = 3,  message = FALSE)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
plot(out)

#chain = getSample(out, coda = T)
gelmanDiagnostics(out, plot = F)

## ----results = 'hide', eval = F-----------------------------------------------
#  settings <- list(initialParticles = iter, iterations= 1)
#  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "SMC", settings = settings)
#  plot(out)

## ----results = 'hide', eval = F-----------------------------------------------
#  settings <- list(initialParticles = iter, iterations= 10)
#  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "SMC", settings = settings)
#  plot(out)

## -----------------------------------------------------------------------------
sampleSize = 30
x <- (-(sampleSize-1)/2):((sampleSize-1)/2)
y <-  1 * x + 1*x^2 + rnorm(n=sampleSize,mean=0,sd=10)
plot(x,y, main="Test Data")

## -----------------------------------------------------------------------------
likelihood1 <- function(param){
    pred = param[1] + param[2]*x + param[3] * x^2
    singlelikelihoods = dnorm(y, mean = pred, sd = 1/(param[4]^2), log = T)
    return(sum(singlelikelihoods))  
}

likelihood2 <- function(param){
    pred = param[1] + param[2]*x 
    singlelikelihoods = dnorm(y, mean = pred, sd = 1/(param[3]^2), log = T)
    return(sum(singlelikelihoods))  
}

## -----------------------------------------------------------------------------
setUp1 <- createBayesianSetup(likelihood1, lower = c(-5,-5,-5,0.01), upper = c(5,5,5,30))

setUp2 <- createBayesianSetup(likelihood2, lower = c(-5,-5,0.01), upper = c(5,5,30))

## ----results = 'hide'---------------------------------------------------------
settings = list(iterations = 15000,  message = FALSE)
out1 <- runMCMC(bayesianSetup = setUp1, sampler = "Metropolis", settings = settings)
#tracePlot(out1, start = 5000)
M1 = marginalLikelihood(out1)
M1

settings = list(iterations = 15000,  message = FALSE)
out2 <- runMCMC(bayesianSetup = setUp2, sampler = "Metropolis", settings = settings)
#tracePlot(out2, start = 5000)
M2 = marginalLikelihood(out2)
M2

## -----------------------------------------------------------------------------
exp(M1$ln.ML - M2$ln.ML)

## -----------------------------------------------------------------------------
exp(M1$ln.ML) / ( exp(M1$ln.ML) + exp(M2$ln.ML))

## -----------------------------------------------------------------------------
DIC(out1)$DIC
DIC(out2)$DIC

## -----------------------------------------------------------------------------
# This will not work, since likelihood1 has no sum argument
# WAIC(out1)

# likelihood with sum argument
likelihood3 <- function(param, sum = TRUE){
    pred <- param[1] + param[2]*x + param[3] * x^2
    singlelikelihoods <- dnorm(y, mean = pred, sd = 1/(param[4]^2), log = T)
    return(if (sum == TRUE) sum(singlelikelihoods) else singlelikelihoods)  
}
setUp3 <- createBayesianSetup(likelihood3, lower = c(-5,-5,-5,0.01), upper = c(5,5,5,30))
out3 <- runMCMC(bayesianSetup = setUp3, sampler = "Metropolis", settings = settings)

WAIC(out3)

