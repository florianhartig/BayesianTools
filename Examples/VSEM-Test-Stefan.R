library(BayesianTools)
library(coda)


# Creating reference data
PAR <- VSEMcreatePAR(1:1000)
refPars   = VSEMgetDefaults()
refPars[12,] = c(0.2, 0.001, 1)
rownames(refPars)[12] = "error-sd"

referenceData <- VSEM(refPars$best[1:11], PAR) 
obs = apply(referenceData, 2, function(x) x + rnorm(length(x), 
                                                    sd = abs(x) * refPars$best[12]))

# Selecting parameters
parSel = c(1:6, 12)


PAR <- BayesianTools::VSEMcreatePAR(1:1000)

likelihood <- function(x, sum = T){
  refPars   = BayesianTools::VSEMgetDefaults()
  refPars[12,] = c(0.2, 0.001, 1)
  rownames(refPars)[12] = "error-sd"
  
  referenceData <- BayesianTools::VSEM(refPars$best[1:11], PAR) 
  obs = apply(referenceData, 2, function(x) x + rnorm(length(x), 
                                                      sd = abs(x) * refPars$best[12]))
  
  # Selecting parameters
  parSel = c(1:6, 12)
  mix = refPars$best
  mix[parSel] = x
  predicted <- BayesianTools::VSEM(mix[1:11], PAR)
  diff = c(predicted[,1:3] - obs[,1:3])
  llValues = dnorm(diff, sd = max(abs(c(predicted[,1:3])),0.0001) * mix[12], log = T) 
  if (sum == F) return(llValues)
  else return(sum(llValues))
}


refPars   = VSEMgetDefaults()
refPars[12,] = c(0.2, 0.001, 1)
rownames(refPars)[12] = "error-sd"

prior <- createUniformPrior(lower = refPars$lower[parSel], upper = refPars$upper[parSel])

# Bayesian Setup
BSVSEM <- createBayesianSetup(likelihood, prior, best = refPars$best[parSel], 
                              names = rownames(refPars)[parSel], parallel = 4)



startvalue = BSVSEM$prior$best

settings = list(iterations = 10000, adapt = F, optimize = F,
                startValue = startvalue, proposalScaling = 1)
test <- runMCMC(BSVSEM, sampler = "Metropolis", settings = settings)

plot(test)


settings = list(iterations = 10000)
test <- runMCMC(BSVSEM, sampler = "DEzs", settings = settings)

plot(test)



settings = list(iterations = 10000, startValue = startvalue, nrChains = 3)
test <- runMCMC(BSVSEM, sampler = "M", settings = settings)

plot(test)


proposalfunction  = test$proposalGenerator$returnProposal
posterior = BSVSEM$posterior$density

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,7))
  chain[1,] = startvalue
  current = posterior(startvalue)
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,], scale = 1)
  
    probab = posterior(proposal) - posterior(chain[i,])
    #print(probab)
    if (log(runif(1)) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}


chain = run_metropolis_MCMC(startvalue, 10000)
plot(mcmc(chain))
