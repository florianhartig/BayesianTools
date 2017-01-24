#' Creates a Metropolis-type MCMC with options for covariance adaptatin, delayed rejection, Metropolis-within-Gibbs, and tempering
#' @author Florian Hartig
#' @param bayesianSetup either an object of class bayesianSetup created by \code{\link{createBayesianSetup}} (recommended), or a log target function 
#' @param settings a list of settings - possible options follow below
#' @param startValue startValue for the MCMC and optimization (if optimize = T). If not provided, the sampler will attempt to obtain the startValue from the bayesianSetup
#' @param optimize logical, determines whether an optimization for start values and proposal function should be run before starting the sampling
#' @param proposalGenerator option proposalgenerator object (see \code{\link{createProposalGenerator}})
#' @param proposalScaling additional scaling parameter for the proposals, needs to be as long as DRlevels. Defaults to 0.5^(- 0:(mcmcSampler$settings$DRlevels -1)
#' @param burnin number of iterations treated as burn-in. These iterations are not recorded in the chain.
#' @param thin thinning parameter. Determines the interval in which values are recorded.
#' @param consoleUpdates integer, determines the frequency with which sampler progress is printed to the console
#' @param adapt logical, determines wheter an adaptive algorithm should be implemented. Default is TRUE.
#' @param adaptationInterval integer, determines the interval  of the adaption if adapt = TRUE.
#' @param adaptationNotBefore integer, determines the start value for the adaption if adapt = TRUE.
#' @param DRlevels integer, determines the number of levels for a delayed rejection sampler. Default is 1, which means no delayed rejection is used.
#' @param temperingFunction function to implement simulated tempering in the algorithm. The function describes how the acceptance rate will be influenced in the course of the iterations.
#' @param gibbsProbabilities vector that defines the relative probabilities of the number of parameters to be changes simultaniously.
#' @param message logical determines whether the sampler's progress should be printed
#' @details The 'Metropolis' function is the main function for all Metropolis based samplers in this package. To call the derivatives from the basic Metropolis-Hastings MCMC, you can either use the corresponding function (e.g. \code{\link{AM}} for an adaptive Metropolis sampler) or use the parameters to adapt the basic Metropolis-Hastings. The advantage of the latter case is that you can easily combine different properties (e.g. adapive sampling and delayed rejection sampling) without changing the function.
#' @import coda
#' @example /inst/examples/MetropolisHelp.R
#' @export
Metropolis <- function(bayesianSetup, 
                       settings = list(startValue = NULL, 
                                       optimize = T, 
                                       proposalGenerator = NULL, 
                                       consoleUpdates=100, 
                                       burnin = 0,
                                       thin = 1, 
                                       parallel = NULL,
                                       adapt = T,
                                       adaptationInterval= 500,
                                       adaptationNotBefore = 3000,
                                       DRlevels = 1 ,
                                       proposalScaling = NULL,
                                       adaptationDepth = NULL,
                                       temperingFunction = NULL,
                                       gibbsProbabilities = NULL,
                                       message = TRUE
                                       )){
  
  ## General setup - this template should be similar for all MCMC algorithms
  
  setup <- checkBayesianSetup(bayesianSetup, parallel = settings$parallel) # calling parallel will check if requested parallelization in settings is provided by the BayesianSetup
  if(is.null(settings$parallel)) settings$parallel = bayesianSetup$parallel # checking back - if no parallelization is provided, we use the parallelization in the BayesianSetup. We could also set parallel = F, but I feel it makes more sense to use the Bayesiansetup as default
  
  settings = applySettingsDefault(settings, sampler = "Metropolis")
  
  if(is.null(settings$startValue)){
    settings$startValue = bayesianSetup$prior$sampler()
  }
  if(is.function(settings$startValue)){
    settings$startValue = settings$startValue()
  }
  
  ## Parameter consistency checks 

  if(is.null(settings$adaptationDepth)){
    settings$adaptationDepth = settings$adaptationNotBefore
  } 
  
  # Decreasing scaling for DRAM by default
  if (is.null(settings$proposalScaling)) settings$proposalScaling = 0.5^(- 0:(settings$DRlevels -1))
  
  tmp <- setupStartProposal(proposalGenerator = settings$proposalGenerator, bayesianSetup = bayesianSetup, settings = settings)
  settings = tmp$settings
  proposalGenerator = tmp$proposalGenerator

  ####### CREATE CHAIN
  
  chain = array(dim = c(1,bayesianSetup$numPars+3))
  chain[1,1:bayesianSetup$numPars] = settings$startValue
  colnames(chain) = c(1:bayesianSetup$numPars, "LP", "LL", "LPr")
  chain[1, (bayesianSetup$numPars+1):(bayesianSetup$numPars+3)] = setup$posterior$density(settings$startValue, returnAll = T)
  
  current = settings$startValue
  currentLP = as.numeric(chain[1, (bayesianSetup$numPars+1)])
  
  ##### Sampling
  
  classFields = list(
    setup = setup,
    settings = settings,
    current = current,
    currentLP = currentLP,
    chain = chain, 
    proposalGenerator = proposalGenerator,
    funEvals = 0,
    acceptanceRate = 0
  )
  
  class(classFields) <- c("mcmcSampler", "bayesianOutput")
  return(classFields)
}


#' gets samples while adopting the MCMC proposal generator
#' @author Florian Hartig
#' @param mcmcSampler an mcmcSampler 
#' @param iterations iterations
#' @description Function to sample with cobinations of the basic Metropolis-Hastings MCMC algorithm (Metropolis et al., 1953), a variation of the adaptive Metropolis MCMC (Haario et al., 2001), the delayed rejection algorithm (Tierney & Mira, 1999), and the delayed rejection adaptive Metropolis algorithm (DRAM, Haario et al), and the Metropolis within Gibbs 
#' @export
sampleMetropolis <- function(mcmcSampler, iterations){
  
  burnin <- mcmcSampler$settings$burnin
  thin <- mcmcSampler$settings$thin
  
  CounterFunEvals = mcmcSampler$funEvals
  CounterAccept = nrow(mcmcSampler$chain)*mcmcSampler$acceptanceRate
  
  if (mcmcSampler$settings$DRlevels > 2) stop("DRlevels > 2 currently not implemented")
  
  # Increase chain
  lastvalue = nrow(mcmcSampler$chain)
  mcmcSampler$chain = rbind(mcmcSampler$chain, array(dim = c(floor((iterations-burnin)/thin),mcmcSampler$setup$numPars+3)))
  
  alpha = rep(NA, mcmcSampler$settings$DRlevels)
  proposalEval = matrix( nrow = mcmcSampler$settings$DRlevels, ncol = 3)
  proposal = matrix( nrow = mcmcSampler$settings$DRlevels, ncol = mcmcSampler$setup$numPars)
  
  # Initialize counter for chain update
  counter <- lastvalue
  
  for (i in lastvalue:(lastvalue+iterations-1)){
  
    accepted = F
    
    if(is.null(mcmcSampler$settings$temperingFunction)) tempering = 1 else tempering = mcmcSampler$settings$temperingFunction(i)
    
    if(tempering < 1) warning("Tempering option < 1. This usually doesn't make sense!")
    
    for (j in 1:mcmcSampler$settings$DRlevels){
      
      proposal[j,] = mcmcSampler$proposalGenerator$returnProposal(x = mcmcSampler$current, scale = mcmcSampler$settings$proposalScaling[j])
      proposalEval[j,] <- mcmcSampler$setup$posterior$density(proposal[j,], returnAll = T)
      CounterFunEvals <- CounterFunEvals+1
      
      # case j = 1 (normal MH-MCMC)
      if (j == 1){
        alpha[j] = metropolisRatio(proposalEval[j,1], mcmcSampler$currentLP, tempering)
        jumpProbab = alpha[1]
        # case j = 2 (first delayed rejection)
      } else if (j == 2 & alpha[j-1] > 0 ){
        alpha[j] = metropolisRatio(proposalEval[j,1], proposalEval[j-1,1], tempering)
        
        temp <- metropolisRatio(mcmcSampler$proposalGenerator$returnDensity(proposal[1,], proposal[2,]), mcmcSampler$proposalGenerator$returnDensity(mcmcSampler$current, proposal[1,]))
        
        jumpProbab = metropolisRatio(proposalEval[j,1], mcmcSampler$currentLP, tempering) * temp * (1.0-alpha[j]) / (1.0-alpha[j-1]) 
      }
      
      if (runif(1) < jumpProbab){
        accepted = T
        mcmcSampler$current = proposal[j,]
        mcmcSampler$currentLP = proposalEval[j,1]
        if((i > (lastvalue+burnin)) && (i %% thin == 0) ){
          counter <- counter+1
        mcmcSampler$chain[counter,] = c(proposal[j,], proposalEval[j,])
        }
        break
      }
    }
    if((accepted == F) && (i > (lastvalue+burnin)) && (i %% thin == 0)){
      counter <- counter +1
      mcmcSampler$chain[counter,] = mcmcSampler$chain[counter-1,]
    } 
    if(accepted == T) CounterAccept <- CounterAccept+1
    
    # Proposal update
    
    if(mcmcSampler$settings$adapt  == T & i > mcmcSampler$settings$adaptationNotBefore & i %% mcmcSampler$settings$adaptationInterval == 0 ){
      start = max(1, counter - mcmcSampler$settings$adaptationDepth)
      mcmcSampler$proposalGenerator = updateProposalGenerator(proposal = mcmcSampler$proposalGenerator, chain = mcmcSampler$chain[start:counter,1:mcmcSampler$setup$numPars], message = F)
    }
    
    # Console update
    
    if(mcmcSampler$settings$message){
    if( i %% mcmcSampler$settings$consoleUpdates == 0 ) cat("\r","Running Metropolis-MCMC, chain ", 
                                                            mcmcSampler$settings$currentChain, "iteration" ,i,"of",iterations,
                                                            ". Current logp: ", mcmcSampler$chain[counter,mcmcSampler$setup$numPars+1]," Please wait!","\r")
    flush.console()
  }
  }
  
  # Make sure chain has right size
  mcmcSampler$chain <- mcmcSampler$chain[1:counter,]
  
  mcmcSampler$codaChain = coda::mcmc(mcmcSampler$chain)
  mcmcSampler$funEvals <- CounterFunEvals
  mcmcSampler$acceptanceRate <- CounterAccept/CounterFunEvals
  return(mcmcSampler)
}