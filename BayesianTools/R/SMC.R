#' SMC sampler
#' @author Florian Hartig
#' @description Sequential Monte Carlo Sampler
#' @param bayesianSetup either an object of class bayesianSetup created by \code{\link{createBayesianSetup}} (recommended), or a log target function 
#' @param initialParticles initial particles - either a draw from the prior, provided as a matrix with the single parameters as columns and each row being one particle (parameter vector), or a numeric value with the number of desired particles. In this case, the sampling option must be provided in the prior of the BayesianSetup. 
#' @param iterations number of iterations
#' @param resampling logical, specifies whether new particles should be created at each iteration
#' @param resamplingSteps how many resampling (MCMC) steps between the iterations
#' @param proposal optional, proposal class
#' @param adaptive logical, should the covariance of the proposal be adapted during sampling?
#' @param proposalScale scaling factor for the proposal generation. Can be adapted if there is too much / too little rejection
#' @details The sampler can be used for rejection sampling as well as for sequential Monte Carlo. For the former case set the iterations to one.
#' 
#' @note The SMC currently assumes that the initial particle is sampled from the prior. If a better initial estimate of the posterior distribution is available, this the sampler should be modified to include this. Currently, however, this is not included in the code, so the appropriate adjustments have to be done by hand. 
#' @export
#' @example /inst/examples/SMCHelp.R
smcSampler <- function(bayesianSetup, initialParticles = 1000, iterations = 10, resampling = T, resamplingSteps = 2, proposal = NULL, adaptive = T, proposalScale = 0.5){
  
  if(resamplingSteps < 1) stop("SMC error, resamplingSteps can't be < 1")
    
  setup <- checkBayesianSetup(bayesianSetup)
  
  info = list()
  info$resamplingAcceptance = matrix(nrow = iterations, ncol = resamplingSteps)
  info$survivingParticles = rep(NA, iterations)
  
  
  if(inherits(initialParticles, "numeric")){
    initialParticles = bayesianSetup$prior$sampler(initialParticles)
  }
  
  if (any(is.infinite(setup$prior$density(initialParticles)))) stop("initialParticles outside prior range")
  
  particles <- initialParticles
  
  rejectionRate = 0
  
  particleSize = nrow(initialParticles)
  
  acceptanceTarget = round(particleSize / 2)
  
  posterior = matrix(nrow = particleSize, ncol = 3)
    
  numPar <- ncol(initialParticles)
  
  if (is.null(proposal)) proposalGenerator = createProposalGenerator(rep(40,numPar))
  
  usedUp = 0
  
  for (i in 1:iterations){
    
    posterior = setup$posterior$density(particles, returnAll = T)
    
    likelihoodValues <- posterior[,2]
    
    # idea - adjust (1/iterations) such that always approx 30% of particles are maintain
    #level = sort(likelihoodValues)[acceptanceTarget]
    #best = likelihoodValues
    
    ll = likelihoodValues - max(likelihoodValues, na.rm = T)
    
    llCutoff = sort(ll)[acceptanceTarget] 
    
    
    
    relativeL = exp(likelihoodValues - max(likelihoodValues, na.rm = T))^(1/iterations)
  
    sel = sample.int(n=length(likelihoodValues), size = length(likelihoodValues), replace = T, prob = relativeL)
    
    info$survivingParticles[i] = length(unique(sel))
    
    particles = particles[sel,]
    
    if (numPar == 1) particles = matrix(particles, ncol = 1)
    
    if (resampling == T){
      
      if (adaptive == T){
        proposalGenerator = updateProposalGenerator(proposalGenerator, particles)
      }
      
      for (j in 1:resamplingSteps){
        particlesProposals = proposalGenerator$returnProposalMatrix(particles, scale = proposalScale)
        
        jumpProb <- exp(setup$posterior$density(particlesProposals) - likelihoodValues[sel])^(i/iterations) * exp(setup$prior$density(particlesProposals)   - setup$prior$density(particles))
        
        accepted <- jumpProb > runif(length(jumpProb), 0 ,1)
        
        rejectionRate = rejectionRate + sum(accepted)
        
        particles[accepted, ] = particlesProposals[accepted, ]        
        
      
      }
    }
  }
  
  info$rejectionRate = rejectionRate / (iterations * resamplingSteps)

  out = list(
    setup = setup,
    initialParticles = initialParticles,
    particles = particles,
    posteriorValues = posterior,
    proposalGenerator = proposalGenerator,
    info = info
  )
  
  class(out) <- c("smcSampler", "bayesianOutput")
  return(out)
  
}
