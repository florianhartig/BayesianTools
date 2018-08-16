#' SMC sampler
#' @author Florian Hartig
#' @description Sequential Monte Carlo Sampler
#' @param bayesianSetup either an object of class bayesianSetup created by \code{\link{createBayesianSetup}} (recommended), or a log target function 
#' @param initialParticles initial particles - either a draw from the prior, provided as a matrix with the single parameters as columns and each row being one particle (parameter vector), or a numeric value with the number of desired particles. In this case, the sampling option must be provided in the prior of the BayesianSetup. 
#' @param iterations number of iterations
#' @param resampling if new particles should be created at each iteration
#' @param resamplingSteps how many resampling (MCMC) steps between the iterations
#' @param proposal optional proposal class
#' @param exponents series of exponents to build the intermediate distributions
#' @param adaptive should the covariance of the proposal be adapted during sampling
#' @param proposalScale scaling factor for the proposal generation. Can be adapted if there is too much / too little rejection
#' @param x Parameter to generate the exponential sequence for building intermediary distributions. Default value from Jeremiah et al. (2012)
#' @param m Parameter to generate the exponential sequence for building intermediary distributions. Default value from Jeremiah et al. (2012)
#' @param sampling Which algorithm to use for particle (re)sampling. Options are "multinomial" (default), "residual" and "systematic"
#' @param ess.limit Threshold value of effective sample size below which resampling is done. By default, the value is set to half the number of particles. To resample at each step, use a value >= the number of particles.
#' @param lastResample Iteration (starting from the end) at which particle resampling is forced. To deactivate this, set to a value < 0
#' @param pars.lower Optional: vector contaiing the mimimum values of calibration parameters. Used to initialize the proposal function
#' @param pars.upper Optional: vector contaiing the maximum values of calibration parameters. Used to initialize the proposal function
#' @param diagnostics an optional function with diagnostics that are calculated on the particles during each SMC iteration

#' @details The sampler can be used for rejection sampling as well as for sequential Monte Carlo. For the former case set the iterations to one.
#' @note The SMC currently assumes that the initial particle is sampled from the prior. If a better initial estimate of the posterior distribution is available, this the sampler should be modified to include this. Currently, however, this is not included in the code, so the appropriate adjustments have to be done by hand. 
#' @export
#' @example /inst/examples/SMCHelp.R
smcSampler <- function(bayesianSetup, 
                       initialParticles = 1000,
                       iterations = 10, 
                       resampling = T, 
                       resamplingSteps = 2, 
                       lastMutateSteps = 5,         # TODO document
                       proposal = NULL, 
                       exponents = NULL, 
                       adaptive = T, 
                       proposalScale = 0.5, 
                       x=3.11, 
                       m=7E-08, 
                       sampling="multinomial",
                       ess.limit=NULL,
                       lastResample = 1,
                       pars.lower=NULL, 
                       pars.upper=NULL, 
                       mutate.method ="Metropolis", # TODO document
                       b=1e-04,                     # TODO document
                       diagnostics = NULL){
  
  ########### SETUP STEPS ########################
  
  if(resamplingSteps < 1) stop("SMC error, resamplingSteps can't be < 1")
  
  info = list()
  info$resamplingAcceptance = as.data.frame(matrix(nrow = iterations, ncol = resamplingSteps)) #Using DF instead of matrix because number of iterations may change due to adaptive algorithm
  info$survivingParticles = rep(NA, iterations)
  info$ess.vec <- vector("numeric", iterations)
  info$dist.vec <- vector("numeric", iterations)
  info$res.of <- rep(FALSE, iterations)
  info$res.ess <- rep(FALSE, iterations)

  setup <- checkBayesianSetup(bayesianSetup)
  
  ### InitialParticles
  #TODO documentation
  
  if(class(initialParticles) == "numeric"){
    initialParticles = bayesianSetup$prior$sampler(initialParticles)
    importanceDensity = bayesianSetup$prior$density
  }
  if(class(initialParticles) == "matrix"){
    importanceDensity = bayesianSetup$prior$density
  }
  if(class(initialParticles) == "list"){
    importanceDensity = initialParticles$density
    initialParticles = as.matrix(initialParticles$particles,ncol=1)
  }
  
  #print(c("density", setup$prior$density(initialParticles)))
  
  if (any(is.infinite(setup$prior$density(initialParticles)))) stop("initialParticles outside prior range")
  
  particles <- initialParticles
  rejectionRate = 0
  particleSize = nrow(initialParticles)
  if(is.null(ess.limit)){ess.limit <- round(particleSize) * 0.5}
  
  weights <- oldweights <- oldInter <- rep(0, particleSize)
  
  posterior = matrix(nrow = particleSize, ncol = 3)
  numPar <- ncol(initialParticles)
  
  if(mutate.method %in% c("Metropolis", "adaptive")){
    if (is.null(proposal)) proposalGenerator = createProposalGenerator(rep(40,numPar))
  }
  
  usedUp = 0

  
  ########### WEIGHT SETUP ########################
    
  # Define an exponential sequence of beta parameters, i.e. the exponents to weight the likelihood against prior.
  # Idea and parameter values from Jeremiah et al. (2012), Environ Modell Softw
  if(is.null(exponents)){
    iters <- seq(0,iterations)
    exponents <- m * (iters*200/iterations)^x
    exponents <- pmin(1,exponents)
    exponents <- exponents[2:length(exponents)]
  }
  print(c("length xponents", length(exponents)))
  print(c("exponents", exponents))
  
  # Initial value for exponent: here the first value of the series is taken.
  # If the situation ESS < E* never occurs, sampler execution is the same as
  # with the non-adaptive algorithm
  estar <- round(0.2*particleSize)
  oldExp <- 0
  curExp <- exponents[1]
  
  icount <- 1 # Iterations counter
  
  # All particles are given equal weight at initialization
  weights[1:length(weights)] <- log(1/particleSize)
  oldweights <- weights
  
  # Initial importance distribution
  importanceValues <- importanceDensity(particles)
  posteriorValues <- setup$posterior$density(particles) 
  #oldInter <- importanceValues
  oldExp <- 0
  
  if(!is.null(pars.lower) & !is.null(pars.lower)){
    sds <- 0.1 * (pars.upper - pars.lower)
  } else{
    sds <- rep(40,numPar)
  }
  if (is.null(proposal)){ proposalGenerator = createProposalGenerator(sds)}
  
  
  ##########################################################
  #                  Loop starts here                      #
  ##########################################################
  
  #for (i in 1:iterations){
  while(icount <= length(exponents)){
    
    #posteriorValues <- setup$posterior$density(particles) 
    #importanceValues <- importanceDensity(particles)
    
    if (numPar == 1) particles = matrix(particles, ncol = 1)
    
    print(c("icount",icount))
    print(c("posteriorValues",head(posteriorValues)))
    print(c("importanceValues",head(importanceValues)))
    
    # Using a while loop instead of for because in the adaptive algorithm, the number of iterations may increase
    
    ###########################################################
    # Mutate
    
    # 1) When to do the resampling
    # 2) HOW - currently adaptive Metropolis, could also do DEzs step 
    
    # Terminology issue: in literature, "resampling" generally describes the replication of particles with high likelihood, i.e. the previous
    # step. This step here is generally referred to as "mutate" or "move".
    # It may be necessary to rename the parameters "resampling" and "resamplingSteps", as well as info$resamplingAcceptance
    
    #if (resampling == T){
    
    
    ####------------- Replaced by mutate function - Start -------------------
    # if(mutate=="Metropolis"){
    #   
    #   if (adaptive == T){
    #     proposalGenerator = updateProposalGenerator(proposalGenerator, particles)
    #   }
    #   
    #   for (j in 1:resamplingSteps){
    #     particlesProposals = proposalGenerator$returnProposalMatrix(particles, scale = proposalScale)
    #     proposalPosteriors <- setup$posterior$density(particlesProposals) 
    #     proposalImportance <- importanceDensity(particlesProposals)
    #     
    #     #jumpProb <- exp(setup$posterior$density(particlesProposals) - likelihoodValues[sel])^(i/iterations) * exp(setup$prior$density(particlesProposals)   - setup$prior$density(particles))
    #     #jumpProb <- exp(setup$posterior$density(particlesProposals) - setup$posterior$density(particles)) * exp(setup$prior$density(particlesProposals)   - setup$prior$density(particles))
    #     jumpProb <- exp(proposalPosteriors - posteriorValues) * exp(setup$prior$density(particlesProposals)   - setup$prior$density(particles))
    #     #jumpProb <- exp(setup$posterior$density(particlesProposals) - posteriorValues) * exp(setup$prior$density(particlesProposals)   - setup$prior$density(particles))
    #     
    #     
    #     print(c("particlesProposals", head(particlesProposals)))
    #     print(c("particles", head(particles)))
    #     accepted <- jumpProb > runif(length(jumpProb), 0 ,1)
    #     rejectionRate = rejectionRate + sum(accepted)
    #     particles[accepted, ] = particlesProposals[accepted, ] 
    #     posteriorValues[accepted] <- proposalPosteriors[accepted]
    #     importanceValues[accepted] <- proposalImportance[accepted]
    #     info$resamplingAcceptance[(icount),j] <- sum(accepted)/particleSize
    #   }
    # } else if(mutate=="DE"){
    #   # Differential evolution
    #   # For now (test): basic DE algorithm
    #   
    #   for(j in 1:resamplingSteps){
    #     
    #     # Loop over particles
    #     for(part in 1:particleSize){
    #       particlesOld <- particles
    #       particleDiff <- rep(0,numPar)
    #       newParts <- rep(0L,2)
    #       
    #       # Sample 2 other particles
    #       while(all(particleDiff==0) | part %in% newParts){
    #         # The sampled particles might be identical (especially after resampling). Therefore, it is checked whether
    #         # the difference between particles is non-zero for at least one parameter. If the particles are identical,
    #         # the sampled particles are discarded and two new particles are sampled.
    #         # Also, the sampled particles should not include the current particle.
    #         newParts <- sample.int(n=particleSize,size=2)
    #         particleDiff <- particlesOld[newParts[1],] - particlesOld[newParts[2],]
    #       }
    #       
    #       particleProp <- particlesOld[part,] + (particleDiff * proposalScale) + runif(numPar,-b,b)
    #       jumpProb <- exp(setup$posterior$density(particleProp) - likelihoodValues[part]) * exp(setup$prior$density(particleProp)   - setup$prior$density(particles[part,]))
    #       
    #       if(jumpProb > runif(length(jumpProb), 0 ,1)){
    #         particles[part,] <- particleProp
    #         if(is.na(info$resamplingAcceptance[(icount-1),j])){
    #           info$resamplingAcceptance[(icount-1),j] <- 1/particleSize
    #         } else{
    #           info$resamplingAcceptance[(icount-1),j] <- info$resamplingAcceptance[(icount-1),j] + 1/particleSize
    #         }  
    #       } else{
    #         particles[part,] <- particlesOld[part,]
    #       }
    #     }
    #   }
    # }
    ####------------- Replaced by mutate function - End -------------------
    
    mutate.out <- mutate(setup = setup, particles = particles, proposalGenerator = proposalGenerator, posteriorValues = posteriorValues, importanceDensity = importanceDensity, method = mutate.method, steps = resamplingSteps, proposalScale = proposalScale, adaptive = adaptive)
    particles <- mutate.out$particles
    posteriorValues <- mutate.out$posteriorValues
    importanceValues <- mutate.out$importanceValues
    info$resamplingAcceptance[(icount),] <- mutate.out$acceptance
    
    
    
    # Reweighting
    
    # idea - adjust (1/iterations) such that always approx 30% of particles are maintain
    # level = sort(likelihoodValues)[acceptanceTarget]
    # best = likelihoodValues
    
    # Jeremia is using beta sequence 
    
    # https://github.com/florianhartig/BayesianTools/issues/23
    
    ## Update weights
    # Create intermediary distribution
    curExp <- exponents[icount]
    
    oldInter <- oldExp * posteriorValues + (1-oldExp) * importanceValues
    interDist <- curExp * posteriorValues + (1-curExp) * importanceValues
    
    # Calculate new weights
    
    weights <- weights + (interDist - oldInter)
    
    # Normalize (log-)weights so that the sum (of non-logs) equals 1

    weights <- weights - BayesianTools:::logSumExp(weights)
    
    # It may occur that the difference between the new and old weights of a particle is so great that an overflow occurs.
    # In this case, the algorithm reduces the difference between the current and next intermediary distribution (i.e. the difference
    # between the current and next value of the exponent). As a result, one more iteration is added.
    # This is the same that happens when the effective sample size becomes *very* low.
    #if(is.infinite(sum(sort(exp(weights))))){
    if(any(is.infinite(weights))){
      print("Infinite weights sum; retrying with new intermediary distribution")
      expdiff <- (exponents[icount] - exponents[(icount-1)]) * 0.5
      newExp <- exponents[(icount-1)] + expdiff
      exponents <- append(exponents,newExp,after=(icount-1))
      weights <- oldweights
      info$res.of[icount] <- TRUE
      next
    }
    
    print("-------")
    print(c("icount", icount))
    print(c("particles", head(particles)))
    print(c("importance", head(importanceValues)))
    print(c("posterior", head(posteriorValues)))
    print(c("interDist", head(interDist)))
    print(c("oldInter", head(oldInter)))
    print(c("weights", head(weights)))
    print(c("oldweights", head(oldweights)))
    
    # posterior = setup$posterior$density(particles, returnAll = T)
    # likelihoodValues <- posterior[,2]
    #     ll = likelihoodValues - max(likelihoodValues, na.rm = T)
    # llCutoff = sort(ll)[acceptanceTarget] 
    # relativeL = exp(likelihoodValues - max(likelihoodValues, na.rm = T))^(1/iterations)
    
    
    #########################################################
    # Re?sampling step 
    
    # Option 1) Always resample (default BT)
    # Option 2) Keep weights and resample if ESS < something
    # => for now, option 1 is activated if the argument ess.limit is set to a value >= the number of particles
    
    # Question: last step complete resampling, so that all particles are the same?
    
    #sel = sample.int(n=length(likelihoodValues), size = length(likelihoodValues), replace = T, prob = relativeL)
    
    # Calculate effective sample size (ESS)
    ess <- 1 / sum(exp(2 * weights))

    if(is.infinite(ess)){print(weights)}
    
    # Determine if resampling is necessary - if yes, resample with probability given by the weight
    # Resample also on the last iteration, or at the iteration (starting from the end) given by the parameter lastResample (output is based
    # on the location of particles in parameter space, the weights are not considered in the output)
    
    print(c("ess", ess))
    print(c("ess.limit", ess.limit))
    print(c("curExp", curExp))

    if(ess < ess.limit | icount == (length(exponents) - lastResample)){
      
      
      # Determine if exponents need to be modified
      if(ess < estar){
        # If yes:
        # - discard newly sampled population of particles
        # - do not increase counter
        # - modify x
        
        expdiff <- (exponents[icount] - exponents[(icount-1)]) * 0.5
        newExp <- exponents[(icount-1)] + expdiff
        print(c("exponents[icount]",exponents[icount]))
        print(c("exponents[icount-1]",exponents[icount-1]))
        print(c("expdiff", expdiff))
        print(c("newExp",newExp))
        
        weights <- oldweights
        
        if(icount==1){
          newExp <- exponents[icount] * 0.5
          exponents <- c(newExp,exponents)
        } else { 
          exponents <- append(exponents,newExp,after=(icount-1))
        }
        
        print(c("length(exponents", length(exponents)))
        print(c("new exponents", exponents))
        info$res.ess[icount] <- TRUE
        next
      } 
      
      oldExp <- curExp

      sel = resample(weights, method = sampling)
            
      particles = particles[sel,]
      posteriorValues <- posteriorValues[sel]
      importanceValues <- importanceValues[sel]
      info$survivingParticles[icount] = length(unique(sel))
      
      # Set all weights equal
      weights[1:length(weights)] <- log(1/particleSize)
      # Normalize (log-)weights so that the sum (of non-logs) equals 1
      weights <- weights - BayesianTools:::logSumExp(weights)
      
      # Get posterior values and distribution for new particles
      #oldInter <- curExp * setup$posterior$density(particles) 
      #          + (1-curExp) * importanceDensity(particles)
      #oldInter <- curExp * posteriorValues + (1-curExp) * importanceValues
      
      print("-- Resampling")
      print(c("posterior", head(posteriorValues)))
      print(c("importance", head(importanceValues)))
      #print(c("oldInter", head(oldInter)))
      

    } else{
      # Save intermediate distribution for the next iteration
      #oldInter <- interDist
    }
    
    oldweights <- weights
    
    #print(c("oldInter saved", head(oldInter)))
    
    info$ess.vec[icount] <- ess
    
    icount <- icount + 1    
    #curExp <- m * (icount*200/iterations)^x
    oldExp <- curExp
    
    if(!is.null(diagnostics)) info$diagnostics[[icount-1]] <- diagnostics(particles)
  }
  
  # Last mutation, to increase diversity between particles
  mutate.out <- mutate(setup = setup, particles = particles, proposalGenerator = proposalGenerator, posteriorValues = posteriorValues, importanceDensity = importanceDensity, method = mutate.method, steps = lastMutateSteps, proposalScale = proposalScale, adaptive = adaptive)
  particles <- mutate.out$particles
  posteriorValues <- mutate.out$posteriorValues
  importanceValues <- mutate.out$importanceValues
  
  info$rejectionRate = rejectionRate / (iterations * resamplingSteps)
  
  settings = list(initialParticles = initialParticles, proposalGenerator = proposalGenerator)
  
  out = list(
    setup = setup,
    settings = settings,
    particles = particles,
    posterior = posterior,
    info = info
  )
  
  class(out) <- c("smcSampler", "bayesianOutput")
  return(out)
  
}


#' Residual Resampling
#' 
#' @keywords internal
residualResampling <- function(weights){
  weights <- exp(weights)
  # Define number of replications for each particle (integer), as well as residuals
  nrep <- weights * length(weights)
  nrep.int <- floor(nrep)
  residuals <- nrep - nrep.int
  
  # Resample remaining particles based on residuals
  length.missing <- length(weights) - sum(nrep.int)
  missing <- vector("numeric", length.missing)
  missing <- sample.int(n=length(weights), size=length.missing, replace = TRUE, prob = residuals)
  
  # Return the indices of particles to be resampled
  new.parts <- rep(1:length(weights), nrep.int)
  new.parts <- c(new.parts, missing)
  return(new.parts)
}

systematicResampling <- function(weights){
  #print(weights)
  weights <- exp(weights)
  # Reorder weights in increasing order
  rank.weights <- rank(weights, ties.method = "first")
  sort.weights <- sort(weights)
  cumu.weights <- c(0,cumsum(sort.weights)[1:(length(weights)-1)])
  #print(cumu.weights)
  n.parts <- length(weights)  # Number of particles
  
  u <- runif(n=1, min=0, max=1/n.parts)
  U <- vector("numeric", n.parts)
  new.parts <- vector("numeric", n.parts)
  new.parts.sort <- vector("numeric", n.parts)
  for(i in 1:n.parts){
    U[i] <- ((i-1)/n.parts) + u
    #print(c(i, U[i]))
    
    new.parts.sort[i] <- tail(which(cumu.weights <= U[i]),1)
    #print(c("new.parts.sort", new.parts.sort[i]))
    # Match sorted ranks back to original order of weights
    new.parts[i] <- which(rank.weights==new.parts.sort[i])
  }
  
  return(new.parts)
}



resample <- function(weights, method = "multinomial"){
  
  particleSize = length(weights)
  
  if(method == "multinomial"){
    sel <- sample.int(n=particleSize, size = particleSize, replace = T, prob = exp(weights))
  } else if(method == "residual"){
    sel <- residualResampling(weights)
  } else if(method == "systematic"){
    sel <- systematicResampling(weights)
  } else{
    stop("Invalid string for resampling argument")
  }
  return(sel)
}

beta.search <- function(ess, target.ess, posteriorValues, importanceValues, oldInter, curWeights, curExp, tol=0.001){
  # A function to dynamically set the next exponent to build the next intermediary distribution.
  # Uses the bisection method. Following Jasra et al. (2001), Scand J Statist, doi: 10.1111/j.1467-9469.2010.00723.x
  
  # Initial exponent - set to 1 (maximum possible value)
  tryDiff <- 1 - curExp
  
  while(abs(ess-target.ess) > tol){
    tryExp <- curExp + tryDiff
    
    tryDist <- tryExp * posteriorValues + (1-tryExp) * importanceValues
    
    tryWeights <- curWeights + (tryDist - oldInter)
    # Normalize (log-)weights so that the sum (of non-logs) equals 1
    tryWeights <- tryWeights - BayesianTools:::logSumExp(tryWeights)
    
    try.ess <- 1 / sum(exp(2 * weights))
    
    if(try.ess - target.ess > tol){
      # Greater ESS than desired -> choose a *larger* exponent in the next iteration (efficiency)
      tryDiff <- tryDiff + tryDiff * 0.5
      
      # Exponent may not be greater than one
      if(tryDiff + curExp >= 1){
        tryExp <- 1
        break
      }
      
    } else if(try.ess - target.ess < -tol | any(is.infinite(tryWeights))){
      # Smaller ESS than desired -> choose a *smaller* exponent in the next iteration (stability)
      # Also includes a failsafe in the case of infinite weights (numerical issue)
      tryDiff <- tryDiff - (tryDiff * 0.5)
      
      # Exponent may not decrease from previous value
      if(tryDiff < 0){
        # This should not occur - leaving it for testing only
        stop("problem defining next exponent")
      }
    }
  }
  
  out <- list(newExp=tryExp, weights=tryWeights, interDist=tryDist, ess=try.ess)
}


mutate <- function(setup, particles, proposalGenerator, posteriorValues, importanceDensity, method, steps, proposalScale, adaptive = TRUE){
  if(is.vector(particles)){particles = matrix(particles, ncol = 1)}
  acceptance <- vector("numeric", length=steps)
  importanceValues <- importanceDensity(particles)
  
  if(method=="Metropolis"){
    if(adaptive){
      proposalGenerator = updateProposalGenerator(proposalGenerator, particles)
    }
    
    for(j in 1:steps){
      particlesProposals = proposalGenerator$returnProposalMatrix(particles, scale = proposalScale)
      proposalPosteriors <- setup$posterior$density(particlesProposals) 
      proposalImportance <- importanceDensity(particlesProposals)
      
      #jumpProb <- exp(setup$posterior$density(particlesProposals) - likelihoodValues[sel])^(i/iterations) * exp(setup$prior$density(particlesProposals)   - setup$prior$density(particles))
      #jumpProb <- exp(setup$posterior$density(particlesProposals) - setup$posterior$density(particles)) * exp(setup$prior$density(particlesProposals)   - setup$prior$density(particles))
      jumpProb <- exp(proposalPosteriors - posteriorValues) * exp(setup$prior$density(particlesProposals)   - setup$prior$density(particles))
      #jumpProb <- exp(setup$posterior$density(particlesProposals) - posteriorValues) * exp(setup$prior$density(particlesProposals)   - setup$prior$density(particles))
      
      
      print(c("particlesProposals", head(particlesProposals)))
      print(c("particles", head(particles)))
      accepted <- jumpProb > runif(length(jumpProb), 0 ,1)
      #rejectionRate = rejectionRate + sum(accepted)
      particles[accepted, ] = particlesProposals[accepted, ] 
      posteriorValues[accepted] <- proposalPosteriors[accepted]
      importanceValues[accepted] <- proposalImportance[accepted]
      acceptance[j] <- sum(accepted)/nrow(particles)
      #info$resamplingAcceptance[(icount),j] <- sum(accepted)/particleSize
    }
  }
  out <- list(particles=particles, posteriorValues = posteriorValues, importanceValues = importanceValues, acceptance = acceptance)
}
