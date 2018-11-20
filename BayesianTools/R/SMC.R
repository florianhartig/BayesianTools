#' SMC sampler
#' @author Florian Hartig, Matthias Speich
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
                       ess.factor = 0.95,           # TODO document
                       lastResample = 1,
                       pars.lower=NULL, 
                       pars.upper=NULL, 
                       mutate.method ="Metropolis", # TODO document
                       b=1e-04,                     # TODO document
                       diagnostics = NULL,
                       reference=NULL){
  
  ########### SETUP STEPS ########################
  
  # Timing: normally done in mcmcRun.R. Timing is included in this function for the purpose of testing the SMC function.
  # Starting the clock
  ptm <- proc.time()
  
  #if(resamplingSteps < 1) stop("SMC error, resamplingSteps can't be < 1")
  
  info = list()
  # The number of iterations is not known at the beginning. Therefore, output vectors (with one value for each iteration)
  # are made large enough that they will probably not need to be grown (growing arrays can be slow).
  if(resamplingSteps > 0) info$resamplingAcceptance = as.data.frame(matrix(nrow = 10000, ncol = resamplingSteps)) #Using DF instead of matrix because number of iterations may change due to adaptive algorithm
  info$survivingParticles = rep(NA, 10000)
  info$ess.vec <-  info$exponents <- vector("numeric", 10000)
  info$diagnostics <- diag.end <- list()
  lastAccept <- vector("numeric", lastMutateSteps)


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
    if(class(initialParticles$particles) == "numeric"){
      initialParticles = as.matrix(initialParticles$particles,ncol=1)
    } else if(class(initialParticles$particles) == "data.frame"){
      initialParticles = as.matrix(initialParticles$particles)
    } else if(class(initialParticles$particles) == "matrix"){
      initialParticles <- initialParticles$particles
    }
  }
  
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
  
  #### TEST
  #proposalScale <- runif(nrow(particles), min = 1e-06, max = 1)
  Z <- particles

  
  ########### WEIGHT SETUP ########################
    
  # Define an exponential sequence of beta parameters, i.e. the exponents to weight the likelihood against prior.
  # Idea and parameter values from Jeremiah et al. (2012), Environ Modell Softw
  # if(is.null(exponents)){
  #   iters <- seq(0,iterations)
  #   exponents <- m * (iters*200/iterations)^x
  #   exponents <- pmin(1,exponents)
  #   exponents <- exponents[2:length(exponents)]
  # }
  
  # Initial value for exponent: here the first value of the series is taken.
  # If the situation ESS < E* never occurs, sampler execution is the same as
  # with the non-adaptive algorithm
  estar <- round(0.2*particleSize)
  oldExp <- 0
  #curExp <- exponents[1]
  curExp <- 0 # For fully adaptive algorithm: initial value of curExp does not matter, but must be < 1
  
  icount <- 1 # Iterations counter
  
  # All particles are given equal weight at initialization
  weights[1:length(weights)] <- log(1/particleSize)
  oldweights <- weights
  
  # Initial importance distribution
  importanceValues <- importanceDensity(particles)
  #print("Calculating initial posteriors")
  posteriorValues <- setup$posterior$density(particles) 
  #print("Initial posteriors done")
  #oldInter <- importanceValues
  oldExp <- 0
  
  if(!is.null(pars.lower) & !is.null(pars.lower)){
    sds <- 0.1 * (pars.upper - pars.lower)
  } else{
    sds <- rep(40,numPar)
  }
  if (is.null(proposal)){ proposalGenerator = createProposalGenerator(sds)}
  
  lastIteration <- FALSE
  
  ## Estimate minimum and maximum value of each parameter (from initial population)
  min.pars <- apply(particles,2,min)
  max.pars <- apply(particles,2,max)
  
  
  ##########################################################
  #                  Loop starts here                      #
  ##########################################################
  
  #for (i in 1:iterations){
  # while(icount <= length(exponents)){
  #while(curExp <= 1 & !lastIteration){ # For fully adaptive algorithm
  #  if(curExp == 1){lastIteration <- TRUE}
  while(curExp < 1){
    
    ### TEMPORARY DEBUG
    if(sum(is.nan(posteriorValues)) > 0 | sum(is.na(posteriorValues)) > 0 | sum(is.infinite(posteriorValues)) > 0){
      print("NaN / NA / Infinite in importanceValues")
      if(sum(is.nan(posteriorValues)) > 1){
        print(c("NaN indices", which(is.nan(posteriorValues))))
        print("Particles with NaN:")
        print(particles[is.nan(posteriorValues),])
      }
      
      if(sum(is.na(posteriorValues)) > 1){
        print(c("NA indices", which(is.na(posteriorValues))))
        print("Particles with NA:")
        print(particles[is.na(posteriorValues),])
      }
      
      if(sum(is.infinite(posteriorValues)) > 1){
        print(c("Infinite indices", which(is.infinite(posteriorValues))))
        print("Particles with Infinite:")
        print(particles[is.infinite(posteriorValues),])
      }
    }
    
    # Update particle min/max (if a particle is more "extreme" than previously recorded)
    cur.min <- apply(particles,2,min)
    cur.max <- apply(particles,2,max)
    min.pars <- pmin(min.pars, cur.min)
    max.pars <- pmax(max.pars, cur.max)

    if (numPar == 1) particles = matrix(particles, ncol = 1)
    
    # Using a while loop instead of for because in the adaptive algorithm, the number of iterations may increase
    
    ###########################################################
    # Mutate
    
    # 1) When to do the resampling
    # 2) HOW - currently adaptive Metropolis, could also do DEzs step 
    
    # Terminology issue: in literature, "resampling" generally describes the replication of particles with high likelihood, i.e. the previous
    # step. This step here is generally referred to as "mutate" or "move".
    # It may be necessary to rename the parameters "resampling" and "resamplingSteps", as well as info$resamplingAcceptance
    
    
    # mutate.out <- mutate(setup = setup, particles = particles, proposalGenerator = proposalGenerator, posteriorValues = posteriorValues, importanceDensity = importanceDensity, method = mutate.method, steps = resamplingSteps, proposalScale = proposalScale, adaptive = adaptive, b=b)
    # particles <- mutate.out$particles
    # posteriorValues <- mutate.out$posteriorValues
    # importanceValues <- mutate.out$importanceValues
    # info$resamplingAcceptance[(icount),] <- mutate.out$acceptance
    
    # Reweighting
    
    # idea - adjust (1/iterations) such that always approx 30% of particles are maintain
    # level = sort(likelihoodValues)[acceptanceTarget]
    # best = likelihoodValues
    
    # Jeremia is using beta sequence 
    
    # https://github.com/florianhartig/BayesianTools/issues/23
    

    ## Update weights
    # Create intermediary distribution
    ess <- 1 / sum(exp(2 * weights))
    oldInter <- oldExp * posteriorValues + (1-oldExp) * importanceValues
    
    #inter.out <- beta.search(ess=ess, target.ess = (ess * ess.factor), posteriorValues = posteriorValues, importanceValues = importanceValues, oldInter = oldInter, curWeights = weights, curExp = curExp)
    inter.out <- beta.search(ess=ess, target.ess = (ess.limit-1), posteriorValues = posteriorValues, importanceValues = importanceValues, oldInter = oldInter, curWeights = weights, curExp = curExp)
    curExp <- inter.out$newExp
    weights <- inter.out$weights
    interDist <- inter.out$interDist
    ess <- inter.out$ess
    doResample <- inter.out$doResample
    

    
    info$ess.vec[icount] <- ess
    info$exponents[icount] <- curExp
    if(abs(1-curExp) < 1e-03){
      curExp <- 1
    }
    
    #########################################################
    # Re?sampling step 

    # Determine if resampling is necessary - if yes, resample with probability given by the weight
    # Resample also on the last iteration, or at the iteration (starting from the end) given by the parameter lastResample (output is based
    # on the location of particles in parameter space, the weights are not considered in the output)

    if(ess < ess.limit | icount == (length(exponents) - lastResample)  | doResample){
 
      oldExp <- curExp

      sel = resample(weights, method = sampling)
      
      particles = particles[sel,]
      
      #print(apply(particles, 2, sd))
      #print(apply(particles, 2, sd)/apply(particles, 2, mean))

      posteriorValues <- posteriorValues[sel]
      importanceValues <- importanceValues[sel]
      info$survivingParticles[icount] = length(unique(sel))
      
      # Set all weights equal
      weights[1:length(weights)] <- log(1/particleSize)
      # Normalize (log-)weights so that the sum (of non-logs) equals 1
      weights <- weights - BayesianTools:::logSumExp(weights)
      


      
      ## Check diversity of particles, and require a mutate step if necessary
      # Scale particles 
      #print(c("nrow unique particles", nrow(unique(particles))))
      #print(c("curExp", curExp))
      #if(nrow(unique(particles)) < 1000){
      scaled.particles <- t(apply(particles, 1, scale.particles, min.pars = min.pars, max.pars = max.pars))
      #diversity <- rao.div(scaled.particles = scaled.particles)
      #print(c("Diversity", diversity))
      #}
      
      if(resamplingSteps > 0){
        mutate.out <- mutate(setup = setup, particles = particles, proposalGenerator = proposalGenerator, posteriorValues = posteriorValues, importanceDensity = importanceDensity, method = mutate.method, steps = resamplingSteps, proposalScale = proposalScale, adaptive = adaptive, b=b, min.pars = min.pars, max.pars = max.pars, Z = NULL)
        particles <- mutate.out$particles
        posteriorValues <- mutate.out$posteriorValues
        importanceValues <- mutate.out$importanceValues
        info$resamplingAcceptance[(icount),] <- mutate.out$acceptance
        if(length(proposalScale) > 1){
          proposalScale <- mutate.out$proposalScale
        }
        
        #plot(density(particles[,1]), xlim=c(0.314, 0.681), main = c(curExp, nrow(unique(particles))))
        
        #Z <- rbind(Z, particles)
        
        # while(nrow(unique(particles)) < nrow(particles) * 0.5){
        #   mutate.out <- mutate(setup = setup, particles = particles, proposalGenerator = proposalGenerator, posteriorValues = posteriorValues, importanceDensity = importanceDensity, method = mutate.method, steps = resamplingSteps, proposalScale = proposalScale, adaptive = adaptive, b=b, min.pars = min.pars, max.pars = max.pars, Z = NULL)
        # 
        #   particles <- mutate.out$particles
        #   posteriorValues <- mutate.out$posteriorValues
        #   importanceValues <- mutate.out$importanceValues
        #   info$resamplingAcceptance[(icount),] <- mutate.out$acceptance
        #   if(length(proposalScale) > 1){
        #     proposalScale <- mutate.out$proposalScale
        #   }
        # 
        #   plot(density(particles[,1]), xlim=c(0.314, 0.681), main = c("oi", curExp, nrow(unique(particles))))
        # }
          
      }
      
    } 
    
    ## Mutate step if required because of nuerical issues or low particle diversity
    # if(doResample){
    #   mutate.out <- mutate(setup = setup, particles = particles, proposalGenerator = proposalGenerator, posteriorValues = posteriorValues, importanceDensity = importanceDensity, method = mutate.method, steps = resamplingSteps, proposalScale = proposalScale, adaptive = adaptive, b=b)
    #   particles <- mutate.out$particles
    #   posteriorValues <- mutate.out$posteriorValues
    #   importanceValues <- mutate.out$importanceValues
    #   info$resamplingAcceptance[(icount),] <- mutate.out$acceptance
    # }
    

    
    oldweights <- weights
    info$ess.vec[icount] <- ess
    
    icount <- icount + 1    
    oldExp <- curExp
    
    if(!is.null(diagnostics)) info$diagnostics[[icount-1]] <- diagnostics(particles, reference)
  }
  
  if (numPar == 1) particles = matrix(particles, ncol = 1)
  
  # Last resampling step, so that particles are distributed according to target distribution
  sel = resample(weights, method = sampling)
  particles = particles[sel,]
  info$survivingParticles[(icount-1)] <- length(unique(sel))
  
  #print("last")
  #print(apply(particles, 2, sd)/apply(particles, 2, mean))
  
  if(resamplingSteps > 0){
    for(mutateStep in 1:lastMutateSteps){
      # Last mutation, to increase diversity between particles
      mutate.out <- mutate(setup = setup, particles = particles, proposalGenerator = proposalGenerator, posteriorValues = posteriorValues, importanceDensity = importanceDensity, method = mutate.method, steps = 1, proposalScale = proposalScale, adaptive = adaptive, b=b, min.pars = min.pars, max.pars = max.pars, Z = NULL)
      particles <- mutate.out$particles
      posteriorValues <- mutate.out$posteriorValues
      importanceValues <- mutate.out$importanceValues
      lastAccept[mutateStep] <- mutate.out$acceptance
      if(!is.null(diagnostics)) diag.end[[mutateStep]] <- diagnostics(particles, reference)
      
      #plot(density(particles[,1]), xlim=c(0.314, 0.681))
      
      #Z <- rbind(Z, particles)
    }

  }
  info$rejectionRate = rejectionRate / (iterations * resamplingSteps)
  
  # Timing: normally done in mcmcRun.R. Timing is included in this function for the purpose of testing the SMC function.
  # Stopping the clock
  elapsed.time <- proc.time() - ptm
  
  # Trim info objects (as their size = number of iterations were not known beforehand, they
  # were allocated a large vector/matrix). Trim up to (icount-1), as icount got incremented
  # at the end of the last iteration.
  if(resamplingSteps > 0) info$resamplingAcceptance <- info$resamplingAcceptance[1:(icount-1),]
  info$survivingParticles <- info$survivingParticles[1:(icount-1)]
  info$exponents <- info$exponents[1:(icount-1)]
  #info$diagnostics <- c(info$diagnostics[[1:(icount-1)]], diag.end)
  info$diagnostics <- c(info$diagnostics, diag.end)
  info$ess.vec <- info$ess.vec[1:(icount-1)]
  if(resamplingSteps > 0) info$lastAccept <- lastAccept
  info$elapsed.time <- elapsed.time
  
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

######################################################
# END of SMC function
######################################################

######################################################
# Auxiliary functions for resampling

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
  weights <- exp(weights)
  # Reorder weights in increasing order
  rank.weights <- rank(weights, ties.method = "first")
  sort.weights <- sort(weights)
  cumu.weights <- c(0,cumsum(sort.weights)[1:(length(weights)-1)])
  n.parts <- length(weights)  # Number of particles
  
  u <- runif(n=1, min=0, max=1/n.parts)
  U <- vector("numeric", n.parts)
  new.parts <- vector("numeric", n.parts)
  new.parts.sort <- vector("numeric", n.parts)
  for(i in 1:n.parts){
    U[i] <- ((i-1)/n.parts) + u

    
    new.parts.sort[i] <- tail(which(cumu.weights <= U[i]), 1)
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

######################################################

# Auxiliary function for adaptive algorithm

beta.search <- function(ess, target.ess, posteriorValues, importanceValues, oldInter, curWeights, curExp, tol=1){
  # A function to dynamically set the next exponent to build the next intermediary distribution.
  # Uses the bisection method. Following Jasra et al. (2011), Scand J Statist, doi: 10.1111/j.1467-9469.2010.00723.x
  
  # Initial exponent - set to 1 (maximum possible value)
  tryDiff <- 1 - curExp
  a <- curExp
  b <- 1
  try.ess <- target.ess + (100*tol) # Dummy initial value that is sure to be different from target value
  tryExp <- 0

  
  while(abs(try.ess-target.ess) > tol & tryExp < 1 & abs(b-a) > 1e-10){
    tryExp <- (a+b) * 0.5

    
    tryDist <- tryExp * posteriorValues + (1-tryExp) * importanceValues
    
    tryWeights <- curWeights + (tryDist - oldInter)
    # Normalize (log-)weights so that the sum (of non-logs) equals 1
    
    ### DEBUG
    if(sum(is.na(tryWeights)) > 0){
      print("DEBUG")
      print("tryWeights")
      print(tryWeights)
      print("tryDist")
      print(tryDist)
      print("oldInter")
      print(oldInter)
      print(c("curExp", curExp))
      print(c("tryExp", tryExp))
      print(c("a", a))
      print(c("b", b))
      print("posteriorValues")
      print(posteriorValues)
      print("importanceValues")
      print(importanceValues)

    }

    tryWeights <- tryWeights - BayesianTools:::logSumExp(tryWeights)
    try.ess <- 1 / sum(exp(2 * tryWeights))
    
    if(try.ess - target.ess > tol & !any(is.infinite(tryWeights))){
      # Greater ESS than desired -> choose a *larger* exponent in the next iteration (efficiency)
      a <- tryExp
    } else if(try.ess - target.ess < -tol | any(is.infinite(tryWeights))){
      # Smaller ESS than desired -> choose a *smaller* exponent in the next iteration (stability)
      # Also includes a failsafe in the case of infinite weights (numerical issue)
      b <- tryExp
    }
  }


  
  if(any(is.infinite(tryWeights)) | is.infinite(try.ess)){
    doResample <- TRUE
    tryExp <- curExp
    tryWeights <- curWeights
    tryDist <- oldInter
    try.ess <- ess
  } else{
    doResample <- FALSE
  }
  
  out <- list(newExp=tryExp, weights=tryWeights, interDist=tryDist, ess=try.ess, doResample = doResample)
  return(out)
}

######################################################
# Auxiliary function for particle mutation


mutate <- function(setup, particles, proposalGenerator, posteriorValues, importanceDensity, method, steps, proposalScale, adaptive = TRUE, b=1E-04, min.pars=NULL, max.pars=NULL, Z = NULL){

  if(is.vector(particles)){particles = matrix(particles, ncol = 1)}
  acceptance <- vector("numeric", length=steps)
  importanceValues <- importanceDensity(particles)
  
  
  if(length(proposalScale) == 1){
    scale.factors <- rep(proposalScale, nrow(particles))
  } else {
    scale.factors <- proposalScale
    scaled.particles <- t(apply(particles, 1, scale.particles, min.pars = min.pars, max.pars = max.pars))
  }
  
  if(is.null(Z)){
    Z <- particles
  }
  
  r_extra <- vector("numeric", length = nrow(particles))

  
    if(adaptive){
      proposalGenerator = updateProposalGenerator(proposalGenerator, particles)
    }
    
    for(j in 1:steps){
      
      if(method=="Metropolis"){
        particlesProposals = proposalGenerator$returnProposalMatrix(particles, scale = proposalScale)
      } else if(method=="DE"){

        particlesProposals <- particles
        
        for(part in 1:nrow(particles)){
          
          if(runif(1) < 0.1) {
            ## Snooker update
            particlesOld <- particles
            particleDiff12 <- rep(0,ncol(particles))
            particleDiff23 <- rep(0,ncol(particles))
            particleDiff13 <- rep(0,ncol(particles))
            newParts <- rep(0L,3)
            
            # Sample 3 other particles
            #while(all(particleDiff12==0) | all(particleDiff23==0) | all(particleDiff13==0) | part %in% newParts){
            while(all(particleDiff12==0) | all(particleDiff23==0) | all(particleDiff13==0)){
              # The sampled particles might be identical (especially after resampling). Therefore, it is checked whether
              # the difference between particles is non-zero for at least one parameter. If the particles are identical,
              # the sampled particles are discarded and two new particles are sampled.
              # Also, the sampled particles should not include the current particle.
              newParts <- sample.int(n=nrow(particles),size=3, replace=FALSE)
              particleDiff12 <- Z[newParts[1],] - Z[newParts[2],]
              particleDiff23 <- Z[newParts[2],] - Z[newParts[3],]
              particleDiff13 <- Z[newParts[1],] - Z[newParts[3],]
            }
            
            x_z <- particles[part,] - Z[newParts[3],]
            D2 <- max(sum(x_z*x_z), 1.0e-300)
            projdiff <- sum((Z[newParts[1],] -Z[newParts[2],]) * x_z)/D2 # inner_product of difference with x_z / squared norm x_z
            particlesProposals[part,] <- particles[part,] + runif(1, min=1.2, max=2.2) * projdiff * x_z
            
            x_z <- particlesProposals[part,] - Z[newParts[3],]
            D2prop <- max(sum(x_z*x_z), 1.0e-300)
            r_extra[part] <- ((ncol(particles)-1)/2) * (log(D2prop) - log(D2))
            
            #if(part %in% 1:10){
            #  print(c("Part snooker", part, r_extra[part]))
            #}
            
          } else{ 
            ## No snooker
            particlesOld <- particles
            particleDiff <- rep(0,ncol(particles))
            newParts <- rep(0L,2)
            
            # Sample 2 other particles
            #while(all(particleDiff==0) | part %in% newParts){
            while(all(particleDiff==0)){
              # The sampled particles might be identical (especially after resampling). Therefore, it is checked whether
              # the difference between particles is non-zero for at least one parameter. If the particles are identical,
              # the sampled particles are discarded and two new particles are sampled.
              # Also, the sampled particles should not include the current particle.
              newParts <- sample.int(n=nrow(particles),size=2, replace=FALSE)
              particleDiff <- Z[newParts[1],] - Z[newParts[2],]
            }
            numPar <- ncol(particles)
            randVector <- runif(numPar,-b,b)
            particlesProposals[part,] <-   particles[part,] + (particleDiff * scale.factors[part]) + randVector  
            r_extra[part] <- 0
            
            # if(part %in% 1:6){
            #   print(c("Part ", part))
            #   print("Particle")
            #   print(particlesOld[part,])
            #   print("Proposal")
            #   print(particlesProposals[part,])
            # }
            #if(part %in% 1:10){
            #  print(c("Part NO snooker", part, r_extra[part]))
            #}
          }
        }  
      }
      

      proposalPosteriors <- setup$posterior$density(particlesProposals) 
      proposalImportance <- importanceDensity(particlesProposals)
      
      #print(c("Priors", head(setup$prior$density(particlesProposals), 10)))
      #print(c("Proposals post", head(proposalPosteriors, 10)))
      #print(c("Particles post", head(posteriorValues, 10)))
      
      #print(c("Infinite", sum(is.infinite(proposalPosteriors))))
      #print(c("Infinite prior", sum(is.infinite(setup$prior$density(particlesProposals)))))
      
      #print(c("Minus Infinite", sum(proposalPosteriors== -Inf)))
      #print(c("Minus Infinite prior", sum(setup$prior$density(particlesProposals)== -Inf)))
      
      ### TODO: why is jumpProb not simply the ratio of posteriors, i.e. why include the prior?
      # This has no influence on the current experiments (with flat priors), but needs to be clarified before other uses.

      #jumpProb <- exp(proposalPosteriors - posteriorValues) * exp(setup$prior$density(particlesProposals)   - setup$prior$density(particles))
      jumpProb <- exp(proposalPosteriors - posteriorValues + r_extra)

      accepted <- jumpProb > runif(length(jumpProb), 0 ,1)
      
      if(sum(is.na(accepted)) > 1){
        print("sum is.na accepted", sum(is.na(accepted)))
        print("proposalPosteriors")
        print(proposalPosteriors)
        print("posteriorValues")
        print(posteriorValues)
        print("r_extra")
        print(r_extra)
        
      }
      
      particles[accepted, ] = particlesProposals[accepted, ] 
      posteriorValues[accepted] <- proposalPosteriors[accepted]
      importanceValues[accepted] <- proposalImportance[accepted]
      acceptance[j] <- sum(accepted)/nrow(particles)
      
      ## Adaptation of scaling factors following Fearnhead & Taylor (2013)
      if(length(proposalScale) > 1){
      
        ### Calculate ESJD for each proposal
        ##Scale particles
        scaled.proposals <- t(apply(particlesProposals, 1, scale.particles, min.pars = min.pars, max.pars = max.pars))
        
        
        #print("scaled.particles")
        #print(head(scaled.particles))
        
        #print("scaled.proposals")
        #print(head(scaled.proposals))
        
        
        #esjd <- vector("numeric", nrow(particles))
        esjd <- vector("numeric", sum(r_extra == 0))
        loop.ind <- 0
        ind.nosnooker <- which(r_extra == 0)
        for(part in 1:nrow(particles)){
          if(r_extra[part] == 0) {
            loop.ind <- loop.ind + 1
            esjd[loop.ind] <- dist(rbind(scaled.particles[part,], scaled.proposals[part,])) * pmin(1,jumpProb[part])
            }
          #print(rbind(scaled.particles[part,], scaled.proposals[part,]))
          #if(r_extra[part] != 0){print("Snooker")}
          #print(c("dist, jumpProb, post, prop",dist(rbind(scaled.particles[part,], scaled.proposals[part,])), jumpProb[part], posteriorValues[part], proposalPosteriors[part]))
        }
        
        #print(c("length esjd", length(esjd)))
        
        esjd[is.infinite(esjd)] <- 0
        #print(c("esjd", head(esjd)))
        #print(c("max esjd", max(esjd)))
        
        # Normalize
        esjd <- esjd/sum(esjd)
        #print(c("esjd norm", head(esjd)))
        #print(c("max esjd norm", max(esjd)))
        
        # Resample and add noise
        sel <- resample(log(esjd), "systematic")
       #print(c("sel", head(sel)))
       #print(c("unique sel", unique(sel)))
        if(length(unique(sel)) <= 5){
          scale.factors <- scale.factors #+ rnorm(length(scale.factors), 0, 0.015)
        } else{
          scale.factors <- scale.factors[ind.nosnooker[sel]] #+ rnorm(length(scale.factors), 0, 0.015)
        }
        
        # Draw scale factors from normal distribution
        mean.sf <- mean(scale.factors)
        sd.sf <- sd(scale.factors)
        
        scale.factors <- rnorm(n = nrow(particles), mean = mean.sf, sd = max(0.05,sd.sf))
        
        scale.factors[scale.factors <= 0] <- 1E-06
        
      }
      
      #print(head(scale.factors, 100))
      #plot(density(scale.factors), xlim=c(0,1))
      
    }
  #}
  if(length(proposalScale) == 1){
    out <- list(particles=particles, posteriorValues = posteriorValues, importanceValues = importanceValues, acceptance = acceptance)
  } else{
    out <- list(particles=particles, posteriorValues = posteriorValues, importanceValues = importanceValues, acceptance = acceptance, proposalScale = scale.factors)
  }
  
  return(out)
}

scale.particles <- function(particle.row, min.pars, max.pars){
  scaled.particles <- (particle.row - min.pars) / (max.pars - min.pars)
  return(scaled.particles)
}

rao.div <- function(scaled.particles){
  # Determine unique particles, their sum and pairwise difference
  unique.particles <- unique(scaled.particles)
  unique.count <- vector("numeric", nrow(unique.particles))
  for(i in 1:nrow(unique.particles)){
    unique.count[i] <- sum(apply(scaled.particles, 1, function(x, cur.unique){
      if(sum(x-cur.unique)==0){
        return(TRUE)
      } else {
        return(FALSE)
      }
    }, cur.unique = unique.particles[i,]))
  }
  
  unique.count <- unique.count/nrow(scaled.particles)
  
  total.dist <- 0
  
  for(i in 1:(nrow(unique.particles)-1)){
    for(j in (i+1):nrow(unique.particles)){
      #total.dist <- total.dist + dist(unique.particles[c(i,j),]) * unique.count[i] * unique.count[j]
      total.dist <- total.dist + sum(abs(unique.particles[i,] - unique.particles[j,])) * unique.count[i] * unique.count[j]
    }
  }
  
  return(total.dist)
}


######
# Auxiliary function to calculate Metropolis ratio when log-likelihood values are very large

# largeRatio <- function(a,b){
#   # Determine the larger of the two numbers
#   
# }