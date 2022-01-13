#' Main wrapper function to start MCMCs, particle MCMCs and SMCs
#' @author Florian Hartig
#' @param bayesianSetup either a BayesianSetup (see \code{\link{createBayesianSetup}}) or a BayesianOutput created by runMCMC. The latter allows to continue a previous MCMC run. See details for how to restart a sampler. 
#' @param sampler sampling algorithm to be run. Default is DEzs. Options are "Metropolis", "AM", "DR", "DRAM", "DE", "DEzs", "DREAM", "DREAMzs", "SMC". For details see the help of the individual functions. 
#' @param settings list with settings for each sampler (see help of sampler for details). If a setting is not provided, defaults (see \code{\link{applySettingsDefault}}) will be used. You can see the default values by running \code{\link{applySettingsDefault}} with the respective sampler name, or in the help of the samplers. 
#' @details The runMCMC function can be started with either one of a) an object of class BayesianSetup with prior and likelihood function (recommended, see \code{\link{createBayesianSetup}}), b) a log posterior or other target function, or c) an object of class BayesianOutput created by runMCMC. The latter allows to continue a previous MCMC run. If a bayesianSetup is provided, check if appropriate parallization options are used - many samplers can make use of parallelization if this option is activated when the class is created.
#' 
#' For details about the different MCMC samplers, make sure you have read the Vignette (run vignette("BayesianTools", package="BayesianTools"). Also, see \code{\link{Metropolis}} for Metropolis based samplers, \code{\link{DE}} and \code{\link{DEzs}} for standard differential evolution samplers, \code{\link{DREAM}} and \code{\link{DREAMzs}} for DREAM sampler, \code{\link{Twalk}} for the Twalk sampler, and \code{\link{smcSampler}} for rejection and Sequential Monte Carlo sampling.\cr
#' 
#' The samplers "AM", "DR", and "DRAM" are special cases of the "Metropolis" sampler and are shortcuts for predefined settings ("AM": adapt=TRUE; "DR": DRlevels=2; "DRAM": adapt=True, DRlevels=2).
#' 
#' The settings list allows to change the settings for the MCMC samplers and some other options. For the MCMC sampler settings, see their help files. Global options that apply for all MCMC samplers are: iterations (number of MCMC iterations), and nrChains (number of chains to run). Note that running several chains is not done in parallel, so if time is an issue it will be better to run the MCMCs individually and then combine them via \code{\link{createMcmcSamplerList}} into one joint object. 
#' 
#' 
#' Startvalues: all samplers allow to provide explicit startvalues. If startvalues are not provided, they are sampled from the prior. Usually, this is a good choice, so don't feel compelled to provide startvalues. 
#' 
#' Note that DE and DREAM variants as well as SMC and T-walk require a population to start, which should be provided as a matrix. Default (NULL) sets the population size for DE to 3 x dimensions of parameters, for DREAM to 2 x dimensions of parameters and for DEzs and DREAMzs to three, sampled from the prior. Note also that the zs variants of DE and DREAM require two populations, the current population and the z matrix (a kind of memory) - if you want to set both, provide a list with startvalue$X and startvalue$Z. 
#' 
#' setting startValue for sampling with nrChains > 1 : if you want to provide different start values for the different chains, provide them as a list
#' 
#' @return The function returns an object of class mcmcSampler (if one chain is run) or mcmcSamplerList. Both have the superclass bayesianOutput. It is possible to extract the samples as a coda object or matrix with \code{\link{getSample}}. 
#' It is also possible to summarize the posterior as a new prior via \code{\link{createPriorDensity}}.
#' @example /inst/examples/mcmcRun.R
#' @seealso \code{\link{createBayesianSetup}} 
#' @export 
runMCMC <- function(bayesianSetup , sampler = "DEzs", settings = NULL){
  
  options(warn = 0)
  
  ptm <- proc.time()  
 
  ####### RESTART ########## 
  
  if("bayesianOutput" %in% class(bayesianSetup)){
    
    # TODO - the next statements should have assertions in case someone overwrites the existing setting or similar
    
    previousMcmcSampler <- bayesianSetup
    
    
    # Catch the settings in case of nrChains > 1
    if(!("mcmcSamplerList" %in% class(previousMcmcSampler) |  "smcSamplerList" %in% class(previousMcmcSampler) )){
      if(is.null(settings)) settings <- previousMcmcSampler$settings
      setup <- previousMcmcSampler$setup
      sampler <- previousMcmcSampler$settings$sampler
      previousSettings <- previousMcmcSampler$settings
    } else{ 
      if(is.null(settings)) settings <- previousMcmcSampler[[1]]$settings
      settings$nrChains <- length(previousMcmcSampler)
      setup <- previousMcmcSampler[[1]]$setup
      sampler <- previousMcmcSampler[[1]]$settings$sampler
      previousSettings <- previousMcmcSampler[[1]]$settings
    }
    
    # Set settings$sampler (only needed if new settings are supplied)
    settings$sampler <- sampler
    
    # overwrite new settings
    for(name in names(settings)) previousSettings[[name]] <- settings[[name]]
    
    settings <- previousSettings
    
    # Check if previous settings will be new default
    
    previousMcmcSampler$settings <- applySettingsDefault(settings = settings, sampler = settings$sampler, check = TRUE)
    
    restart <- TRUE

    
  ## NOT RESTART STARTS HERE ###################
    
  }else if(class(bayesianSetup) == "BayesianSetup"){
    restart <- FALSE

    if(is.null(settings$parallel)) settings$parallel <- bayesianSetup$parallel
    if(is.numeric(settings$parallel)) settings$parallel <- TRUE      

    setup <- checkBayesianSetup(bayesianSetup, parallel = settings$parallel)
    settings <- applySettingsDefault(settings = settings, sampler = sampler, check = TRUE)
  } else stop("runMCMC requires a class of type BayesianOutput or BayesianSetup")
  
  ###### END RESTART ##############
  
  
  # TODO - the following statement should be removed once all further functions access settings$sampler instead of sampler
  # At the moment only the same sampler can be used to restart sampling.
  sampler = settings$sampler
  
  #### Assertions
  if(!restart && setup$numPars == 1) if(!getPossibleSamplerTypes()$univariate[which(getPossibleSamplerTypes()$BTname == settings$sampler)]) stop("This sampler can not be applied to a univariate distribution")

  if(restart == T) if(!getPossibleSamplerTypes()$restartable[which(getPossibleSamplerTypes()$BTname == settings$sampler)]) stop("This sampler can not be restarted")
  
  ########### Recursive call in case multiple chains are to be run
  if(settings$nrChains >1){

    # Initialize output list
    out<- list()
    
    # Run several samplers
    for(i in 1:settings$nrChains){
      
      settingsTemp <- settings
      settingsTemp$nrChains <- 1 # avoid infinite loop
      settingsTemp$currentChain <- i
      
      if(restart){
        out[[i]] <- runMCMC(bayesianSetup = previousMcmcSampler[[i]], settings = settingsTemp)
      }else{
        if(is.list(settings$startValue)) settingsTemp$startValue = settings$startValue[[i]]
        out[[i]] <- runMCMC(bayesianSetup = setup, sampler = settings$sampler, settings = settingsTemp)
      }
    }
    if(settings$sampler == "SMC") class(out) = c("smcSamplerList", "bayesianOutput")
    else class(out) = c("mcmcSamplerList", "bayesianOutput")
    return(out)
    
  ######### END RECURSIVE CALL 
  # MAIN RUN FUNCTION HERE  
  }else{
    
    # check start values 
    setup$prior$checkStart(settings$startValue)
    
    
    if (sampler == "Metropolis" || sampler == "AM" || sampler == "DR" || sampler == "DRAM"){
      if(restart == FALSE){
        mcmcSampler <- Metropolis(bayesianSetup = setup, settings = settings)
        mcmcSampler <- sampleMetropolis(mcmcSampler = mcmcSampler, iterations = settings$iterations)
      } else {
        mcmcSampler <- sampleMetropolis(mcmcSampler = previousMcmcSampler, iterations = settings$iterations) 
      }
    }
    
    

    ############## Differential Evolution #####################
    if (sampler == "DE"){
      
      if(restart == F) out <- DE(bayesianSetup = setup, settings = settings)
      else out <- DE(bayesianSetup = previousMcmcSampler, settings = settings)
      
      #out <- DE(bayesianSetup = bayesianSetup, settings = list(startValue = NULL, iterations = settings$iterations, burnin = settings$burnin, eps = settings$eps, parallel = settings$parallel, consoleUpdates = settings$consoleUpdates,
       #            blockUpdate = settings$blockUpdate, currentChain = settings$currentChain))
      
      mcmcSampler = list(
        setup = setup,
        settings = settings,
        chain = out$Draws,
        X = out$X,
        sampler = "DE"
      )
    }
    
    ############## Differential Evolution with snooker update
    if (sampler == "DEzs"){
      # check z matrix
      if(!is.null(settings$Z)) setup$prior$checkStart(settings$Z,z = TRUE)
      
      if(restart == F) out <- DEzs(bayesianSetup = setup, settings = settings)
      else out <- DEzs(bayesianSetup = previousMcmcSampler, settings = settings)
      
      mcmcSampler = list(
        setup = setup,
        settings = settings,
        chain = out$Draws,
        X = out$X,
        Z = out$Z,
        sampler = "DEzs"
      )
    }
    
    ############## DREAM   
    if (sampler == "DREAM"){
      
      if(restart == F) out <- DREAM(bayesianSetup = setup, settings = settings)
      else out <- DREAM(bayesianSetup = previousMcmcSampler, settings = settings)

      mcmcSampler = list(
        setup = setup,
        settings = settings,
        chain = out$chains,
        pCR = out$pCR,
        sampler = "DREAM",
        lCR = out$lCR,
        X = out$X,
        delta = out$delta
      )
    }
    
    ############## DREAMzs   
    if (sampler == "DREAMzs"){
        # check z matrix
        if(!is.null(settings$Z)) setup$prior$checkStart(settings$Z,z = TRUE)
      
        if(restart == F) out <- DREAMzs(bayesianSetup = setup, settings = settings)
        else out <- DREAMzs(bayesianSetup = previousMcmcSampler, settings = settings)
      
        mcmcSampler = list(
        setup = setup,
        settings = settings,
        chain = out$chains,
        pCR = out$pCR,
        sampler = "DREAMzs",
        JumpRates = out$JumpRates,
        X = out$X,
        Z = out$Z
        )
      
    }
    
    if(sampler == "Twalk"){
      warning("At the moment using T-walk is discouraged: numeric instability")
      if(!restart){
        if(is.null(settings$startValue)){
          settings$startValue = bayesianSetup$prior$sampler(2)
        }
        mcmcSampler <- Twalk(bayesianSetup = setup, settings = settings)
      }else{
        mcmcSampler <- Twalk(bayesianSetup = previousMcmcSampler, settings = settings)
      }
      mcmcSampler$setup <- setup
      mcmcSampler$sampler <- "Twalk"
    }
    
    
    if ((sampler != "SMC")){
      class(mcmcSampler) <- c("mcmcSampler", "bayesianOutput")
    }
    
    ############# SMC #####################
    
    if (sampler == "SMC"){

      mcmcSampler <- smcSampler(bayesianSetup = bayesianSetup, initialParticles = settings$initialParticles, iterations = settings$iterations, resampling = settings$resampling, resamplingSteps = settings$resamplingSteps, proposal = settings$proposal, adaptive = settings$adaptive, proposalScale = settings$proposalScale )
      mcmcSampler$settings = settings
    }

    mcmcSampler$settings$runtime = mcmcSampler$settings$runtime + proc.time() - ptm
    if(is.null(settings$message) || settings$message == TRUE){
      message("runMCMC terminated after ", mcmcSampler$settings$runtime[3], "seconds")
    }
    return(mcmcSampler)
  }
}


#bayesianSetup = bayesianSetup, initialParticles = settings$initialParticles, iterations = settings$iterations, resampling = settings$resampling, resamplingSteps = settings$resamplingSteps, proposal = settings$proposal, adaptive = settings$adaptive, parallel = settings$parallel


#' Provides the default settings for the different samplers in runMCMC
#' @author Florian Hartig
#' @param settings optional list with parameters that will be used instead of the defaults
#' @param sampler one of the samplers in \code{\link{runMCMC}} 
#' @param check logical determines whether parameters should be checked for consistency
#' @details see \code{\link{runMCMC}} 
#' @export
applySettingsDefault<-function(settings=NULL, sampler = "DEzs", check = FALSE){
  
  if(is.null(settings)) settings = list()
  
  if(!is.null(sampler)){
    if(!is.null(settings$sampler)) {
      # TODO: this is a bit hacky. The best would prabably be to change the Metropolis function to allow AM, DR and DRAM
      #       arguments and call applySettingsDefault for those
      if (settings$sampler %in% c("AM", "DR", "DRAM") && sampler == "Metropolis") {
        sampler <- settings$sampler
      }
      if(settings$sampler != sampler) {
        warning("sampler argument overwrites an existing settings$sampler in applySettingsDefault. This only makes sense if one wants to take over settings from one sampler to another")
      }
    }
    settings$sampler = sampler
  }
  
  if(!settings$sampler %in% getPossibleSamplerTypes()$BTname) stop("trying to set values for a sampler that does not exist")
  
  if (settings$sampler == "AM") {
    defaultSettings <- getMetropolisDefaultSettings()
    defaultSettings$adapt <- TRUE # RB: is already T
  }
  
  if (settings$sampler == "DR") {
    defaultSettings <- getMetropolisDefaultSettings()
    defaultSettings$DRlevels <- 2
  }
  
  if (settings$sampler == "DRAM") {
    defaultSettings <- getMetropolisDefaultSettings()
    defaultSettings$adapt <- TRUE # RB: is already T
    defaultSettings$DRlevels <- 2
  }
  
  if (settings$sampler == "Metropolis"){
    defaultSettings <- getMetropolisDefaultSettings()
  }
  
  if (settings$sampler == "DE"){
    defaultSettings <- getDEDefaultSettings()
    defaultSettings$f <- -2.38
    defaultSettings$blockUpdate <- list("none", k = NULL, h = NULL, pSel = NULL, pGroup = NULL, 
                                       groupStart = 1000, groupIntervall = 1000)
  }
  
  if (settings$sampler == "DEzs"){
    defaultSettings <- getDEDefaultSettings()
    defaultSettings$f <- 2.38
    defaultSettings$blockUpdate <- list("none", k = NULL, h = NULL, pSel = NULL, pGroup = NULL, 
                                        groupStart = 1000, groupIntervall = 1000)
    defaultSettings <- c(defaultSettings, list(Z = NULL,
                                               zUpdateFrequency = 1,
                                               pSnooker = 0.1,
                                               pGamma1 = 0.1,
                                               eps.mult =0.2,
                                               eps.add = 0))
  }
  

  if (settings$sampler == "DREAM"){
    defaultSettings <- getDREAMDefaultSettings()
    defaultSettings$pCRupdate <- TRUE
  }
    
  if (settings$sampler == "DREAMzs"){
    defaultSettings <- getDREAMDefaultSettings()
    defaultSettings$pCRupdate <- FALSE
    defaultSettings$Z = NULL
    defaultSettings$ZupdateFrequency = 10
    defaultSettings$pSnooker = 0.1
  }
  
  if (settings$sampler == "SMC"){
    defaultSettings = list(iterations = 10, 
                           resampling = T, 
                           resamplingSteps = 2, 
                           proposal = NULL, 
                           adaptive = T, 
                           proposalScale = 0.5, 
                           initialParticles = 1000
                           )
  }
  
  
  if (settings$sampler == "Twalk"){
    defaultSettings = list(iterations = 10000, 
                           consoleUpdates=100, 
                           at = 6, 
                           aw = 1.5, 
                           pn1 = NULL, 
                           Ptrav = 0.4918, 
                           Pwalk = NULL, 
                           Pblow = NULL, 
                           startValue = NULL,
                           burnin = 0,
                           thin = 1,
                           message = TRUE)
  }
  
  
  ## CHECK DEFAULTS

  if(check){
    nam = c(names(defaultSettings), "sampler", "nrChains",
            "runtime", "sessionInfo", "parallel")
    
    ind <- which((names(settings) %in% nam == FALSE))
    
    nam_n <- names(settings)[ind]
    for(i in 1:length(nam_n)) nam_n[i] <- paste(nam_n[i], " ")
    
    if(length(ind) > 0){
      message("Parameter(s) ", nam_n , " not used in ", settings$sampler, "\n")
    }
  }
  
  defaultSettings$nrChains = 1
  defaultSettings$runtime = 0
  defaultSettings$sessionInfo = utils::sessionInfo()
  
  nam = names(defaultSettings)
  
  for (i in 1:length(defaultSettings)){
    if(! nam[i] %in% names(settings)){
      addition = list( defaultSettings[[i]])
      names(addition) = nam[i]
      settings = c(settings, addition)
    }
  }    

  
  if (! is.null(settings$burnin)){
    if (settings$burnin > settings$iterations) stop("BayesianToools::applySettingsDefault - setting burnin cannnot be larger than setting iteration")
    if (! is.null(settings$adaptationNotBefore)){
      if (settings$burnin >= settings$adaptationNotBefore) stop("BayesianToools::applySettingsDefault - setting burnin cannnot be larger than setting adaptationNotBefore") 
    }    
  }

  return(settings)  
}


#' Help function to find starvalues and proposalGenerator settings
#' @author Florian Hartig
#' @param proposalGenerator proposal generator
#' @param bayesianSetup either an object of class bayesianSetup created by \code{\link{createBayesianSetup}} (recommended), or a log target function 
#' @param settings list with settings
#' @keywords internal
setupStartProposal <- function(proposalGenerator = NULL, bayesianSetup, settings){
  
  # Proposal
  range = (bayesianSetup$prior$upper - bayesianSetup$prior$lower) / 50
  
  if(is.null(settings$startValue)) settings$startValue = (bayesianSetup$prior$upper + bayesianSetup$prior$lower) / 2
  
  if (length(range) != bayesianSetup$numPars) range = rep(1,bayesianSetup$numPars)
  
  if(is.null(proposalGenerator)){
    proposalGenerator = createProposalGenerator(range, gibbsProbabilities = settings$gibbsProbabilities)
  }

  ####### OPTIMIZATION
  
  if (settings$optimize == T){
    if(is.null(settings$message) || settings$message == TRUE){
    cat("BT runMCMC: trying to find optimal start and covariance values", "\b")
    }
    
    target <- function(x){
      out <- bayesianSetup$posterior$density(x)
      if (out == -Inf) out =  -1e20 # rnorm(1, mean = -1e20, sd = 1e-20)
      return(out)
    }
    
    try( {
      if(bayesianSetup$numPars > 1) optresul <- optim(par=settings$startValue,fn=target, method="Nelder-Mead", hessian=F, control=list("fnscale"=-1, "maxit" = 10000))
      else optresul <- optim(par=settings$startValue,fn=target, method="Brent", hessian=F, control=list("fnscale"=-1, "maxit" = 10000), lower = bayesianSetup$prior$lower, upper = bayesianSetup$prior$upper)      
      settings$startValue = optresul$par
      hessian = numDeriv::hessian(target, optresul$par)
      
     
      proposalGenerator$covariance = as.matrix(Matrix::nearPD(MASS::ginv(-hessian))$mat)
      #proposalGenerator$covariance = MASS::ginv(-optresul$hessian)
      
      # Create objects for startValues and covariance to add space between values
      startV <-covV <- character()
      
      for(i in 1:length(settings$startValue)){
        startV[i] <- paste(settings$startValue[i], "")
      } 
      for(i in 1:length( proposalGenerator$covariance)){
        covV[i] <- paste( proposalGenerator$covariance[i], "")
      } 
      
      if(is.null(settings$message) || settings$message == TRUE){
      message("BT runMCMC: Optimization finished, setting startValues to " , 
              startV, " - Setting covariance to " , covV)
      }
      
      proposalGenerator = updateProposalGenerator(proposalGenerator)
      
    }
    , silent = FALSE)
  }  
  out = list(proposalGenerator = proposalGenerator, settings = settings)
  return(out)
}

#' Returns possible sampler types
#' @export
#' @author Florian Hartig
getPossibleSamplerTypes <- function(){
  
  out = list(
    BTname = c("AM", "DR", "DRAM", "Metropolis", "DE", "DEzs", "DREAM", "DREAMzs", "Twalk", "SMC"),
    possibleSettings = list() ,
    possibleSettingsName = list() ,
    
    univariatePossible = c(T, T, T, T, T, T, T, T, T, F),
    restartable = c(T, T, T, T, T, T, T, T, T, F)
  )

  return(out)
} 

#' Returns Metropolis default settings
#' @author Tankred Ott
#' @keywords internal
getMetropolisDefaultSettings <- function () {
  defaultSettings <- getDefaultSettings()
  defaultSettings <- c(defaultSettings, list(optimize = T,
                                             proposalGenerator = NULL,
                                             adapt = T,
                                             adaptationInterval = 500,
                                             adaptationNotBefore = 3000,
                                             DRlevels = 1 ,
                                             proposalScaling = NULL,
                                             adaptationDepth = NULL,
                                             temperingFunction = NULL,
                                             proposalGenerator = NULL,
                                             gibbsProbabilities = NULL))
  return(defaultSettings)
}

#' Returns DREAM default settings
#' @author Robert Bosek
#' @keywords internal
getDREAMDefaultSettings <- function () {
  defaultSettings <- getDefaultSettings()
  defaultSettings$consoleUpdates <- 10
  defaultSettings$eps <- 0
  defaultSettings <- c(defaultSettings, list(nCR = 3,
                                             gamma = NULL,
                                             e = 5e-2,  # TODO check
                                             updateInterval = 10,
                                             adaptation = 0.2,
                                             DEpairs = 2))
  
  return(defaultSettings)
}


#' Returns Differential-Evolution default settings
#' @author Robert Bosek
#' @keywords internal
getDEDefaultSettings <- function () {
  defaultSettings <- getDefaultSettings()
  defaultSettings$eps <- 0
  return(defaultSettings)
}

#' Returns basic default settings
#' @author Robert Bosek
#' @keywords internal
getDefaultSettings <- function () {
  defaultSettings <- list(startValue = NULL,
                         iterations = 10000,
                         burnin = 0,
                         thin = 1,
                         consoleUpdates = 100,
                         currentChain = 1,
                         parallel = NULL,
                         message = TRUE)
  return(defaultSettings)
}


