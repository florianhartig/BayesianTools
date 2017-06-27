# NOTE: The functions in this class are just templates that are to be implemented for all subclasses of BayesianOutput. They are not functional. 


#' Extracts the sample from a bayesianOutput
#' @author Florian Hartig
#' @param sampler an object of class mcmcSampler, mcmcSamplerList, smcSampler, smcSamplerList
#' @param parametersOnly if F, likelihood, posterior and prior values are also provided in the output
#' @param coda works only for mcmc classes - provides output as a coda object. Note: if mcmcSamplerList contains mcmc samplers such as DE that have several chains, the internal chains will be collapsed. This may not be the desired behavior for all applications. 
#' @param start for mcmc samplers start value in the chain. For SMC samplers, start particle
#' @param end for mcmc samplers end value in the chain. For SMC samplers, end particle
#' @param thin thinning parameter. Either an integer determining the thinning intervall (default is 1) or "auto" for automatic thinning.
#' @param numSamples sample size (only used if thin = 1)
#' @param whichParameters possibility to select parameters by index
#' @param reportDiagnostics logical, determines whether settings should be included in the output
#' @param ... further arguments
#' @example /inst/examples/getSampleHelp.R
#' @details If thin is greater than the total number of samples in the sampler object the first and the last element (of each chain if a sampler with multiples chains is used) are sampled. If numSamples is greater than the total number of samples all samples are selected. In both cases a warning is displayed.
#' @export
getSample <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = 1, numSamples = NULL, whichParameters = NULL, includesProbabilities = F, reportDiagnostics = FALSE, ...) UseMethod("getSample")


getSample.matrix <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", whichParameters = NULL, includesProbabilities = F, reportDiagnostics = F, ...){
    if(is.null(end)) end = nrow(sampler)

    if(includesProbabilities) nPars = ncol(sampler) - 3 else nPars = ncol(sampler)
    
    if(parametersOnly == T | includesProbabilities == F) {
      out = sampler[start:end,1:nPars] 
      if(class(out) == "numeric") out = as.matrix(sampler) # case 1 parameter
    } else {
      out = out[start:end,] 
      #if(!is.null(sampler$setup$names)) colnames(out) = c(sampler$setup$names, "Lposterior", "Llikelihood", "Lprior")
    }
    
    ########################
    # THINNING
    if (thin == "auto"){
      thin = max(floor(nrow(out) / 5000), 1)
    }
    if(is.null(thin) || thin == F || thin < 1 || is.nan(thin)) thin = 1
    if (thin > nrow(sampler)) warning("thin is greater than the total number of samples!")
    if (! thin == 1){
      sel = seq(1,dim(out)[1], by = thin )
      out = out[sel,]
    }
    #############
    
    if (!is.matrix(out)) out <- matrix(out, nrow = 1)
        
    if (!is.null(whichParameters)) out = out[,whichParameters]
    if(coda == T) out = makeObjectClassCodaMCMC(out, start = start, end = end, thin = thin)

    if(reportDiagnostics == T){
      return(list(chain = out, start = start, end = end, thin = thin))
    } else return(out)
}


#' @author Tankred Ott
#' @export
# TODO: This is right now only a helper function for getSample.mcmc. It is needed to return a vector istead of a matrix, if 
#       the mcmc object passed to getSample.mcmc contains a vector.
getSample.double <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", whichParameters = NULL, includesProbabilities = F, reportDiagnostics = F, ...){
  if(is.null(end)) end = length(sampler)
  
  nTotalSamples <- length(sampler)
  
  thin = correctThin(nTotalSamples, thin)
  
  sel = seq(1, nTotalSamples, by = thin)
  return(sampler[sel])
}


#' @author Tankred Ott
#' @export
# TODO: This is right now only a helper function for getSample.mcmc. It is needed to return a vector istead of a matrix, if 
#       the mcmc object passed to getSample.mcmc contains a vector.
getSample.integer <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", whichParameters = NULL, includesProbabilities = F, reportDiagnostics = F, ...){
  if(is.null(end)) end = length(sampler)
  
  nTotalSamples <- length(sampler)
  
  thin = correctThin(nTotalSamples, thin)
  
  sel = seq(1, nTotalSamples, by = thin)
  return(sampler[sel])
}


getSample.data.frame <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", whichParameters = NULL, includesProbabilities = F, reportDiagnostics = F, ...){
  getSample(matrix(sampler), parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, whichParameters = whichParameters, includesProbabilities = includesProbabilities, reportDiagnostics = reportDiagnostics)
}


# The following two S3 implementations make getSample compatible with coda::mcmc and coda::mcmc.list

#' @author Tankred Ott
#' @export
getSample.mcmc <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", whichParameters = NULL, includesProbabilities = F, reportDiagnostics = F, ...){

  # TODO: implement handling of wrong inputs?
  
  if(coda == T){
    # mcmc objects can contain matrices or vectors
    if (is.matrix(sampler)) {
      nTotalSamples <- nrow(sampler)
    } else {
      nTotalSamples <- length(sampler)
    }
    
    if (is.null(end)) end = nTotalSamples
    
    # check/correct thin
    thin <- correctThin(nTotalSamples, thin)
    
    # see http://svitsrv25.epfl.ch/R-doc/library/coda/html/window.mcmc.html
    # for coda's window implementation
    return(window(sampler, start = start, end = end, thin = thin))
    
  } else if(coda == F){
    # mcmc objects can contain matrices or vectors
    if (is.matrix(sampler)) {
      out <- getSample(as.matrix(sampler), parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, whichParameters = whichParameters, includesProbabilities = includesProbabilities, reportDiagnostics = reportDiagnostics)
    } else {
      out <- getSample(as.vector(sampler), parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, whichParameters = whichParameters, includesProbabilities = includesProbabilities, reportDiagnostics = reportDiagnostics)
    }
    return(out)
  }
}


#' @author Tankred Ott
#' @export
getSample.mcmc.list <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", whichParameters = NULL, includesProbabilities = F, reportDiagnostics = F, ...){
  
  # TODO: implement handling of wrong inputs?
  
  if(coda == T){
    
    if (is.matrix(sampler[[1]])) {
      nTotalSamples <- nrow(sampler[[1]])
    } else {
      nTotalSamples <- length(sampler[[1]])
    }
    
    if (is.null(end)) end = nTotalSamples
    
    # check/correct thin
    thin <- correctThin(nTotalSamples, thin)
    
    # see http://svitsrv25.epfl.ch/R-doc/library/coda/html/window.mcmc.html
    # for coda's window implementation
    return(window(sampler, start = start, end = end, thin = thin))
    
  } else if(coda == F){
    if(is.matrix(sampler[[1]])) {
      return(getSample(combineChains(sampler), parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, whichParameters = whichParameters, includesProbabilities = includesProbabilities, reportDiagnostics = reportDiagnostics))
    } else {
      return(as.vector(getSample(combineChains(sampler), parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, whichParameters = whichParameters, includesProbabilities = includesProbabilities, reportDiagnostics = reportDiagnostics)))
    }
  }
}
  

