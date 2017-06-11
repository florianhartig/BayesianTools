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
getSample <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = 1, numSamples = NULL, whichParameters = NULL, reportDiagnostics = FALSE, ...) UseMethod("getSample")


getSample.matrix <- function(mat, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", whichParameters = NULL, includesProbabilities = F, reportDiagnostics = F, ...){

    if(is.null(end)) end = nrow(mat)

    if(includesProbabilities) nPars = ncol(mat) - 3 else nPars = ncol(mat)
    
    if(parametersOnly == T | includesProbabilities == F) {
      out = mat[start:end,1:nPars] 
      if(class(out) == "numeric") out = as.matrix(out) # case 1 parameter
    } else {
      out = out[start:end,] 
      #if(!is.null(sampler$setup$names)) colnames(out) = c(sampler$setup$names, "Lposterior", "Llikelihood", "Lprior")
    }
    
    ########################
    # THINNING
    if (thin == "auto"){
      thin = max(floor(nrow(out) / 5000),1)
    }
    if(is.null(thin) || thin == F || thin < 1) thin = 1
    if (thin > nrow(mat)) warning("thin is greater than the total number of samples!")
    if (! thin == 1){
      sel = seq(1,dim(out)[1], by = thin )
      out = out[sel,]
    }
    #############
    
    if (!is.null(whichParameters)) out = out[,whichParameters]
    if(coda == T) out = makeObjectClassCodaMCMC(out, start = start, end = end, thin = thin)

    if(reportDiagnostics == T){
      return(list(chain = out, start = start, end = end, thin = thin))
    } else return(out)
}


getSample.data.frame <- function(data, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", whichParameters = NULL, includesProbabilities = F, reportDiagnostics = F, ...){
  getSample(matrix(data), parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, whichParameters = whichParameters, includesProbabilities = includesProbabilities, reportDiagnostics = reportDiagnostics)
}

# The following two S3 implementations make getSample compatible with coda::mcmc and coda::mcmc.list


getSample.mcmc <- function(data, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", whichParameters = NULL, includesProbabilities = F, reportDiagnostics = F, ...){
  
  if(coda == T){
    
    # coda objects can contain matrices or vectors
    if (is.matrix(data)) {
      nTotalSamples <- nrow(data)
    } else {
      nTotalSamples <- length(data)
    }
    
    thin <- correctThin(nTotalSamples, thin)
    
    if (is.null(end)) end = nTotalSamples
    
    # see http://svitsrv25.epfl.ch/R-doc/library/coda/html/window.mcmc.html
    # for coda's window implementation
    return(window(data, start = start, end = end, thin = thin))
    
  } else if(coda == F){
    return(getSample(as.matrix(data), parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, whichParameters = whichParameters, includesProbabilities = includesProbabilities, reportDiagnostics = reportDiagnostics))
  }
}

getSample.mcmc.list <- function(data, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", whichParameters = NULL, includesProbabilities = F, reportDiagnostics = F, ...){
  
  if(coda == T){
    
    # TODO - if the object is already coda, don't convert, but just select, potentially using http://svitsrv25.epfl.ch/R-doc/library/coda/html/window.mcmc.html
    
  } else if(coda == F){
    
    # TODO - this doesn't work properly, the chains are either ordered one after the other  (as.matrix), or with as.array as 3-dim array, but they should be shuffled according to BT Logic
    
    # see http://svitsrv25.epfl.ch/R-doc/library/coda/html/mcmc.convert.html
    
    getSample(as.array(data), parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, whichParameters = whichParameters, includesProbabilities = includesProbabilities, reportDiagnostics = reportDiagnostics)
  }
}
  

