#' Convenience function to create an object of class SMCSamplerList from a list of mcmc samplers
#' @param ... a list of MCMC samplers
#' @return a list of class smcSamplerList with each object being an smcSampler
#' @export
createSmcSamplerList <- function(...){
  smcList <- list(...)
  for (i in 1:length(smcList)){
    if (! ("mcmcSampler" %in% class(smcList[[i]])) ) stop("list objects are not of class mcmcSampler")
  }
  class(smcList) = c("smcSamplerList", "bayesianOutput")
  return(smcList)
}


#' @method summary smcSamplerList
#' @export
summary.smcSamplerList <- function(object, ...){
  sample = getSample(object, parametersOnly = T, ...)
  summary(sample)
}

#' @method print smcSamplerList
#' @export
print.smcSamplerList <- function(x, ...){
  print("smcSamplerList - you can use the following methods to summarize, plot or reduce this class:")
  print(methods(class ="smcSamplerList"))
}

#' @method plot smcSamplerList
#' @export
plot.smcSamplerList <- function(x, ...){
  marginalPlot(x, ...)
}

#' @export
getSample.smcSamplerList <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = 1,
                                     numSamples = NULL, whichParameters = NULL, reportDiagnostics = FALSE, ...){
  
  out = list()

  for (i in 1:length(sampler)){
    
    out[[i]] = getSample(sampler[[i]], parametersOnly = parametersOnly, whichParameters = whichParameters, start = start, end = end, thin = thin,
                         numSamples = numSamples, coda = F, reportDiagnostics = F)
      
  }
  out = combineChains(out, merge =F)
  
  return(out)
}





