#' calculates the Maxiumum APosteriori value (MAP)
#' @author Florian Hartig
#' @param bayesianOutput an object of class BayesianOutput (mcmcSampler, smcSampler, or mcmcList)
#' @param ... optional values to be passed on the the getSample function 
#' @details Currently, this function simply returns the parameter combination with the highest posterior in the chain. A more refined option would be to take the MCMC sample and do additional calculations, e.g. use an optimizer, a kernel density estimator, or some other tool to search / interpolate around the best value in the chain.
#' @seealso \code{\link{WAIC}}, \code{\link{DIC}}, \code{\link{marginalLikelihood}}
#' @export
MAP <- function(bayesianOutput, ...){
  
  samples = getSample(bayesianOutput, parametersOnly = F, ...)
  
  if("mcmcSamplerList" %in% class(bayesianOutput)) nPars <- bayesianOutput[[1]]$setup$numPars
  else nPars = bayesianOutput$setup$numPars
  
  best = which.max(samples[,nPars + 1])
  
  return(list(parametersMAP = samples[best, 1:nPars], valuesMAP = samples[best, (nPars + 1):(nPars + 3)] ))

}
  