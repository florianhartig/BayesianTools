#' Runs Gelman Diagnotics over an BayesianOutput
#' @author Florian Hartig
#' @param sampler an object of class mcmcSampler or mcmcSamplerList
#' @param thin parameter determining the thinning intervall. Either an integer or "auto" (default) for automatic thinning.
#' @param plot should a Gelman plot be generated
#' @param ... further arguments passed to \code{\link{getSample}}
#' 
#' @details The function calls the coda package to calculate Gelman diagnostics and plots
#' 
#' The original idea is that this function is applied to the outcome of several independent MCMC runs. Technically and practically, it can also be applied to a single MCMC run that has several internal chains, such as DE, DEzs, DREAM, DREAMzs or T-Walk. As argued in ter Braak et al. (2008), the internal chains should be independent after burn-in. While this is likely correct, it also means that they are not completely independent before, and we observed this behavior in the use of the algorithms (i.e. that internal DEzs chains are more similar to each other than the chains of independent DEzs algorithms). A concern is that this non-independence could lead to a failure to detect that the sampler hasn't converged yet. We would therefore recommend to run several DEzs and check convergence with those, instead of running only one.  
#' 
#' ter Braak, Cajo JF, and Jasper A. Vrugt. "Differential evolution Markov chain with snooker updater and fewer chains." Statistics and Computing 18.4 (2008): 435-446.
#' 
#' 
#' @export
gelmanDiagnostics <- function(sampler, thin = "auto", plot = F, ...){
  sample = getSample(sampler, coda = T, ...)
  if(! ("mcmc.list" == class(sample))) stop("Trying to apply gelmanDiagnostics to an object that doesn't return an mcmc.list. Make sure you have a sampler that runs several chains, or an mcmcSamlerList")
  diag = coda::gelman.diag(sample)
  gP <- function(...,start, end, parametersOnly, coda, numSamples, whichParameters, includesProbabilities, reportDiagnostics, thin, plot, sampler) coda::gelman.plot(sample, ...)
  if(plot == T) do.call(gP, as.list(match.call()))
  return(diag)
}


