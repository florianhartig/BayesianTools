#' Gelman Diagnostics
#' 
#' Runs Gelman Diagnotics for an object of class BayesianOutput
#' 
#' @author Florian Hartig
#' @param sampler an object of class mcmcSampler or mcmcSamplerList
#' @param thin parameter determining the thinning intervall. Either an integer or "auto" (default) for automatic thinning.
#' @param plot should a Gelman plot be generated
#' @param ... further arguments passed to \code{\link{getSample}}
#' 
#' @details The function calls [coda::gelman.diag] to calculate Gelman-Rubin diagnostics [coda::gelman.plot] to produce the plots. 
#' 
#' The idea of these diagnostics is to compare withing and between chain variance of several independent MCMC runs (Gelman & Rubin, 1992). The ratio of the 2 is called the potential scale reduction factor (psfr, also called Rhat). If psfr = 1, this suggest that the independent MCMC runs are essentially identical, and which in turn suggests that they have converged. In practice, values < 1.05, or sometimes < 1.1 for all parameters are considered acceptable. 
#' 
#' To obtain reliable Gelman-Rubin diagnostics, the independent MCMCs should be started at different points of the parameter space, ideally overdispersed.
#' 
#' The diagnostics also calculate a multivariate version of the psrf (mpsrf, Brooks & Gelman 1998). In practice, values < 1.1 or < 1.2 are often considered acceptable. While useful as an overview, mpsrf < 1.1 does not necessarily mean that all individual psrf < 1.05, and thus I would in doubt recommend looking at the individual psrf and decide on a case-by-case basis if a lack of convergence for a particular parameter is a concern. 
#' 
#' Also, note that convergence is a continuum, and different aspects of a posterior estimation converge with different speed. The rules about 1.05 were obtained by looking at the error of the posterior median / mean. If the goal for the inference is a posterior quantity that is more unstable than the mean, for example tail probabilities or the DIC, one should try to obtain large posterior samples with smaller psrf values. 
#' 
#' **Note on the use of Gelman diagnostics for population MCMCs, in particular the DE sampler family**: the Gelman diagnostics were originally designed for being applied to the outcome of several independent MCMC runs. Technically and practically, it can also be applied to a single population MCMC run that has several internal chains, such as DE, DEzs, DREAM, DREAMzs or T-Walk. As argued in ter Braak et al. (2008), the internal chains should be independent after burn-in. While this is likely correct, it also means that they are not completely independent before, and we observed this behavior in the use of the algorithms (i.e. that internal DEzs chains are more similar to each other than the chains of independent DEzs algorithms), see for example [BT issue 226](https://github.com/florianhartig/BayesianTools/issues/226). A concern is that this non-independence could lead to a failure to detect that the sampler hasn't converged yet, due to a wrong burn-in. We would therefore recommend to run several DEzs and check convergence with those, instead of running only one.  
#' 
#' @references
#' 
#' Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.
#' 
#' Brooks, SP. and Gelman, A. (1998) General methods for monitoring convergence of iterative simulations. Journal of Computational and Graphical Statistics, 7, 434-455.
#' 
#' ter Braak, Cajo JF, and Jasper A. Vrugt. "Differential evolution Markov chain with snooker updater and fewer chains." Statistics and Computing 18.4 (2008): 435-446.
#' 
#' @export
gelmanDiagnostics <- function(sampler, thin = "auto", plot = F, ...){
  sample = getSample(sampler, coda = T, ...)
  if(! ("mcmc.list" == class(sample))) stop("Trying to apply gelmanDiagnostics to an object that doesn't return an mcmc.list. Make sure you have a sampler that runs several chains, or an mcmcSamlerList")
  pars = ncol(sample[[1]])
  diag = NULL
  try({diag  = coda::gelman.diag(sample)}, silent = T)
  if(is.null(diag)){
    message("gelmanDiagnostics could not be calculated, possibly there is not enoug variance in your MCMC chains. Try running the sampler longer")
    diag = list()
    diag$psrf = matrix(nrow = pars, ncol = 2)
    rownames(diag$psrf) = colnames(sample)
    diag$mpsrf = NA  
  } 
  if(pars == 1) diag$mpsrf = NA  # fixes #221
  if(plot == T & ! is.na(diag$mpsrf)){
    # Wrapper around the gelman.plot to filter out getSample arguments from ...
    gP <- function(...,start, end, parametersOnly, coda, numSamples, whichParameters, reportDiagnostics, thin, plot, sampler) coda::gelman.plot(sample, ...)
    do.call(gP, as.list(match.call()))
  } 
  return(diag)
}

