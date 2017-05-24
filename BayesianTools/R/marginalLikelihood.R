
# Motivation for this functions from 
# https://radfordneal.wordpress.com/2008/08/17/the-harmonic-mean-of-the-likelihood-worst-monte-carlo-method-ever/
# https://gist.github.com/gaberoo/4619102


#  ' @export
#marginalLikelihood <- function(x,lik,V,sampler$setup$likelihood$density,sampler$setup$prior$density,..., num.samples=1000,log=TRUE) UseMethod("marginalLikelihood")

#' Calcluated the marginal likelihood from a set of MCMC samples
#' @export
#' @author Florian Hartig
#' @param sampler an object that implements the getSample function, i.e. a mcmc / smc Sampler (list)
#' @param numSamples number of samples to use. How this works, and if it requires recalculating the likelihood, depends on the method
#' @param method method to choose. Currently available are "Chib" (default), the harmonic mean "HM", and sampling from the prior "prior". See details
#' @param ... further arguments passed to \code{\link{getSample}}
#' @details The function currently implements three ways to calculate the marginal likelihood.\cr
#'  The recommended way is the method "Chib" (Chib and Jeliazkov, 2001). which is based on MCMC samples, but performs additional calculations. 
#'  Despite being the current recommendation, note there are some numeric issues with this algorithm that may limit reliability for larger dimensions.
#'   
#'  The harmonic mean approximation,
#'   is implemented only for comparison. Note that the method is numerically 
#'   unrealiable and usually should not be used. \cr
#' 
#' The third method is simply sampling from the prior. While in principle unbiased,
#'  it will only converge for a large number of samples, and is therefore
#'   numerically inefficient. \cr
#'    
#' @example /inst/examples/marginalLikelihoodHelp.R
#' @references Chib, Siddhartha, and Ivan Jeliazkov. "Marginal likelihood from the Metropolis-Hastings output." Journal of the American Statistical Association 96.453 (2001): 270-281.
#' @seealso \code{\link{WAIC}}, \code{\link{DIC}}, \code{\link{MAP}}
marginalLikelihood <- function(sampler, numSamples = 1000, method = "Chib", ...){
  
  if (method == "Chib"){
    
    chain = getSample(sampler = sampler, parametersOnly = F, ...)
    
    if(class(sampler)[1] %in% c("mcmcSamplerList", "smcSamplerList")) sampler = sampler[[1]]
    
    x = chain[,1:sampler$setup$numPars]
    
    lik = chain[,sampler$setup$numPars + 2]
    MAPindex = which.max(chain[,sampler$setup$numPars + 1])
    
    #propGen = createProposalGenerator(covariance = cov(x))
    
    V = cov(x)
    
    # calculate reference parameter 
    
    theta.star <- x[MAPindex,]
    lik.star <- lik[MAPindex]
    
    # get samples from posterior
    
    g <- sample.int(nrow(x),numSamples,replace=TRUE) # should replace really be true?
    q.g <- mvtnorm::dmvnorm(x[g,],mean=theta.star,sigma=V,log=FALSE)
    lik.g <- lik[g]
    alpha.g <- sapply(lik.g,function(l) min(1,exp(lik.star-l))) # Metropolis Ratio
    
    #lik.g <- apply(theta.g,1,sampler$setup$likelihood$density,...)
    
    
    # get samples from proposal
    theta.j <-mvtnorm::rmvnorm(numSamples,mean=theta.star,sigma=V)
    lik.j <- apply(theta.j,1,sampler$setup$likelihood$density)
    alpha.j <- sapply(lik.j,function(l) min(1,exp(l-lik.star)))  # Metropolis Ratio
    
    # Prior 
    pi.hat <- mean(alpha.g*q.g)/mean(alpha.j)
    pi.star <- 0
    
    if (!is.null(sampler$setup$prior$density)) pi.star <- sampler$setup$prior$density(theta.star)
    ln.m <- lik.star + pi.star - log(pi.hat)
    
    out = list(marginalLikelihod=ln.m,ln.lik.star=lik.star,ln.pi.star=pi.star,ln.pi.hat=log(pi.hat), method = "Chib")
    
  } else if (method == "HM"){
    
    chain = getSample(sampler = sampler, parametersOnly = F, ...)
    lik = chain[,sampler$setup$numPars + 2]
    ml = log(1/mean(1/exp(lik)))
    # ml = 1 / logSumExp(-lik, mean = T) function needs to be adjusted
    out <- list(marginalLikelihod=ml, method ="HM")
    
  } else if (method == "Prior"){
    samples = sampler$setup$prior$sampler(numSamples)
    likelihoods = sampler$setup$likelihood$density(samples)
    ml = logSumExp(likelihoods, mean = T)
    out <- list(marginalLikelihod=ml, method ="Prior")
  }
  
  warning("Note to the user: be aware that marginal likelihood calculations are notoriously prone to numerical stability issues. Especially in high-dimensional parameter spaces, there is no guarantee that the algorithms implemented in this function converge in all cases. Proceed at your own risk!")
  
  return(out)
}  







