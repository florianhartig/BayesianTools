
# Motivation for this functions from 
# https://radfordneal.wordpress.com/2008/08/17/the-harmonic-mean-of-the-likelihood-worst-monte-carlo-method-ever/
# https://gist.github.com/gaberoo/4619102


#  ' @export
#marginalLikelihood <- function(x,lik,V,sampler$setup$likelihood$density,sampler$setup$prior$density,..., num.samples=1000,log=TRUE) UseMethod("marginalLikelihood")

#' Calcluated the marginal likelihood from a set of MCMC samples
#' @export
#' @author Florian Hartig
#' @param sampler an MCMC or SMC sampler or list, or for method "Prior" also a BayesianSetup
#' @param numSamples number of samples to use. How this works, and if it requires recalculating the likelihood, depends on the method
#' @param method method to choose. Currently available are "Chib" (default), the harmonic mean "HM", sampling from the prior "Prior", and bridge sampling "Bridge". See details
#' @param ... further arguments passed to \code{\link{getSample}}
#' @details The marginal likelihood is the average likelihood across the prior space. It is used, for example, for Bayesian model selection and model averaging. 
#' 
#' It is defined as \deqn{ML = \int L(\Theta) p(\Theta) d\Theta}
#' 
#' Given that MLs are calculated for each model, you can get posterior weights (for model selection and/or model averaging) on the model by 
#' 
#' \deqn{P(M_i|D) = ML_i * p(M_i) / (\sum_i ML_i * p(M_i) )}
#' 
#' In BT, we return the log ML, so you will have to exp all values for this formula. 
#' 
#' It is well-known that the ML is VERY dependent on the prior, and in particular the choice of the width of uninformative priors may have major impacts on the relative weights of the models. It has therefore been suggested to not use the ML for model averaging / selection on uninformative priors. If you have no informative priors, and option is to split the data into two parts, use one part to generate informative priors for the model, and the second part for the model selection. See help for an example. 
#' 
#' The marginalLikelihood function currently implements four ways to calculate the marginal likelihood. Be aware that marginal likelihood calculations are notoriously prone to numerical stability issues. Especially in high-dimensional parameter spaces, there is no guarantee that any of the implemented algorithms will converge reasonably fast. The recommended (and default) method is the method "Chib" (Chib and Jeliazkov, 2001), which is based on MCMC samples, with a limited number of additional calculations. Despite being the current recommendation, note there are some numeric issues with this algorithm that may limit reliability for larger dimensions.
#'   
#'  The harmonic mean approximation, is implemented only for comparison. Note that the method is numerically unrealiable and usually should not be used. 
#' 
#' The third method is simply sampling from the prior. While in principle unbiased, it will only converge for a large number of samples, and is therefore numerically inefficient. 
#' 
#' The Bridge method uses bridge sampling as implemented in the R package "bridgesampling". It is potentially more exact than the Chib method, but might require more computation time. However, this may be very dependent on the sampler.
#' 
#' @return A list with log of the marginal likelihood, as well as other diagnostics depending on the chose method
#'    
#' @example /inst/examples/marginalLikelihoodHelp.R
#' @references 
#' 
#' Chib, Siddhartha, and Ivan Jeliazkov. "Marginal likelihood from the Metropolis-Hastings output." Journal of the American Statistical Association 96.453 (2001): 270-281.
#' 
#' Dormann et al. 2018. Model averaging in ecology: a review of Bayesian, information-theoretic, and tactical approaches for predictive inference. Ecological Monographs
#' 
#' @seealso \code{\link{WAIC}}, \code{\link{DIC}}, \code{\link{MAP}}
marginalLikelihood <- function(sampler, numSamples = 1000, method = "Chib", ...){
  

  if ((class(sampler)[1] %in% c("mcmcSamplerList", "smcSamplerList"))) {
    setup <- sampler[[1]]$setup
    posterior = sampler[[1]]$setup$posterior$density 
  } else if ((class(sampler)[1] %in% c("mcmcSampler", "smcSampler"))) {
    setup <- sampler$setup
    posterior = sampler$setup$posterior$density 
  } else if ((class(sampler)[1] %in% c("BayesianSetup"))) {
    setup <- sampler
    posterior = sampler$posterior$density 
  } else stop("sampler must be a sampler or a BayesianSetup")
  
  
  if (method == "Chib"){
    
    chain <- getSample(sampler = sampler, parametersOnly = F, ...)
    
    if(class(sampler)[1] %in% c("mcmcSamplerList", "smcSamplerList")) sampler <- sampler[[1]]
    
    x <- chain[,1:sampler$setup$numPars,drop=F] #  RB: drop=F
    
    lik <- chain[,sampler$setup$numPars + 2]
    MAPindex <- which.max(chain[,sampler$setup$numPars + 1])
    
    #propGen = createProposalGenerator(covariance = cov(x))
    
    V <- cov(x)
    
    # calculate reference parameter 
    
    theta.star <- x[MAPindex,,drop=F]
    lik.star <- lik[MAPindex]
    
    # get samples from posterior
    
    g <- sample.int(nrow(x), numSamples, replace=TRUE) # should replace really be true?
    q.g <- mvtnorm::dmvnorm(x[g,,drop=F], mean = theta.star, sigma = V, log = FALSE) # RB: drop=F
    lik.g <- lik[g]
    alpha.g <- sapply(lik.g, function(l) min(1, exp(lik.star - l))) # Metropolis Ratio
    
    #lik.g <- apply(theta.g,1,sampler$setup$likelihood$density,...)
    
    
    # get samples from proposal
    theta.j <- mvtnorm::rmvnorm(numSamples, mean = theta.star, sigma = V)
    lik.j <- apply(theta.j, 1, sampler$setup$likelihood$density)
    alpha.j <- sapply(lik.j, function(l) min(1, exp(l - lik.star)))  # Metropolis Ratio
    
    # Prior 
    pi.hat <- mean(alpha.g * q.g) / mean(alpha.j)
    pi.star <- 0
    
    if (!is.null(sampler$setup$prior$density)) pi.star <- sampler$setup$prior$density(theta.star)
    ln.m <- lik.star + pi.star - log(pi.hat)
    
    out <- list(ln.ML = ln.m, ln.lik.star = lik.star, ln.pi.star = pi.star, ln.pi.hat = log(pi.hat), method = "Chib")
    
  } else if (method == "HM"){
    
    warning("The Harmonic Mean estimator is notoriously unstable. It's only implemented for comparison. We strongly advice against using it for research!")
    
    chain <- getSample(sampler = sampler, parametersOnly = F, ...)
    lik <- chain[, setup$numPars + 2]
    ml <- log(1 / mean(1 / exp(lik)))
    # ml = 1 / logSumExp(-lik, mean = T) function needs to be adjusted
    out <- list(ln.ML=ml, method ="HM")
    
  } else if (method == "Prior"){
    
    samples <- setup$prior$sampler(numSamples)
    likelihoods <- setup$likelihood$density(samples)
    
    ml <- logSumExp(likelihoods, mean = T)
    out <- list(ln.ML=ml, method ="Prior")
    
  } else if (method == "Bridge") {
    
    chain <- getSample(sampler = sampler, parametersOnly = F, numSamples = numSamples, ...)
    
    nParams <- setup$numPars
    lower <- setup$prior$lower
    upper <- setup$prior$upper
    
    
    out <- list(ln.ML = bridgesample(chain ,nParams, lower, upper, posterior)$logml, method ="Bridge")
    
  } else if ("NN") {
    
    # TODO: implement nearest neighbour method:
    # https://arxiv.org/abs/1704.03472
    stop("Not yet implemented")
    
  } else {
    stop(paste(c("\"", method, "\" is not a valid method parameter!"), sep = " ", collapse = ""))
  }
  
  return(out)
}  


#' Calculates the marginal likelihood of a chain via bridge sampling
#' @export
#' @author Tankred Ott
#' @param chain a single mcmc chain with samples as rows and parameters and posterior density as columns.
#' @param nParams number of parameters
#' @param lower optional - lower bounds of the prior
#' @param upper optional - upper bounds of the prior
#' @param posterior posterior density function
#' @param ... arguments passed to bridge_sampler
#' @details This function uses "bridge_sampler" from the package "bridgesampling".
#' @example /inst/examples/bridgesampleHelp.R
#' @seealso \code{\link{marginalLikelihood}}
#' @keywords internal
bridgesample <- function (chain, nParams, lower = NULL, upper = NULL, posterior, ...) {
  # TODO: implement this without bridgesampling package
  # https://github.com/quentingronau/bridgesampling
  if (is.null(lower)) lower <- rep(-Inf, nParams)
  if (is.null(upper)) upper <- rep(Inf, nParams)
  
  names(lower) <- names(upper) <- colnames(chain[, 1:nParams])
  
  log_posterior = function(x, data){
    return(posterior(x))
  }
  
  out <- bridgesampling::bridge_sampler(
    samples = chain[, 1:nParams],
    log_posterior = log_posterior,
    data = chain,
    lb = lower,
    ub = upper,
    ...
  )
  
  return(out)
}





