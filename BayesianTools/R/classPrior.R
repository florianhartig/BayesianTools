#' Creates a standardized prior class
#' @param density Prior density
#' @param sampler Sampling function for density (optional)
#' @param lower vector with lower bounds of parameters
#' @param upper vector with upper bounds of parameter
#' @param best vector with "best" parameter values
#' @note min and max don't update the absolute values of the prior density. Hence, if absolute values are a concern, trucated distributions should be used for the density
#' @export
#' @seealso \code{\link{createPriorDensity}} \cr
#'          \code{\link{createBetaPrior}} \cr
#'          \code{\link{createUniformPrior}} \cr
#'          \code{\link{createTruncatedNormalPrior}}\cr
#'          \code{\link{createBayesianSetup}}\cr
#' @example /inst/examples/createPrior.R
createPrior <- function(density = NULL, sampler = NULL, lower = NULL, upper = NULL, best = NULL){
  
  # case density is a Bayesian Posterior
  if(inherits(density,"bayesianOutput")) return(createPriorDensity(density, lower = lower, upper = upper, best = best))
  
  if(! is.null(lower) & ! is.null(upper)) if(any(lower > upper)) stop("prior with lower values > upper")
  
  if(is.null(best) & ! is.null(lower) & ! is.null(upper)) best = (upper + lower) / 2
  
  # if no density is provided 
  if (is.null(density)){
    density <- function(x){
      return(0) 
    }
  }

  catchingPrior <- function(x){
    
    # check if bounds are respected
    if(!is.null(lower)){
      if (any(x < lower)) return(-Inf)
    } 
    if(!is.null(upper)){
      if (any(x > upper)) return(-Inf)
    }
    
    # calculate prior density within try-catch statement
    out <- tryCatch(
    {
      density(x)
    },
    error=function(cond) {
      warning("Problem in the prior", cond)
      return(-Inf)
    }
    )
    # extra check
    if (out == Inf) stop("Inf encountered in prior")

    return(out)
  }
  
  parallelDensity<- function(x){
    if (is.vector(x)) return(catchingPrior(x))
    else if(is.matrix(x)) return(apply(x, 1, catchingPrior))
    else stop("parameter must be vector or matrix")
  }
  
  # Check and parallelize the sampler
  if(!is.null(sampler)){
    npar <- length(sampler())
    parallelSampler <- function(n=NULL){
      if(is.null(n)) out = sampler()
      else{
        if (npar == 1) out = matrix(replicate(n, sampler()))
        else if (npar >1) out = t(replicate(n, sampler(), simplify = T))
        else stop("sampler provided doesn't work")
      } 
      return(out)
    }
  }
  else parallelSampler = function(n = NULL){
   stop("Attept to call the sampling function of the prior, although this function has not been provided in the Bayesian setup. A likely cause of this error is that you use a function or sampling algorithm that tries to sample from the prior. Either change the settings of your function, or provide a sampling function in your BayesianSetup (see ?createBayesianSetup, and ?createPrior)") 
  }
  

  
  out<- list(density = parallelDensity, sampler = parallelSampler, lower = lower, upper = upper, best = best, originalDensity = density)
  class(out) <- "prior"
  return(out)
}


#' Convenience function to create a simple uniform prior distribution
#' @param lower vector of lower prior range for all parameters
#' @param upper vector of upper prior range for all parameters
#' @param best vector with "best" values for all parameters
#' @note for details see createPrior
#' @seealso \code{\link{createPriorDensity}}, \code{\link{createPrior}}, \code{\link{createBetaPrior}}, \code{\link{createTruncatedNormalPrior}}, \code{\link{createBayesianSetup}} 
#' @example /inst/examples/createUniformPrior.R
#' @export
createUniformPrior<- function(lower, upper, best = NULL){
  len = length(lower)
  density <- function(x){
    if (length(x) != len) stop("parameter vector does not match prior")
    else return(sum(dunif(x, min = lower, max = upper, log = T)))
  }
  sampler <- function() runif(length(lower), lower, upper)
  
  out <- createPrior(density = density, sampler = sampler, lower = lower, upper = upper, best = best)
  return(out)
}


#' Convenience function to create a truncated normal prior
#' 
#' @param mean best estimate for each parameter
#' @param sd sdandard deviation
#' @param lower vector of lower prior range for all parameters
#' @param upper vector of upper prior range for all parameters
#' @note for details see createPrior
#' @seealso \code{\link{createPriorDensity}} \cr
#'          \code{\link{createPrior}} \cr
#'          \code{\link{createBetaPrior}} \cr
#'          \code{\link{createUniformPrior}} \cr
#'          \code{\link{createBayesianSetup}} \cr
#' @export
#' @example /inst/examples/createTruncatedNormalPrior.R

createTruncatedNormalPrior<- function(mean, sd, lower, upper){
  len = length(mean)
  density <- function(x){
    if (length(x) != len) stop("parameter vector does not match prior")
    else return(sum(msm::dtnorm(x, mean = mean, sd = sd, lower = lower, upper = upper, log = T)))
  }
  sampler <- function(){
    msm::rtnorm(n = length(mean), mean = mean, sd = sd, lower = lower, upper = upper)
  }
  out <- createPrior(density = density, sampler = sampler, lower = lower, upper = upper)
  return(out)
}


#' Convenience function to create a beta prior
#' 
#' @param a shape1 of the beta
#' @param b shape2 of the beta
#' @param upper upper values for the parameters
#' @param lower lower values for the parameters
#' @note for details see createPrior
#' @details This creates a beta prior, assuming that lower / upper values for parameters are are fixed. The beta is the calculated relative to this lower / upper space. 
#' @seealso \code{\link{createPriorDensity}} \cr
#'          \code{\link{createPrior}} \cr
#'          \code{\link{createTruncatedNormalPrior}} \cr
#'          \code{\link{createUniformPrior}} \cr
#'          \code{\link{createBayesianSetup}} \cr
#' @export
createBetaPrior<- function(a, b, lower=0, upper=1){
  len = length(a)
  if (! any(upper > lower)) stop("wrong values in beta prior")
  range = upper - lower
  density <- function(x){
    x = (x - lower) / range
    if (length(x) != len) stop("parameter vector does not match prior")
    else return(sum( dbeta(x, shape1 = a, shape2 = b, log=T) ))
  }
  sampler <- function(){
    out = rbeta(n = len, shape1 = a, shape2 = b)
    out = (out * range) + lower
    return(out)
  }
  out <- createPrior(density = density, sampler = sampler, lower = NULL, upper = NULL)
  return(out)
}


#' Fits a density function to a multivariate sample
#' @export
#' @param x a matrix or BayesianOutput
#' @param method method to generate prior - default and currently only option is multivariate
#' @param eps numerical
#' @param lower vector with lower bounds of parameter
#' @param upper vector with upper bounds of parameter
#' @param best vector with "best" values of parameter
#' @seealso \code{\link{createPrior}} \cr
#'          \code{\link{createBetaPrior}} \cr
#'          \code{\link{createTruncatedNormalPrior}} \cr
#'          \code{\link{createUniformPrior}} \cr
#'          \code{\link{createBayesianSetup}} \cr
#' @example /inst/examples/createPriorDensity.R
createPriorDensity <- function(x, method = "multivariate", eps = 1e-10, lower = NULL, upper = NULL, best = NULL){
  
  if(inherits(x,"bayesianOutput")) x = getSample(x)
  
  
  if(method == "multivariate"){
    nPars = ncol(x)
    covariance = cov(x)
    covar = as.matrix(Matrix::nearPD(covariance + diag(eps, nPars))$mat)
    
    density = function(par){
      dens = mvtnorm::dmvnorm(x = par, sigma = covar, log = T)
      return(dens)
    }
    
    sampler = function(n=1){
      par = as.vector(mvtnorm::rmvnorm(n = n, sigma = covar))
      return(par)
    }
    
    out <- createPrior(density = density, sampler = sampler, lower = lower, upper = upper, best = best)
    return(out)
  }
}


