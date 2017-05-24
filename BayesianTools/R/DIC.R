
# Modified from https://gist.github.com/gaberoo/4619102

##############################################################################
# Estimate Deviance Information Criterion (DIC)
#
# References:
#   Bayesian Data Analysis.
#   Gelman, A., Carlin, J., Stern, H., and Rubin D.
#   Second Edition, 2003
#
#   Bayesian predictive information criterion for the evaluation of 
#     hierarchical Bayesian and empirical Bayes models. 
#   Ando, T.
#   Biometrika, 2007 
#
# Input:
#   x       : matrix of posterior samples
#   lik     : vector of the likelihood of the posterior samples
#   lik.fun : function that calculates the likelihood
#   ...     : other parameters that are passed to 'lik.fun'
#
# Output:
#   list()
#     DIC   : Deviance Information Criterion
#     IC    : Bayesian Predictive Information Criterion
#     pD    : Effective number of parameters (pD = Dbar - Dhat)
#     pV    : Effective number of parameters (pV = var(D)/2)
#     Dbar  : Expected value of the deviance over the posterior
#     Dhat  : Deviance at the mean posterior estimate
##############################################################################


#' Deviance information criterion
#' @author Florian Hartig
#' @param sampler An object of class bayesianOutput (mcmcSampler, smcSampler, or mcmcList) 
#' @param ... further arguments passed to \code{\link{getSample}}
#' @references Spiegelhalter, D. J.; Best, N. G.; Carlin, B. P. & van der Linde, A. (2002) Bayesian measures of model complexity and fit. J. Roy. Stat. Soc. B, 64, 583-639.\cr\cr
#' Gelman, A.; Hwang, J. & Vehtari, A. (2014) Understanding predictive information criteria for Bayesian models. Statistics and Computing, Springer US, 24, 997-1016-.
#' @details Output:
#'   list with the following elements: \cr
#'     DIC   : Deviance Information Criterion \cr
#'     IC    : Bayesian Predictive Information Criterion \cr
#'     pD    : Effective number of parameters (pD = Dbar - Dhat) \cr
#'     pV    : Effective number of parameters (pV = var(D)/2) \cr
#'     Dbar  : Expected value of the deviance over the posterior \cr
#'     Dhat  : Deviance at the mean posterior estimate \cr
#' @seealso \code{\link{WAIC}}, \code{\link{MAP}}, \code{\link{marginalLikelihood}}
#' @export
DIC <- function(sampler, ...){
  
  draw = getSample(sampler, parametersOnly = F, ...)
  
  if(class(sampler)[1] %in% c("mcmcSamplerList", "smcSamplerList")) sampler = sampler[[1]]
  
  x = draw[,1:sampler$setup$numPars]
  lik = draw[,sampler$setup$numPars+2]
  lik.fun = sampler$setup$likelihood$density
  
  D.bar <- -2*mean(lik)
  if(is.vector(x)) theta.bar = mean(x) else theta.bar <- apply(x,2,mean)
  D.hat <- -2*lik.fun(theta.bar)
  pD <- D.bar - D.hat
  pV <- var(-2*lik)/2
  list(DIC=pD+D.bar,IC=2*pD+D.bar,pD=pD,pV=pV,Dbar=D.bar,Dhat=D.hat)
}
