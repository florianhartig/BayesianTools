#' Simulation-based calibration tests
#' 
#' This function data averaged posterior
#' 
#' @param posteriorList a list with posterior samples. List items must be of a class that is supported by \code{\link{getSample}}
#' @param priorDraws a matrix with parameter values, drawn from the prior, that were used to simulate the data underlying the posteriorList. If colnames are provided, these will be used in the plots
#' @param ... arguments to be passed to \code{\link{getSample}}. Consider in particular the thinning option. 
#'
#' @details The purpose of this function is to evaluate the results of a simulation-based calibration of an MCMC analysis. 
#' 
#' Briefly, the idea is to repeatedly
#' 
#' 1. sample parameters from the prior, 
#' 2. simulate new data based on these parameters, 
#' 3. calculate the posterior for these data
#' 
#' If the sampler and the likelihood are implemented correctly, the average of over all the posterior distribution should then again yield the prior (e.g. Cook et al., 2006). 
#' 
#' To test if this is the case, we implement the methods suggested by Talts et al., which is to calculate the rank statistics between the parameter draws and the posterior draws, which we then formally evaluate with a qq unif plot, and a ks.test
#' 
#' We speculate that a ks.test between the two distribution would likely give an identical result, but this is not noted in Talts et al.
#' 
#' 
#' Cook, S. R., Gelman, A. and Rubin, D. B. (2006). Validation of Software for Bayesian Models Using Posterior Quantiles. J. Comput. Graph. Stat. 15 675-692.
#' 
#' Talts, Sean, Michael Betancourt, Daniel Simpson, Aki Vehtari, and Andrew Gelman. "Validating Bayesian Inference Algorithms with Simulation-Based Calibration." arXiv preprint arXiv:1804.06788 (2018).
#' 
#' @export
#' 
calibrationTest <- function(posteriorList, priorDraws, ...){
  
  x = mergeChains(posteriorList, ...)
  
  nPar <- ncol(x)
  
  oldPar <- par(mfrow = getPanels(nPar*3))
  
  res = numeric(nPar)
  names(res) = colnames(priorDraws)

  for(i in 1:nPar){
    
    lim = range(x[,i], priorDraws[,i])
    
    hist(x[,i], breaks = 50, freq = F, col = "#99000020", main = colnames(priorDraws)[i])
    hist(priorDraws[,i], breaks = 50, freq = F, col = "#00990020", add = T)
    
    cDens = ecdf(x[,i])
    rankDist <- cDens(priorDraws[,i]) 
    hist(rankDist, breaks = 50, freq = F)
    abline(h = 1, col = "red")
    
    gap::qqunif(rankDist,pch=2,bty="n", logscale = F, col = "black", cex = 0.6, main = colnames(priorDraws)[i], cex.main = 1)
    
    res[i] = ks.test(x[,i], priorDraws[,i])$p.value
    legend("topleft", c(paste("KS test: p=", round(res[i], digits = 5)), paste("Deviation ", ifelse(res[i] < 0.05, "significant", "n.s."))), text.col = ifelse(res[i] < 0.05, "red", "black" ), bty="n")     
  }

  par(oldPar)  
  
  out = list()
  out$statistic = NULL
  out$method = "ks.test on rank statistics posterior / parameters"
  out$alternative = "two.sided"
  out$p.value = res
  class(out) = "htest"
  

  
  return(out)
  
}


