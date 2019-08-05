#' Performs a one-factor-at-a-time sensitivity analysis for the posterior of a given bayesianSetup within the prior range.
#' @author Florian Hartig
#' @param bayesianSetup An object of class BayesianSetup
#' @param selection indices of selected parameters
#' @param equalScale if T, y axis of all plots will have the same scale
#' @note This function can also be used for sensitivity analysis of an arbitrary output - just create a BayesianSetup with this output. 
#' @example /inst/examples/plotSensitivityHelp.R
#' @export
plotSensitivity <- function(bayesianSetup, selection = NULL, equalScale = T){
  
  if (is.null(selection)) selection = 1:bayesianSetup$numPars
  n = length(selection)
  
  post = list()
  lowS = bayesianSetup$prior$lower[selection]
  upS = bayesianSetup$prior$upper[selection]
  refPar = bayesianSetup$prior$best[selection]
  names = bayesianSetup$names[selection]
  fullRefPar <- bayesianSetup$prior$best
  
  minR = Inf
  maxR = -Inf
  
  for (j in 1:n){
    post[[j]] <- data.frame(par = seq(lowS[j], upS[j], len = 20), resp = rep(NA, 20))
    
    for (i in 1:20){
      parS <- refPar
      parS[j] = post[[j]]$par[i]
      parS2 = fullRefPar
      parS2[selection] = parS
      post[[j]]$resp[i] = bayesianSetup$posterior$density(parS2)
    }
    minR = min(minR, post[[j]]$resp)
    maxR = max(maxR, post[[j]]$resp)
  }
  
  oldPar = par(mfrow = BayesianTools:::getPanels(n))
  
  
  for (i in 1:n){
    if(equalScale == T) plot(resp~par, xlab = names[i], type = "l", col = "red", data = post[[i]], ylim = c(minR, maxR), ylab = "Response")
    else plot(resp~par, xlab = names[i], type = "l", col = "red", data = post[[i]], ylab = "Response")
    
    abline(v = refPar[i])
  }
  
  names(post) = names
  post$reference = refPar
  
  par(oldPar)
  return(post)
}




