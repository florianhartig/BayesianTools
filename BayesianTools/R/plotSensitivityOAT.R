#' Performs a one-factor-at-a-time sensitivity analysis for the posterior of a given bayesianSetup within the prior range.
#' @author Florian Hartig
#' @param bayesianSetup An object of class BayesianSetup
#' @param selection indices of selected parameters
#' @param equalScale if T, y axis of all plots will have the same scale
#' @note This function can also be used for sensitivity analysis of an arbitrary output - just create a BayesianSetup with this output. 
#' @example #' @example /inst/examples/plotSensitivityHelp.R
#' @export
plotSensitivity <- function(bayesianSetup, selection = NULL, equalScale = T){
  
  if (is.null(selection)) selection = 1:bayesianSetup$numPars
  n = length(selection)
  
  if(type == "OAT"){
    post = list()
    lowS = bayesianSetup$prior$lower[selection]
    upS = bayesianSetup$prior$upper[selection]
    refPar = bayesianSetup$prior$best
    
    minR = Inf
    maxR = -Inf
    
    for (j in 1:n){
      post[[j]] <- data.frame(par = seq(lowS[j], upS[j], len = 20), resp = rep(NA, 20))
      for (i in 1:20){
        parS <- refPar
        parS[n] = post[[j]]$par[i]
        post[[j]]$resp[i] = bayesianSetup$posterior$density(parS)
      }
      minR = min(minR, post[[j]]$resp)
      maxR = max(maxR, post[[j]]$resp)
    }
    
    oldPar = par(mfrow = BayesianTools:::getPanels(n))
    
    
    for (i in 1:n){
      plot(resp~par, xlab = bayesianSetup$names[n], type = "l", col = "red", data = post[[j]])
      abline(v = refPar[n])
    }

    par(oldPar)
    return(post)
  }

}




