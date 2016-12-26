#' Performs a one-factor-at-a-time sensitivity analysis for the posterior of a given bayesianSetup within the prior range.
#' @param bayesianSetup An object of class BayesianSetup
#' @param selection indices of selected parameters
#' @note This function can also be used for sensitivity analysis of an arbitrary output - just create a BayesianSetup with this output. 
#' @export
plotSensitivity <- function(bayesianSetup, selection = NULL){
  
  if (is.null(selection)) selection = 1:bayesianSetup$numPars
  n = length(selection)
  post = list()
  for (j in 1:n){
    lowS = bayesianSetup$prior$lower[selection]
    upS = bayesianSetup$prior$upper[selection]
    parSen <- seq(lowS[n], upS[n], len = 20)
    refPar <- bayesianSetup$prior$best[selection]
    post[[j]] <- rep(NA, 20)
    for (i in 1:20){
      parS <- refPar
      parS[n] = parSen[i]
      post[[j]][i] = bayesianSetup$posterior$density(parS)
    }
  }
}



# plot(x = parSen, y = post, main = rownames(refPars)[parSel][n], type = "l", col = "red")
#   abline(v = refPar[n])
# } 
# par(mfrow = c(2,3))
# for (i in 1:6) plotSensitivity(i)
