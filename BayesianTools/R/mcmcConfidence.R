#' Calculate confidence region from an MCMC or similar sample
#' @author Florian Hartig
#' @param sampleMatrix matrix of outcomes. Could be parameters or predictions
#' @param quantiles quantiles to be calculated
#' @export
#' @seealso \code{\link{getPredictiveDistribution}} \cr
#'          \code{\link{getPredictiveIntervals}} \cr
#' 
getCredibleIntervals <- function(sampleMatrix, quantiles = c(0.025, 0.975)){

  x = matrix (ncol = ncol(sampleMatrix), nrow = length(quantiles))
  rownames(x) = quantiles
    
  for (i in 1:length(quantiles)){
    x[i,] = apply(sampleMatrix,2,function(x)quantile(x,probs=quantiles[i]))  
  } 
  return(x)
}


#' Calculates predictive distribution based on the parameters
#' @author Florian Hartig
#' @param parMatrix matrix of parameter values
#' @param model model / function to calculate predictions. Outcome should be a vector
#' @param numSamples number of samples to be drawn 
#' @details If numSamples is greater than the number of rows in parMatrix, or NULL, or FALSE, or less than 1 all samples in parMatrix will be used.
#' @export
#' @seealso \code{\link{getPredictiveIntervals}} \cr
#'          \code{\link{getCredibleIntervals}} \cr
getPredictiveDistribution<-function(parMatrix, model, numSamples = 1000){
  
  # Do thinning if wanted and neccessary
  if (numSamples != F && nrow(parMatrix) > 2*numSamples && !is.null(numSamples) && numSamples > 0){
    sel = round(seq(1,nrow(parMatrix), len = numSamples ))
    parMatrixSel = parMatrix[sel,]
  }else{
    parMatrixSel = parMatrix
  }
  
  # calculate predictions
  
  run1 = model(parMatrixSel[1,])
  
  out = matrix(NA, ncol = length(run1), 
               nrow = nrow(parMatrixSel))
  
  out[1,] = run1
  
  for (i in 2:nrow(parMatrixSel)){
    out[i,] = model(parMatrixSel[i,])
  }
  return(out)
}


#' Calculates Bayesian credible (confidence) and predictive intervals based on parameter sample
#' @author Florian Hartig
#' @param parMatrix matrix of parameter values
#' @param model model / function to calculate predictions. Outcome should be a vector
#' @param numSamples number of samples to be drawn
#' @param quantiles quantiles to calculate
#' @param error function with signature f(mean, par) that generates error expectations from mean model predictions. Par is a vector from the matrix with the parameter samples (full length). f needs to know which of these parameters are parameters of the error function. If supplied, will calculate also predictive intervals additional to credible intervals
#' @details If numSamples is greater than the number of rows in parMatrix, or NULL, or FALSE, or less than 1 all samples in parMatrix will be used.
#' @export
#' @seealso \code{\link{getPredictiveDistribution}} \cr
#'          \code{\link{getCredibleIntervals}} \cr
getPredictiveIntervals<-function(parMatrix, model, numSamples = 1000, quantiles = c(0.025, 0.975), error = NULL){
  out = list()
    
  # Posterior predictive credible interval
  pred = getPredictiveDistribution(parMatrix, model = model, numSamples = numSamples)
  out$posteriorPredictiveCredibleInterval = getCredibleIntervals(sampleMatrix = pred, quantiles = quantiles)
  
  # Posterior predictive prediction interval
  # Posterior predictive simulations
  if(!is.null(error)){
    
    predDistr = pred
    for (i in 1:nrow(predDistr)){
      predDistr[i,] = error(mean = pred[i,], par = parMatrix[i,]) 
    }
    predInt = getCredibleIntervals(sampleMatrix = predDistr, quantiles = quantiles)   
    out$posteriorPredictivePredictionInterval = predInt
    out$posteriorPredictiveSimulations = predDistr
  }
  return(out)
}




