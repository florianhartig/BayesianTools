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
#' @param thin thinning 
#' @export
#' @seealso \code{\link{getPredictiveIntervals}} \cr
#'          \code{\link{getCredibleIntervals}} \cr
getPredictiveDistribution<-function(parMatrix, model, thin = 1000){
  
  # Do thinning if wanted and neccessary
  if (thin != F & nrow(parMatrix) > 2*thin){
    sel = round(seq(1,nrow(parMatrix), len = thin ))
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
#' @param thin thinning 
#' @param quantiles quantiles to calculate
#' @param error function that accepts a vector as in observation. If supplied, will calculate also predictive intervals additional to credible intervals
#' @export
#' @seealso \code{\link{getPredictiveDistribution}} \cr
#'          \code{\link{getCredibleIntervals}} \cr
getPredictiveIntervals<-function(parMatrix, model, thin = 1000, quantiles = c(0.025, 0.975), error = NULL){
  pred = getPredictiveDistribution(parMatrix, model = model, thin = thin)
  out = getCredibleIntervals(sampleMatrix = pred, quantiles = quantiles)
  
  if(!is.null(error)){
    
    predDistr = pred
    for (i in 1:nrow(predDistr)){
      predDistr[i,] = error(mean = pred[i,], par = parMatrix[i,]) 
    }

    predInt = getCredibleIntervals(sampleMatrix = predDistr, quantiles = quantiles)   
    out = rbind(out, predInt)
  }
  
  return(out)
}




