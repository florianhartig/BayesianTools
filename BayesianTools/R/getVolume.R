#' Calculate posterior volume
#' 
#' @param sampler an object of superclass bayesianOutput or any other class that has the getSample function implemented (e.g. Matrix)
#' @param prior schould also prior volume be calculated
#' @param method method for volume estimation. Currently, the only option is "MVN"
#' @param ... additional parameters to pass on to the \code{\link{getSample}}
#' @details The idea of this function is to provide an estimate of the "posterior volume", i.e. how "broad" the posterior is. One potential application is to the overall reduction of parametric uncertainty between different data types, or between prior and posterior.  
#' 
#' Implemented methods for volume estimation:
#' 
#' Option "MVN" - in this option, the volume is calculated as the determinant of the covariance matrix of the prior / posterior sample. 
#' 
#' @example /inst/examples/getVolume.R
#' @export
getVolume <- function(sampler, prior = F, method = "MVN", ...){
    
    x = getSample(sampler, ...)
    
    if(method == "MVN"){
      nPars = ncol(x)
      postVol = det(cov(x))
    }else stop("BayesianTools: unknown method argument in getVolume")
      
    if(prior == T){
      x = out$setup$prior$sampler(5000)
      
      if(method == "MVN"){
        nPars = ncol(x)
        priorVol = det(cov(x))
      }else stop("BayesianTools: unknown method argument in getVolume")
      return(list(priorVol = priorVol, postVol = postVol))
    }else return(postVol)
}
  

