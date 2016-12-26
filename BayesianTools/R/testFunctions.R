#' Multivariate normal likelihood 
#' @author Florian Hartig
#' @description Generates a 3 dimensional multivariate normal likelihood function.
#' @param mean vector with the three mean values of the distribution
#' @param sigma either a correlation matrix, or "strongcorrelation", or "no correlation"
#' @param sample should the function create samples
#' @param n number of samples to create
#' @param throwErrors parameter for test purpose
#' @import mvtnorm
#' @import emulator
#' @details 3-d multivariate normal density function with mean 2,4,0 and either strong correlation (default), or no correlation.  
#' @export
#' @seealso \code{\link{testDensityBanana}} \cr
#'          \code{\link{testLinearModel}}
#' @example inst/examples/generateTestDensityMultiNormalHelp.R
#' 
# @param x TODO not fully implemented yet! either the parameter vector if the function is used for density, or the number of replicates when sampling

generateTestDensityMultiNormal <- function(mean = c(0,0,0), sigma = "strongcorrelation", sample = F, n = 1, throwErrors = -1){
  if (sigma == "strongcorrelation"){
    m <- c(0.2, 0.3, 0.3503)
    sigma = emulator::corr.matrix(cbind(m),scales=1)    
  }else if (sigma == "no correlation"){
    sigma = diag(rep(1,3))
  }
  if (sample == F) out <- function(x) mvtnorm::dmvnorm(x, mean = mean, sigma = sigma, log=T)
  else out <- function(n) mvtnorm::rmvnorm(n=n, mean = mean, sigma = sigma)
  
  # for test purposes
  if(runif(1) < throwErrors ){
    stop("TestError")
  }
  
  return(out)
}


#' Banana-shaped density function 
#' @author Florian Hartig
#' @param p 2-dim parameter vector
#' @import mvtnorm
#' @note inspired from package FMEmcmc, sees to go back to Laine M (2008). Adaptive MCMC Methods with Applications in Environmental and Models. Finnish Meteorological Institute Contributions 69. ISBN 978-951-697-662-7.
#' examples()
#' @export
#' @seealso \code{\link{generateTestDensityMultiNormal}} \cr
#'          \code{\link{testLinearModel}}
testDensityBanana <- function (p){
  P <- c(p[1], p[2] - (p[1]^2+1))
  Cov <- matrix(nrow = 2, data = c(1, 0.9, 0.9, 1))
  return(mvtnorm::dmvnorm(P, mean = rep(0, length(P)), sigma = Cov, log = T))
}


#' Normal likelihood
#' @author Florian Hartig
#' @param x a parameter vector of arbitrary length
#' @param sum if likelihood should be summed or not
#' @export
testDensityNormal <- function(x, sum = T){
  if(sum == T) return(sum(dnorm(x, log = T)))
  else return(dnorm(x, log = T))
}


#' Fake model, returns a ax + b linear response to 2-param vector
#' @author Florian Hartig
#' @param x 2-dim parameter vector
#' @param env optional, environmental covariate
#' @example /inst/examples/testLinearModel.R
#' @export
#' @seealso \code{\link{generateTestDensityMultiNormal}} \cr
#'          \code{\link{testDensityBanana}}
testLinearModel <- function(x, env = NULL){
  if (is.null(env)) env = seq(-3,3,len = 20)
  x[1] * env + x[2]
}


#' Test function infinity ragged
#' @author Florian Hartig
#' @param x 2-dim parameter vector
#' @param error should error or infinity be returned
#' @export
#' @seealso \code{\link{generateTestDensityMultiNormal}} \cr
#'          \code{\link{testDensityBanana}}
testDensityInfinity <- function(x, error = F){
  if( error == F) return(ifelse ((3*sum(x^2) %% 3) < 0.3, -Inf, 1))
  else return(ifelse ((3*sum(x^2) %% 3) < 0.3, stop("testerror"), 1))
}
