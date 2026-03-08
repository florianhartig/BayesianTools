
#' Multivariate normal likelihood 
#' @author Florian Hartig
#' @description Generates a 3 dimensional multivariate normal likelihood function.
#' @param mean vector with the three mean values of the distribution
#' @param sigma either a correlation matrix, or "strongcorrelation", or "no correlation"
#' @param sample logical, should the function create samples?
#' @param n number of samples to create
#' @param throwErrors parameter for test purpose. Between 0 and 1 for proportion of errors
#' @details 3-d multivariate normal density function with mean 2,4,0 and either strong correlation (default), or no correlation.  
#' @export
#' @seealso \code{\link{testDensityBanana}} \cr
#'          \code{\link{testLinearModel}}
#' @example inst/examples/generateTestDensityMultiNormalHelp.R
#' 
# @param x TODO not fully implemented yet! either the parameter vector if the function is used for density, or the number of replicates when sampling

generateTestDensityMultiNormal <- function(mean = c(0,0,0), sigma = "strongcorrelation", sample = F, n = 1, throwErrors = -1){
  # for test purposes
  if(runif(1) < throwErrors ){
    stop("TestError")
  }
  
  if (sigma == "strongcorrelation"){
    m <- c(0.2, 0.3, 0.3503)
    sigma = emulator::corr.matrix(cbind(m),scales=1)    
  }else if (sigma == "no correlation"){
    sigma = diag(rep(1,3))
  }
  if (sample == F){
    out <- function(x) mvtnorm::dmvnorm(x, mean = mean, sigma = sigma, log=T)
    return(out)
  }else{
     out_sample <- function(n) mvtnorm::rmvnorm(n=n, mean = mean, sigma = sigma)
     return(out_sample)
  }
}


#' Banana-shaped density function 
#' @author Florian Hartig
#' @param p 2-dim parameter vector
#' @note inspired from package FMEmcmc, seems to go back to Laine M (2008). Adaptive MCMC Methods with Applications in Environmental and Models. Finnish Meteorological Institute Contributions 69. ISBN 978-951-697-662-7.
#' @export
#' @seealso \code{\link{generateTestDensityMultiNormal}} \cr
#'          \code{\link{testLinearModel}}
testDensityBanana <- function (p){
  P <- c(p[1], p[2] - (p[1]^2+1))
  Cov <- matrix(nrow = 2, data = c(1, 0.9, 0.9, 1))
  return(mvtnorm::dmvnorm(P, mean = rep(0, length(P)), sigma = Cov, log = T))
}

#' GelmanMeng test function
#' 
#' @param x parameter vector
#' @param A function parameter
#' @param B function parameter
#' @param C1 function parameter
#' @param C2 function parameter
#' @param log log
#' 
#' A non-elliptical, bivariate density function proposed by Gelman and Meng (1991). 
testDensityGelmanMeng <- function(x, A = 1, B = 0, C1 = 3, C2 = 3, log = TRUE) {
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  r <- -0.5 * (A * x[,1]^2 * x[,2]^2 + x[,1]^2 + x[,2]^2
      - 2 * B * x[,1] * x[,2] - 2 * C1 * x[,1] - 2 * C2 * x[,2])
  if (!log) r <- exp(r)
  as.vector(r)
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


#' 3d Mutivariate Normal likelihood
#' @param x a parameter vector of arbitrary length
#' @param sigma either a correlation matrix, or "strongcorrelation", or "no correlation"
#' @export
testDensityMultiNormal <- function(x, sigma = "strongcorrelation"){
  if (sigma == "strongcorrelation"){
    m <- c(0.2, 0.3, 0.3503)
    sigma = emulator::corr.matrix(cbind(m),scales=1)    
  }else if (sigma == "no correlation"){
    sigma = diag(rep(1,3))
  }
  return(mvtnorm::dmvnorm(x, mean = c(0,0,0), sigma = sigma, log=T))

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
