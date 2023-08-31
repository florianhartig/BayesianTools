#' Creates a standardized likelihood class#'
#' @author Florian Hartig
#' @param likelihood log likelihood density
#' @param names parameter names (optional)
#' @param parallel parallelization , either i) no parallelization --> F, ii) native R parallelization --> T / "auto" will select n-1 of your available cores, or provide a number for how many cores to use, or iii) external parallelization --> "external". External means that the likelihood is already able to execute parallel runs in the form of a matrix. 
#' @param catchDuplicates logical, determines whether unique parameter combinations should only be evaluated once. This is only applicable when the likelihood accepts a matrix with parameters as columns. 
#' @param parallelOptions a list containing two lists. First, "packages" specifies the R packages necessary to run the likelihood function. Second, "objects" contains the objects in the global environment needed to run the likelihood function (for details see \code{\link{createBayesianSetup}}).
#' @param sampler sampler
#' @seealso \code{\link{likelihoodIidNormal}} \cr
#'          \code{\link{likelihoodAR1}} \cr
#' @export
createLikelihood <- function(likelihood, names = NULL, parallel = F, catchDuplicates=T, 
                             sampler = NULL, parallelOptions = NULL){
  
  # check if point-wise likelihood available
  pwLikelihood = if ("sum" %in% names(as.list(args(likelihood)))) TRUE else FALSE
  
  catchingLikelihood <- function(x, ...){
    out <- tryCatch(
    {
      y = likelihood(x, ...)
      if (any(y == Inf | is.nan(y) | is.na(y) | !is.numeric(y))){
        message(paste("BayesianTools warning: positive Inf or NA / nan values, or non-numeric values occured in the likelihood. Setting likelihood to -Inf.\n Original value was", y, "for parameters", x, "\n\n "))
        y[is.infinite(y) | is.nan(y) | is.na(y) | !is.numeric(y)] = -Inf
      }
      y 
    },
    error=function(cond){
      cat(c("Parameter values ", x, "\n"))
      message("Problem encountered in the calculation of the likelihood with parameter ", x, "\n Error message was", cond, "\n set result of the parameter evaluation to -Inf ", "ParameterValues ")
      return(-Inf)
    }
        )
    return(out)
  }

  # initalize cl 
  cl <- NULL
  
  if (parallel == T | parallel == "auto" | is.numeric(parallel)) {
    tmp <- generateParallelExecuter(likelihood, parallel, parallelOptions) 
    parallelLikelihood <- tmp$parallelFun
    cl <- tmp$cl
    parallel = T
  }

  
  parallelDensity<- function(x, ...){
    if (is.vector(x)) return(catchingLikelihood(x, ...))
    else if(is.matrix(x)){
      if(catchDuplicates == TRUE){
        # Check for the rows that are not duplicated
        wn <- which(!duplicated(x))
        if(length(wn) <2) {
          return(parallelLikelihood(x, ...)) }
        else {
        # Define a output vector 
        out1 <- rep(0,length=nrow(x))
       
        # Run the likelihood function for unique values
        if (parallel == "external"){ 
          out1[wn]<-likelihood(x[wn,], ...)
          }
        else{
          if (parallel == T){ 
            out1[wn]<-parallelLikelihood(x[wn,], ...)
          }
            else{
              out1[wn]<-apply(x[wn,], 1, likelihood, ...)   
            }
        }
        # Copy the values for the duplicates
        for(i in 1:length(out1)){
          if(out1[i] != 0) next
          else{
            same <- numeric()
            for(k in 1:length(out1)){
              if(all(x[k,]== x[i,])){
                same <- c(same,k)
              }
            }
            out1[same[-1]] <- out1[same[1]]
          }
        }
      
        return(out1)
        }}
      
      else{
      if (parallel == "external") return(likelihood(x, ...))
      else if (parallel == T){
      return(parallelLikelihood(x, ...))}
      else return(apply(x, 1, likelihood, ...))   
 
      }
    }
    else stop("parameter must be vector or matrix")
  }
  out<- list(density = parallelDensity, sampler = sampler, cl = cl, pwLikelihood = pwLikelihood, parNames = names)
  class(out) <- "likelihood"
  return(out)
}



#library(mvtnorm)
#library(sparseMVN)

#' Normal / Gaussian Likelihood function
#' @author Florian Hartig
#' @param predicted vector of predicted values
#' @param observed vector of observed values
#' @param sd standard deviation of the i.i.d. normal likelihood
#' @export
likelihoodIidNormal <- function(predicted, observed, sd){
  notNAvalues = !is.na(observed)
  if (sd <= 0) return (-Inf)
  else return(sum(dnorm(predicted[notNAvalues], mean = observed[notNAvalues], sd = sd, log = T)))
}

# TODO - gibbs sample out the error terms 

#' AR1 type likelihood function
#' @author Florian Hartig
#' @param predicted vector of predicted values
#' @param observed vector of observed values
#' @param sd standard deviation of the iid normal likelihood
#' @param a temporal correlation in the AR1 model
#' @note The AR1 model considers the process: \cr y(t) = a y(t-1) + E \cr e = i.i.d. N(0,sd) \cr |a| < 1 \cr At the moment, no NAs are allowed in the time series.
#' @export
likelihoodAR1 <- function(predicted, observed, sd, a){
  if (any(is.na(observed))) stop("AR1 likelihood cannot work with NAs included, split up the likelihood")
  if (sd <= 0) return (-Inf)
  if (abs(a) >= 1) return (-Inf)
  
  n = length(observed)
  
  res = predicted - observed
  
  # this calculates the unconditiona LL for this data, see e.g. http://stat.unicas.it/downloadStatUnicas/seminari/2008/Julliard0708_1.pdf
  
  ll =  0.5 * (  - n * log(2*pi)
                 - n * log(sd^2) 
                 + log( 1- a^2 )
                 - (1- a^2) / sd^2 * res[1]^2
                 - 1 / sd^2 * sum( (res[2:n] - a * res[1:(n-1)])^2)
                )
  return(ll)
}
# Tests
# library(stats)
# data<-arima.sim(n=1000,model = list(ar=0.9))
# x <- ar(data, aic = F, order.max = 1)
# opt <- function(par){
#   -likelihoodAR1(data, rep(0,1000), sd = par[1], a = par[2] )
# }
# optim(c(1.1,0.7), opt  )









