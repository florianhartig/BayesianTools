#' Creates a standardized collection of prior, likelihood and posterior functions, including error checks etc.
#' @param likelihood log likelihood density function
#' @param prior log prior density
#' @param priorSampler prior sampling function
#' @param lower vector with lower prior limits
#' @param upper vector with upper prior limits
#' @param best vector with best values
#' @param names optional vector with parameter names
#' @param parallel parallelization option. Default is F. Other options are T, or "external". See details.
#' @param parallelOptions list containing three lists. First "packages" determines the R packages necessary to run the likelihood function. Second "variables" the objects in the global envirnment needed to run the likelihood function and third "dlls" the DLLs needed to run the likelihood function (see Details and Examples). 
#' @param catchDuplicates Logical, determines whether unique parameter combinations should only be evaluated once. Only used when the likelihood accepts a matrix with parameter as columns. 
#' @details For parallelization, option T means that an automatic parallelization via R is attempted, or "external", in which case it is assumed that the likelihood is already parallelized. In this case it needs to accept a matrix with parameters as columns.
#' Further you can specify the packages, objects and DLLs that are exported to the cluster. 
#' By default a copy of your workspace is exported. However, depending on your workspace this can be very inefficient.
#'
#' 
#' @export
#' @seealso \code{\link{checkBayesianSetup}} \cr
#'          \code{\link{createLikelihood}} \cr
#'          \code{\link{createPrior}} \cr
#' @example /inst/examples/classBayesianSetup.R
#' 
#' 
#@param model TODO
createBayesianSetup <- function(likelihood, 
                                prior = NULL, 
                                priorSampler = NULL, 
                                parallel = FALSE,
                                lower= NULL, 
                                upper = NULL, 
                                best = NULL, 
                                names = NULL, 
                                parallelOptions = NULL, 
                                catchDuplicates = FALSE){
  
  # TODO implement parameter "model" (function that makes predictions from the model)
  model <- NULL
  
  
   if(is.null(upper) && is.null(lower) && is.null(prior)) stop("Either boundaries or prior density and prior sampler must be provided.")
  
  
  if(is.null(parallelOptions)) parallelOptions <- list(variables = "all", packages = "all", dlls = "all")
  
  # PRIOR CHECKS
  if(inherits(prior,"bayesianOutput")){
    priorClass = createPriorDensity(prior, lower = lower, upper = upper, best = best)
  } else if(! "prior" %in% class(prior)){
    if(is.null(prior) & !is.null(lower) & !is.null(upper)){
      priorClass = createUniformPrior(lower = lower,upper = upper, best = best)
    }else if(is.null(prior) || class(prior) == "function"){
      priorClass = createPrior(density = prior, sampler = priorSampler, lower = lower, upper = upper, best)
    }else{
      stop("wrong input in prior")
    }
  } else{
    priorClass = prior 
    if("function" != class(prior)){
    if( !is.null(lower) & !is.null(prior$lower)) warning("Boundary values provided in prior and Bayesiansetup, the latter will be ignored")
    if( !is.null(upper) & !is.null(prior$upper)) warning("Boundary values provided in prior and Bayesiansetup, the latter will be ignored")
  }
  }
  
  if(! "likelihood" %in% class(likelihood)){
    if(class(likelihood) == "function"){
      likelihoodClass = createLikelihood(likelihood, parallel= parallel,parallelOptions = parallelOptions ,catchDuplicates = catchDuplicates)
    }else{
      stop("likelihood must be a function")
    }
  }else{
    likelihoodClass = likelihood
  }
  
  # TODO - macht das Sinn numPars hier zu definieren?
  if(is.null(priorSampler)) priorSampler <- function(x) return(NULL) # Avoids error in calculation of numPars
  
  if("function" != class(prior)){
    numPars = max(length(best), length(lower), length(upper), length(prior$lower), length(prior$upper), length(names),
                  length(priorSampler()))
  } else {
    numPars = max(length(best), length(lower), length(upper),  length(names),
                  length(priorSampler()))
    
  }
  
  
  if (is.null(names) & numPars > 0) names =  paste("par", 1:numPars)
  #else numPars = NULL

  
  posteriorClass = createPosterior(priorClass,likelihoodClass)

  out <- list(prior = priorClass, likelihood = likelihoodClass, posterior = posteriorClass, names = names, numPars = numPars, model = model, parallel = parallel)
  class(out) <- "BayesianSetup"
  
  return(out)
}
# 
# #' Generates initial sample TODO
# #' @param n TODO
# #' @param checkInf TODO
# #' @param overdispersed TODO
# #' @param maxIterations TODO
# #' @export 
# generateInitialSamples <- function(n, checkInf = T, overdispersed = F, maxIterations = 5){
#     if(is.null(sampler)) stop("sampling not implemented")
#     done = F
#     
#     stop("to implement")
#     
#     # check infinity of likelihood / create overdispersion
#     
#   }

#TODO: FH I wonder if we should keep this function option alive - seems better to me to explicitly do 
  # this with the createBayesianSetup 
  
#' Checks if an object is of class 'BayesianSetup'
#' @description Function used to assure that an object is of class 'BayesianSetup'. If you pass a function, it is coverted to an object of class 'BayesianSetup' (using \code{\link{createBayesianSetup}}) before it is returned.
#' @param bayesianSetup either object of class bayesianSetup or a log posterior function
#' @param parallel if bayesianSetup is a function, this will set the parallelization option for the class BayesianSetup that is created internally. If bayesianSetup is already a BayesianSetup, then this will check if parallel = T is requested but not supported by the BayesianSetup. This option is for internal use in the samplers
#' @note The recommended option to use this function in the samplers is to have parallel with default NULL in the samplers, so that checkBayesianSetup with a function will create a bayesianSetup without parallelization, while it will do nothing with an existing BayesianSetup. If the user sets parallelization, it will set the approriate parallelization for a function, and check in case of an existing BayesianSetup. The checkBayesianSetup call in the samplers should then be followed by a check for parallel = NULL in sampler, in which case paralell can be set from the BayesianSetup
#' @seealso \code{\link{createBayesianSetup}}
#' @export
checkBayesianSetup <- function(bayesianSetup, parallel = F){
  
  if(class(bayesianSetup) == "function"){
    if(is.null(parallel)) parallel = F
    bayesianSetup = createBayesianSetup(bayesianSetup, parallel = parallel)
  } 
  else if(class(bayesianSetup) == "BayesianSetup"){
    if(!is.null(parallel)) if(parallel == T & bayesianSetup$parallel == F) stop("parallel = T requested in sampler but BayesianSetup does not support parallelization. See help of BayesianSetup on how to enable parallelization")
  } 
  else stop("bayesianSetup must be class BayesianSetup or a function")
  
  return(bayesianSetup)
}


#' Function to close cluster in BayesianSetup
#' @description Function closes 
#' the parallel executer (if available)
#' @param  bayesianSetup object of class BayesianSetup
#' @export
stopParallel <- function(bayesianSetup){

  ## Stop cluster
  try(parallel::stopCluster(bayesianSetup$likelihood$cl), silent = TRUE)
  
  ## Remove object
 # pos <- -1
#  if(is.null(envir)) envir <- as.environment(pos)
  
 # .Internal(remove(deparse(substitute(bayesianSetup)), envir = envir, inherits = FALSE))
  
}




