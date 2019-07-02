#' Creates a standardized collection of prior, likelihood and posterior functions, including error checks etc.
#' @author Florian Hartig, Tankred Ott
#' @param likelihood log likelihood density function
#' @param prior either a prior class (see \code{\link{createPrior}}) or a log prior density function
#' @param priorSampler if a prior density (and not a prior class) is provided to prior, the optional prior sampling function can be provided here
#' @param lower vector with lower prior limits
#' @param upper vector with upper prior limits
#' @param best vector with best prior values
#' @param names optional vector with parameter names
#' @param parallel parallelization option. Default is F. Other options are T, or "external". See details.
#' @param parallelOptions list containing three lists. First "packages" determines the R packages necessary to run the likelihood function. Second "variables" the objects in the global envirnment needed to run the likelihood function and third "dlls" the DLLs needed to run the likelihood function (see Details and Examples). 
#' @param catchDuplicates Logical, determines whether unique parameter combinations should only be evaluated once. Only used when the likelihood accepts a matrix with parameter as columns. 
#' @param plotLower vector with lower limits for plotting
#' @param plotUpper vector with upper limits for plotting
#' @param plotBest vector with best values for plotting
#' @details If prior is of class prior (e.g. create with \code{\link{createPrior}}), priorSampler, lower, upper and best will be ignored.\cr If prior is a function (log prior density), priorSampler (custom sampler), or lower/upper (uniform sampler) is required.\cr If prior is NULL, and lower and upper are passed, a uniform prior (see \code{\link{createUniformPrior}}) will be created with boundaries lower and upper.
#' 
#' For parallelization, option parallel = T means that an automatic parallelization via a standard R socket cluster is attempted. By default, a copy of your workspace, including DLLs and objects are exporte to the cluster workers. Because this can be very inefficient, you can explictly specify the packages, objects and DLLs that are to be exported via parallelOptions.
#' 
#' Using parallel = T requires that the function to be parallelized is well encapsulate, i.e. can run on a shared memory / shared hard disk machine in parallel without interfering with each other. For some functions and programs, this is not the case, so that a custom-programmed parallelization is required. 
#' 
#' In this case, and only in this case, you should specify parallel = "external". In this case, it is assumed that the likelihood is programmed such that it accepts a matrix with parameters as columns and the different model runs as rows. It is then up to the user if and how to parallelize this function. This option gives most flexibility to the user, in particular for complicated parallel architecture or shared memory problems.
#'
#' For more details on parallelization, make sure to read the vignette (run vignette("BayesianTools", package="BayesianTools"))
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
                                parallelOptions = list(variables = "all", packages = "all", dlls = NULL), 
                                catchDuplicates = FALSE,
                                plotLower = NULL,
                                plotUpper = NULL,
                                plotBest = NULL
){
  
  # TODO implement parameter "model" (function that makes predictions from the model)
  model <- NULL
  
  
  # INPUTS CHECKS
  if(is.null(upper) && is.null(lower) && is.null(prior)) stop("Either boundaries or prior density and prior sampler must be provided.")
  # if(!is.null(lower) || !is.null(upper) || !is.null(best)) print("DEPRECATED: lower/upper/best arguments for createBayesianSetup are deprecated and will be removed in a future update. Pass those arguments in the info parameter instead or use createUnformPrior.")
  if(("prior" %in% class(prior)) && (!is.null(lower) || !is.null(upper))) warning("Prior object and boundary values provided to createBayesiansetup, the latter will be ignored")
  if(("prior" %in% class(prior)) && (!is.null(priorSampler))) warning("Prior object and priorSampler provided to createBayesiansetup, the latter will be ignored")
  
  if(is.null(parallelOptions)) parallelOptions <- list(variables = "all", packages = "all", dlls = "all")
  
  
  # PRIOR CHECKS
  priorClass = NULL
  if ("prior" %in% class(prior)) {
    priorClass = prior
    
  } else if (inherits(prior,"bayesianOutput")) {
    priorClass = createPriorDensity(prior)
    
  } else if ("function" %in% class(prior)) {
    if ("function" %in% class(priorSampler)) priorClass = createPrior(prior, priorSampler)
    else if (!is.null(lower) && !is.null(upper)) priorClass = createPrior(prior, lower=lower, upper=upper, best=best)
    else stop("If prior is a function, priorSampler or lower/upper is required")
    
  } else if (is.null(prior)) {
    # TODO: deprecate this
    # checks for NULL for lower/upper are already done at begin of function
    priorClass = createUniformPrior(lower = lower, upper = upper, best = best)
    
  } else stop("wrong input for prior")
  
  
  # LIKELIHOOD CHECKS
  if ("likelihood" %in% class(likelihood)) {
    likelihoodClass = likelihood
  } else if ("function" %in% class(likelihood)) {
    likelihoodClass = createLikelihood(likelihood, parallel = parallel, parallelOptions = parallelOptions, catchDuplicates = catchDuplicates)
  } else {
    stop("likelihood must be an object of class likelihood or a function")
  }
  pwLikelihood = likelihoodClass$pwLikelihood
  
  # GET NUMBER OF PARAMETERS
  numPars = length(priorClass$sampler())
  
  # CREATE POSTERIOR
  posteriorClass = createPosterior(priorClass,likelihoodClass)
  
  # CHECK FOR PLOTTING PARAMETERS
  if (is.null(plotLower)) plotLower <- priorClass$lower
  if (is.null(plotUpper)) plotUpper <- priorClass$upper
  if (is.null(plotBest)) plotBest <- priorClass$best
  
  if (is.null(plotLower) | is.null(plotUpper) | is.null(plotBest))
    print("Info is missing upper/lower/best. This can cause plotting and sensitivity analysis functions to fail. If you want to use those functions provide (plot)upper/lower/best either in createBayesianSetup or prior")
  
  # CHECK NAMES
  if (is.null(names)) {
    if (!is.null(priorClass$parNames)) names = priorClass$parNames
    else if (!is.null(likelihoodClass$parNames)) names = likelihoodClass$parNames
    else if (numPars > 0) names = paste("par", 1:numPars)
  }
  
  # CONSTRUCT OUTPUT
  info <- list(priorLower = priorClass$lower, priorUpper = priorClass$upper, priorBest = priorClass$best,
               plotLower = plotLower, plotUpper = plotUpper, plotBest = plotBest,
               parNames = names, numPars = numPars)
  out <- list(prior = priorClass, likelihood = likelihoodClass, posterior = posteriorClass,
              names = names, numPars = numPars, model = model, parallel = parallel, pwLikelihood = pwLikelihood, info = info)
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
#' @author Florian Hartig
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
#' @author Stefan Paul
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


#' @author Maximilian Pichler
#' @export

print.BayesianSetup <- function(x, ...){
  cat('BayesianSetup: \n\n')
  
  bayesianSetup = x
  info = c( "priorLower", "priorUpper", "plotLower", "plotUpper")
  parInfo = data.frame(matrix(NA, ncol = 4, nrow = bayesianSetup$info$numPars))
  colnames(parInfo) = info
  rownames(parInfo) = bayesianSetup$info$parNames
  for(i in 1:4) if(!is.null(bayesianSetup$info[[info[i]]])) parInfo[,i] <- bayesianSetup$info[[info[i]]]
  print(parInfo)
  
}
