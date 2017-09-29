
#' Convert coda::mcmc objects to BayesianTools::mcmcSampler
#' @description Function is used to make the plot and diagnostic functions 
#' available for coda::mcmc objects
#' @param sampler An object of class mcmc or mcmc.list 
#' @param names vector giving the parameter names (optional)
#' @param info matrix (or list with matrices for mcmc.list objects) with three coloumns containing log posterior, log likelihood and log prior of the sampler for each time step (optional; but see Details)
#' @param likelihood likelihood function used in the sampling (see Details)
#' @details The parameter 'likelihood' is optional for most functions but can be needed e.g for 
#' using the \code{\link{DIC}} function.
#' 
#' Also the parameter info is optional for most uses. However for some functions (e.g. \code{\link{MAP}})
#' the matrix or single coloumns (e.g. log posterior) are necessary for the diagnostics.
#' @export

convertCoda <- function(sampler, names = NULL, info = NULL, likelihood = NULL){

  likelihood <- list(density = likelihood)
  
  if(class(sampler) == "mcmc"){
    
    if(is.null(names)){
      names <- paste("Par",1:ncol(sampler))
    }
    setup <- list(names = names, numPars = ncol(sampler), likelihood = likelihood)
    
    if(is.null(info)) info <- matrix(NA, nrow = nrow(sampler), ncol = 3)
    out <- list(chain = cbind(sampler,info), setup = setup)
    class(out) = c("mcmcSampler", "bayesianOutput")
    
    
  }else{ if(class(sampler) == "mcmc.list"){
    
    if(is.null(names)){
      names <- paste("Par",1:ncol(sampler[[1]]))
    }
    setup <- list(names = names, numPars = ncol(sampler[[1]]), likelihood = likelihood)
    
    if(is.null(info)){
      info <- list()
      for(i in 1:length(sampler)) info[[i]] <- matrix(NA, nrow = nrow(sampler[[1]]), ncol = 3)
    } 
    
    chain <- list()
    for(i in 1:length(sampler)){
      chain[[i]] <- cbind(sampler[[i]], info[[i]])
    }
    class(chain) = "mcmc.list"
    out <- list(chain  = chain, setup = setup)
    class(out) = c("mcmcSampler", "bayesianOutput")
  }else stop("sampler must be of class 'coda::mcmc' or 'coda::mcmc.list'")
  }
  return(out)
  
}