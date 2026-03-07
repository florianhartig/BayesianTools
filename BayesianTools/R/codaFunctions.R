#' Function to combine chains
#' 
#' @param x a list of MCMC chains
#' @param merge logical value (T or F) to determine if chains should be merged 
#' @return combined chains
#' @note to combine several chains to a single McmcSamplerList, see [createMcmcSamplerList]
#' @keywords internal
combineChains <- function(x, merge = T){
  
  if(merge == T){
    temp1 = as.matrix(x[[1]])
    
    names = colnames(temp1)
    
    sel = seq(1, by = length(x), len = nrow(temp1) )
    
    out = matrix(NA, nrow = length(x) * nrow(temp1), ncol = ncol(temp1))
    out[sel, ] = temp1
    if (length(x) > 1){
      for (i in 2:length(x)){
        out[sel+i-1, ] = as.matrix(x[[i]])
      }   
    } 
    
    colnames(out) = names
    
  } else{
    
    out = as.matrix(x[[1]])
    if (length(x) > 1){
      for (i in 2:length(x)){
        out = rbind(out, as.matrix(x[[i]]))
      }   
    } 
  }
  
  return(out)
}



#' Helper function to change an object to a coda mcmc class,
#'
#' @param chain mcmc Chain
#' @param start for MCMC samplers, the initial parameter vector in the chain. For SMC samplers, the initial particle swarm
#' @param end for MCMC samplers, the last parameter vector in the chain. For SMC samplers, end final particle swarm.
#' @param thin thinning parameter
#' @return an object of class coda::mcmc
#' @details   Very similar to coda::mcmc but with less overhead
#' @keywords internal

makeObjectClassCodaMCMC <- function (chain, start = 1, end = numeric(0), thin = 1){
  attr(chain, "mcpar") <- c(start, end, thin)
  attr(chain, "class") <- "mcmc"
  chain
}


#' Convert coda::mcmc objects to BayesianTools::mcmcSampler
#' @description Function to support plotting and diagnostic functions for coda::mcmc objects.
#' @param sampler an object of class mcmc or mcmc.list 
#' @param names a vector with parameter names (optional)
#' @param info a matrix (or list with matrices for mcmc.list objects) with three columns containing log posterior, log likelihood and log prior of the sampler for each time step (optional; but see Details)
#' @param likelihood likelihood function used for sampling (see Details)
#' @details The parameter 'likelihood' is optional for most functions but can be needed e.g for  \code{\link{DIC}} function.
#' Also, the parameter information is typically optional for most uses. However, for certain functions (e.g. \code{\link{MAP}}), the matrix or single columns (e.g. log posterior) are necessary for diagnostics.
#' @export

convertCoda <- function(sampler, names = NULL, info = NULL, likelihood = NULL){
  
  likelihood <- list(density = likelihood)
  
  if(inherits(sampler, "mcmc")){
    
    if(is.null(names)){
      names <- paste("Par",1:ncol(sampler))
    }
    setup <- list(names = names, numPars = ncol(sampler), likelihood = likelihood)
    
    if(is.null(info)) info <- matrix(NA, nrow = nrow(sampler), ncol = 3)
    out <- list(chain = cbind(sampler,info), setup = setup)
    class(out) = c("mcmcSampler", "bayesianOutput")
    
    
  }else{ if(inherits(sampler, "mcmc.list")){
    
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

