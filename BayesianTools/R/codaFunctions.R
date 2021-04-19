#' Function to combine chains
#' 
#' @param x a list of MCMC chains
#' @param merge logical determines whether chains should be merged
#' @return combined chains
#' 
#' @note to combine several chains to a single McmcSamplerList, see \code{\link{createMcmcSamplerList}}
#' 
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
#' @param start for mcmc samplers start value in the chain. For SMC samplers, start particle
#' @param end for mcmc samplers end value in the chain. For SMC samplers, end particle
#' @param thin thinning parameter
#' @return object of class coda::mcmc
#' @details   Very similar to coda::mcmc but with less overhead
#' @keywords internal
makeObjectClassCodaMCMC <- function (chain, start = 1, end = numeric(0), thin = 1){
  attr(chain, "mcpar") <- c(start, end, thin)
  attr(chain, "class") <- "mcmc"
  chain
}


