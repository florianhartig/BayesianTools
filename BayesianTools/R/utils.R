#' Function to get the setup from a bayesianOutput
#' @param x bayesianOutput
#' @return bayesianSetup
#' @author Tankred Ott
#' @keywords internal
getSetup <- function(x) {
  classes <- class(x)
  if (any(c('mcmcSampler', 'smcSampler') %in% classes)) x$setup
  else if (any(c('mcmcSamplerList', 'smcSamplerList') %in% classes)) x[[1]]$setup
  else stop('Can not get setup from x')
}

#' Function to thin matrices
#' @param mat matrix to thin
#' @param thin thinning parameter
#' @return thinned matrix
#' @keywords internal
thinMatrix <- function(mat, thin = "auto"){
  if (thin == "auto"){
    thin = max(floor(nrow(mat) / 5000),1)
  }
  
  if (! is.null(thin) & ! thin == F){
    sel = seq(1,dim(mat)[1], by = thin )
    mat = mat[sel,]
  }
  return(mat)
}


#' Function to scale matrices
#' @param mat matrix to scale
#' @param min minimum value
#' @param max maximum value
#' @return sclaed matrix
#' @keywords internal
scaleMatrix <- function(mat, min, max){
  if(class(mat)[1] %in% c("matrix", "data.frame")){
    for(i in 1:ncol(mat)){
      mat[,i] <- (mat[,i] - min[i]) / (max[i] - min[i])
    }    
  }else if (is.vector(mat)){
    mat = (mat - min) / (max - min)
  }else stop("wrong class")
  return(mat)
}


#' Function to calculate the metropolis ratio
#' @author Florian Hartig
#' @param LP2 log posterior old position
#' @param LP1 log posterior of proposal
#' @param tempering value for tempering
#' @keywords internal
metropolisRatio <- function(LP2, LP1, tempering = 1){
  # this catches two -Inf cases / I wonder if we should throw a warning in this case
  if( is.na(LP2 - LP1)) out = -Inf
  else out =   min(exp( (LP2 - LP1) / tempering), 1)
  return(out)
} 


#' getPanels
#' 
#' Calculates the argument x for par(mfrow = x) for a desired number of panels
#' 
#' @author Florian Hartig
#' @param x the desired number of panels 
#' @export
getPanels <- function(x){
  if (x <= 0) stop("number can't be < 1")
  
  lower = floor(sqrt(x))
  upper = ceiling(sqrt(x))
  
  if (lower == upper) return(c(lower, lower))
  else{
    if (lower * upper >= x) return(c(lower, upper))
    else return(c(upper, upper))    
  }
}  

#' Gets n equally spaced samples (rows) from a matrix or vector
#' @author Tankred Ott
#' @param x matrix or vector
#' @param numSamples number of samples (rows) to be drawn
#' @details Gets n equally spaced samples (rows) from a matrix and returns a new matrix (or vector) containing those samples
#' @keywords internal
sampleEquallySpaced <- function(x, numSamples) {
  # wrong input: x is neither vector nor matrix
  if (!is.matrix(x) && !is.vector(x)) {
    stop("Expected matrix or vector for x!")
  }
  # wrong input: numSamples is not single numeric value
  if (!is.vector(numSamples) || !is.numeric(numSamples) || length(numSamples) > 1) {
    stop("Expected a single numeric value for numSamples!")
  }
  
  len = 0
  if (is.matrix(x)) {
    len = nrow(x)
  } else {
    len = length(x)
  }
  
  if (len == 1) {
    return(x)
  }
  
  # wrong input: numSamples > total number of samples
  if (numSamples > len) {
    numSamples = len
    warning("numSamples is greater than the total number of samples! All samples were selected.")
  # wrong input: numsaples 0 or negative
  } else if (numSamples < 1) { 
    numSamples = 1;
    warning("numSamples is less than 1! Only the first sample was selected.")
  }
  
  sel <- seq(1, len, len = numSamples) # RB: what does len = numSample do? pass as argument to seq? -> seq doesnt use it
  if (is.matrix(x)) {
    out <- x[sel, , drop=F]
    # if x has only a single col, x[sel,] is a vector and needs to be converted ## RB: add drop=F before and remove if statement to avoid loss of colnames
#    if(!is.matrix(out)) {
#      out <- matrix(out, ncol = ncol(x))
#    }
  } else {
    out <- x[sel]
  }
  
  return(out)
}

#' Checks if thin is conistent with nTotalSamples samples and if not corrects it.
#' @author Tankred Ott
#' @param nTotalSamples total number of rows/samples 
#' @param thin thinning
#' @param autoThinFraction fraction of the data that will be sampled when thin is set to "auto". E.g. 0.5 means thin will be nTotalSamples * 0.5. The resulting thin value is rounded down to the next integer.
#' @details Checks if the thin argument is consistent with the data consisting of nTotalSamples samples/rows and corrects thin if not.
#' @author Tankred Ott
# #' @export
#' @keywords internal
correctThin <- function(nTotalSamples, thin, autoThinFraction = 0.001) {
  if (autoThinFraction > 1 || autoThinFraction <= 0) {
    stop("autoThinFraction must be greater than 0 and less than 1!")
  }
  
  if (thin == "auto"){
    thin = max(floor(nTotalSamples * autoThinFraction), 1)
  } else if (is.null(thin) || thin == F || thin < 1 || is.nan(thin)) {
    thin = 1
  } else if (thin > nTotalSamples) {
    warning("thin is greater than the total number of samples! Only the first sample/row was selected.")
    thin = nTotalSamples
  }
  
  return(thin)
}

#' @title Rescale
#' @description Rescales values in the interval "from" (lower, upper) to the new interval "to" (lower, upper).
#' @param x vector of values tp be scaled
#' @param from vector of length 2, original interval (lower, upper)
#' @param to vector of length 2, target interval (lower, upper)
#' 
#' @keywords internal
#' @author Tankred Ott
rescale <- function (x, from, to) {
  # scale x from 0 to 1
  x <- (x - from[1]) / (from[2] - from[1])
  # scale to new interval
  return(x * (to[2] - to[1]) + to[1])
}


