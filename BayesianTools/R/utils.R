# #' Function to thin matrices
# #' @param mat matrix to thin
# #' @param thin thinning parameter
# #' @return thinned matrix

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


# #' Function to scale matrices
# #' @param mat matrix to scale
# #' @param min minimum value
# #' @param max maximum value
# #' @return sclaed matrix

scaleMatrix <- function(mat, min, max){
  if(class(mat) %in% c("matrix", "data.frame")){
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
#' @export
metropolisRatio <- function(LP2, LP1, tempering = 1){
  # this catches two -Inf cases / I wonder if we should throw a warning in this case
  if( is.na(LP2 - LP1)) out = -Inf
  else out =   min(exp( (LP2 - LP1) / tempering), 1)
  return(out)
} 


#' Calculates the panel combination for par(mfrow = )
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

#' Gets n equally spaced samples (rows) from a matrix
#' @author Tankred Ott
#' @param x matrix or vector
#' @param numSamples number of samples (rows) to be drawn
#' @details Gets n equally spaced samples (rows) from a matrix and returns a new matrix (or vector) containing those samples
#' @export
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
    if (len == 1) {
      return(x)
    }
  } else {
    len = length(x)
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
  
  sel <- seq(1, len, len = numSamples)
  if (is.matrix(x)) {
    out <- x[sel,]
    # if x is only a single row x[sel,] is a vector and needs to
    # be converted
    if(is.matrix(out)) {
      return(out)
    } else {
      return(matrix(out, byrow = FALSE))
    }
  } else {
    return(x[sel])
  }
}

#' Checks if thin is conistent with nTotalSamples samples and if not corrects it.
#' @author Tankred Ott
#' @param nTotalSamples total number of rows/samples 
#' @param thin thinning
#' @param autoThinFraction fraction of the data that will be sampled when thin is set to "auto". E.g. 0.5 means thin will be nTotalSamples * 0.5.
#' @details Checks if the thin argument is consistent with the data consisting of nTotalSamples samples/rows and corrects thin if not.
#' @author Tankred Ott
#' @export
correctThin <- function(nTotalSamples, thin, autoThinFraction = 0.01) {
  if (autoThinFraction > 1 || autoThinFraction <= 0) {
    stop("autoThinFraction must be greater than 0 and less than 1!")
  }
  
  if (thin == "auto"){
    thin = max(floor(nTotalSamples * 0.1), 1)
  } else if (is.null(thin) || thin == F || thin < 1 || is.nan(thin)) {
    thin = 1
  } else if (thin > nTotalSamples) {
    warning("thin is greater than the total number of samples! Only the first sample/row was selected.")
    thin = nTotalSamples
  }
  
  return(thin)
}
