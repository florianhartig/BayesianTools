#' Function to thin matrices
#' @param mat matrix to thin
#' @param thin thinning parameter
#' @return thinned matrix

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


#' Funktion to calculate the metropolis ratio
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
  


