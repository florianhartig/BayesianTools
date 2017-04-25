#' @export
marginalPlot <- function(x, ...) UseMethod("marginalPlot")

#' Plot MCMC marginals
#' @author Florian Hartig
#' @param mat object of class "bayesianOutput" or a matrix or data frame of variables 
#' @param thin thinning of the matrix to make things faster. Default is to thin to 5000 
#' @param scale should the results be scaled. Value can be either NULL (no scaling), T, or a matrix with upper / lower bounds as columns. If set to T, attempts to retrieve the scaling from the input object mat (requires that this is of class BayesianOutput)
#' @param best if provided, will draw points at the given values (to display true / default parameter values). Value can be either NULL (no drawing), a vector with values, or T, in which case the function will attempt to retrieve the values from a BayesianOutput
#' @param ... additional parameters to pass on to the \code{\link{getSample}}
#' @export
#' @seealso \code{\link{plotTimeSeries}} \cr
#'          \code{\link{tracePlot}} \cr
#'          \code{\link{correlationPlot}}
#' @example /inst/examples/plotMarginals.R
marginalPlot <- function(mat, thin = "auto", scale = NULL, best = NULL, ...){
  
 
  try(attachNamespace("sm"), silent = T) # this fixes dirty import in vioplot, to avoid error Error in do.call("sm.density", c(list(data, xlim = est.xlim), args)) : could not find function "sm.density"
  
  if(is.null(scale)) scale <- FALSE
  if(is.null(best))best <- FALSE
  
  
  if(inherits(mat,"bayesianOutput")){
    if(is.logical(scale)){
      if(scale == TRUE & !is.null(mat$setup$prior$lower) & !is.null(mat$setup$prior$upper)) scale = cbind(mat$setup$prior$lower, mat$setup$prior$upper)
    }  
    if(is.logical(best)){
      if(best == TRUE & !is.null(mat$setup$prior$best)) best = mat$setup$prior$best
    } 
    mat = getSample(mat, thin = thin)
  } 
  
  numPars = ncol(mat)
  names = colnames(mat)

  # TODO this is a hack to make the is.numeric(scale) below work for data.frame inputs, e.g. marginalPlot(out, scale = refPars[parSel, 2:3], best = refPars[parSel,1], start = 5000) . In general, the type conversion in this function should be cleaned up.
  if(is.data.frame(scale)) scale = as.matrix(scale)  
  
  if(is.numeric(scale)) {
    min = scale[,1]
    max = scale[,2]
    mat = scaleMatrix(mat, min, max)
    if(is.numeric(best))  best = scaleMatrix(best, min, max)
  }


  # TODO ... add names
  
  if(is.logical(scale)){
    if(scale == FALSE){
    main = "Marginal parameter uncertainty, unscaled"
    xlim = range(mat)
    }
  }else{
    main = "Marginal parameter uncertainty\n scaled to min/max values provided"
    xlim = c(0,1)
  }
  
  plot(NULL, ylim = c(0,numPars +1), type="n", yaxt="n", xlab="", ylab="", xlim = xlim, main = main)
  for (i in 1:numPars){
    vioplot::vioplot(mat[,i], at = i, add = T, col = "orangered", horizontal = T)
    axis(side = 2,at=i,labels = names[i],las=1)
  }
  if(is.numeric(best)) points(best,1:length(best), cex = 3, pch = 4, lwd = 2)
  

}





