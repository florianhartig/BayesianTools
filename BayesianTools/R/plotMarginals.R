#' @export
marginalPlot <- function(x, ...) UseMethod("marginalPlot")

#' Plot MCMC marginals
#' @author Florian Hartig
#' @param mat object of class "bayesianOutput" or a matrix or data frame of variables 
#' @param thin thinning of the matrix to make things faster. Default is to thin to 5000 
#' @param scale should the results be scaled. Value can be either NULL (no scaling), T, or a matrix with upper / lower bounds as columns. If set to T, attempts to retrieve the scaling from the input object mat (requires that this is of class BayesianOutput)
#' @param best if provided, will draw points at the given values (to display true / default parameter values). Value can be either NULL (no drawing), a vector with values, or T, in which case the function will attempt to retrieve the values from a BayesianOutput
#' @param histogram Logical, determining whether a violin plot or a histogram should be plotted
#' @param plotPrior Logical, determining whether the prior should be plotted in addition to the posteriors. Only applicable if mat is an object of class "bayesianOutput"
#' @param nDrawsPrior Integer, number of draws from the prior, when plotPrior is active 
#' @param ... additional parameters to pass on to the \code{\link{getSample}}
#' @export
#' @references Internally, this function uses an adapted version of the function vioplot from the vioplot R package (Copyright (c) 2004, Daniel Adler)
#' @seealso \code{\link{plotTimeSeries}} \cr
#'          \code{\link{tracePlot}} \cr
#'          \code{\link{correlationPlot}}
#' @example /inst/examples/plotMarginals.R
marginalPlot <- function(mat, thin = "auto", scale = NULL, best = NULL, histogram = FALSE, plotPrior = FALSE, nDrawsPrior = 1000, ...){
  priorMat <- NULL
  
  if (plotPrior == TRUE) {
    if(("mcmcSamplerList" %in% class(mat)) || ("smcSamplerList" %in% class(mat))) {
      priorMat <- mat[[1]]$setup$prior$sampler(nDrawsPrior)
    } else if (("mcmcSampler" %in% class(mat)) || ("smcSampler" %in% class(mat))) {
      priorMat <- mat$setup$prior$sampler(nDrawsPrior)
    } else {
      warning("Parameter 'mat' is not of class 'bayesianOutput', set plotPrior to FALSE.")
      plotPrior <- FALSE
    }
  }
  
  if(is.null(scale)) scale <- FALSE
  if(is.null(best)) best <- FALSE
  
  
  if(inherits(mat,"bayesianOutput")){
    if(("mcmcSamplerList" %in% class(mat)) || ("smcSamplerList" %in% class(mat))) {
      lower <- mat[[1]]$setup$prior$lower
      upper <- mat[[1]]$setup$prior$upper
      if (best) best <- mat[[1]]$setup$prior$best
    } else if (("mcmcSampler" %in% class(mat)) || ("smcSampler" %in% class(mat))) {
      lower <- mat$setup$prior$lower
      upper <- mat$setup$prior$upper
      if (best) best <- mat$setup$prior$best
    }
    if(is.logical(scale)){
      if(scale == TRUE & !is.null(lower) & !is.null(upper)) scale = cbind(lower, upper)
    }  
    if(is.logical(best)){
      if(best == TRUE & !is.null(best)) best = best
    } 
    mat = getSample(mat, thin = thin)
  } else if ("matrix" %in% class(mat)) {
    if(scale==TRUE) scale <- t(apply(x, 2, range))
  }
  
  numPars = ncol(mat)
  names = colnames(mat)

  # TODO this is a hack to make the is.numeric(scale) below work for data.frame inputs, e.g. marginalPlot(out, scale = refPars[parSel, 2:3], best = refPars[parSel,1], start = 5000) . In general, the type conversion in this function should be cleaned up.
  if(is.data.frame(scale)) scale = as.matrix(scale)  
  
  if(is.numeric(scale)) {
    min = scale[,1]
    max = scale[,2]
    mat = BayesianTools:::scaleMatrix(mat, min, max)
    if (plotPrior) priorMat <- BayesianTools:::scaleMatrix(priorMat, min, max)
    if(is.numeric(best)) best <- BayesianTools:::scaleMatrix(best, min, max)
  }


  # TODO ... add names
  xlim = if (plotPrior) range(mat, priorMat) else range(mat)
  if(is.logical(scale)){
    if(scale == FALSE){
      main = "Marginal parameter uncertainty, unscaled"
      xlim = if (plotPrior) range(mat, priorMat) else range(mat)
    }
  }else{
    main = "Marginal parameter uncertainty\n scaled to min/max values provided"
    xlim = c(0,1)
  }
  
  plot(NULL, ylim = c(0,numPars +1), type="n", yaxt="n", xlab="", ylab="", xlim = xlim, main = main)
  for (i in 1:numPars){
    if (plotPrior) {
      vioplot(priorMat[,i], at = i, add = T, col = "#4682B4A0", horizontal = T, style = "bottom")
      vioplot(mat[,i], at = i, add = T, col = "orangered", horizontal = T, style = "top")
    } else {
      vioplot(mat[,i], at = i, add = T, col = "orangered", horizontal = T)
    }
  
    axis(side = 2,at=i,labels = names[i],las=1)
  }
  if(is.numeric(best)) points(best,1:length(best), cex = 3, pch = 4, lwd = 2)
  if(plotPrior) legend("topright", c("posterior", "prior"), col = c("orangered", "#4682B4A0"), pch = 15, cex = 0.8, bty = "n")

}


# Adapted from the R package vioplot
# Copyright (c) 2004, Daniel Adler
vioplot <- function(x,...,range=1.5,h=NULL,ylim=NULL,names=NULL, horizontal=FALSE, 
                    col="magenta", border="black", lty=1, lwd=1, rectCol="black", colMed="white", pchMed=19, at, add=FALSE, wex=1, 
                    drawRect=TRUE, style="both")
{
  # process multiple datas
  datas <- list(x,...)
  n <- length(datas)
  
  if(missing(at)) at <- 1:n
  
  # pass 1
  #
  # - calculate base range
  # - estimate density
  #
  
  # setup parameters for density estimation
  upper  <- vector(mode="numeric",length=n)
  lower  <- vector(mode="numeric",length=n) 
  q1     <- vector(mode="numeric",length=n)
  q3     <- vector(mode="numeric",length=n)
  med    <- vector(mode="numeric",length=n)
  base   <- vector(mode="list",length=n)
  height <- vector(mode="list",length=n)
  baserange <- c(Inf,-Inf)
  
  # global args for sm.density function-call   
  args <- list(display="none")
  
  if (!(is.null(h)))
    args <- c(args, h=h)
  
  
  for(i in 1:n) {
    
    data<-datas[[i]]
    
    # calculate plot parameters
    #   1- and 3-quantile, median, IQR, upper- and lower-adjacent
    
    data.min <- min(data)
    data.max <- max(data)
    q1[i]<-quantile(data,0.25)
    q3[i]<-quantile(data,0.75)
    med[i]<-median(data)
    iqd <- q3[i]-q1[i]
    upper[i] <- min( q3[i] + range*iqd, data.max )
    lower[i] <- max( q1[i] - range*iqd, data.min )
    
    
    #   strategy:
    #       xmin = min(lower, data.min))
    #       ymax = max(upper, data.max))
    #
    
    est.xlim <- c( min(lower[i], data.min), max(upper[i], data.max) ) 
    
    # estimate density curve
    
    smout <- do.call(sm::sm.density, c( list(data, xlim=est.xlim), args ) )
    
    
    # calculate stretch factor
    #
    #  the plots density heights is defined in range 0.0 ... 0.5 
    #  we scale maximum estimated point to 0.4 per data
    #
    
    hscale <- 0.4/max(smout$estimate) * wex
    
    
    # add density curve x,y pair to lists
    
    base[[i]]   <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    
    
    # calculate min,max base ranges
    
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1],t[1])
    baserange[2] <- max(baserange[2],t[2])
    
  }
  
  # pass 2
  #
  # - plot graphics
  
  # setup parameters for plot
  
  if(!add){
    xlim <- if(n==1) 
      at + c(-.5, .5)
    else 
      range(at) + min(diff(at))/2 * c(-1,1)
    
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  } else {
    label <- names
  }
  
  boxwidth <- 0.05 * wex
  
  
  # setup plot
  
  if(!add)
    plot.new()
  if(!horizontal) {
    if(!add){
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1,at = at, label=label )
    }  
    
    box()
    
    for(i in 1:n) {
      
      # plot left/right density curve
      
      polygon( c(at[i]-height[[i]], rev(at[i]+height[[i]])), 
               c(base[[i]], rev(base[[i]])),
               col = col, border=border, lty=lty, lwd=lwd)
      
      
      if(drawRect){
        # plot IQR
        lines( at[c( i, i)], c(lower[i], upper[i]) ,lwd=lwd, lty=lty)
        
        # plot 50% KI box
        rect( at[i]-boxwidth/2, q1[i], at[i]+boxwidth/2, q3[i], col=rectCol)
        
        # plot median point
        points( at[i], med[i], pch=pchMed, col=colMed )
      }
    }
    
  }
  else {
    if(!add){
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2,at = at, label=label )
    }
    
    box()
    for(i in 1:n) {
      
      # plot left/right density curve
      yUp <- yDown <- NA
      if (style=="top") {
        yUp <- rev(at[i]+height[[i]])
        yDown <- rep(at[i], length(height[[i]]))
      } else if (style=="bottom") {
        yUp <- rep(at[i], length(height[[i]]))
        yDown <- at[i]-height[[i]]
      } else {
        if (style != "both") {
          warning("style not recognized. Set style to 'both'.")
          style <- "both"  
        }
        yUp <- rev(at[i]+height[[i]])
        yDown <- at[i]-height[[i]]
      }
      polygon( c(base[[i]], rev(base[[i]])),
               # c(at[i]-height[[i]], rev(at[i]+height[[i]])),
               c(yDown, yUp),
               col = col, border=border, lty=lty, lwd=lwd)
      
      
      if(drawRect){
        # plot IQR
        
        if(style == "both"){
          boxTop <- at[i] + boxwidth/2
          boxBottom <- at[i] - boxwidth/2
          lineY <- at[c(i,i)]
        } else if (style == "bottom") {
          boxTop <- at[i]
          boxBottom <- at[i] - boxwidth
          lineY <- at[c(i,i)] - boxwidth / 2
        } else if (style == "top") {
          boxTop <- at[i] + boxwidth
          boxBottom <- at[i]
          lineY <- at[c(i,i)] + boxwidth / 2
        }
        
        # lines( c(lower[i], upper[i]), lineY ,lwd=lwd, lty=lty)
        
        # plot 50% KI box
        rect( q1[i], boxBottom, q3[i], boxTop,  col=rectCol)
        
        
        # plot median point
        points( med[i], lineY[1], pch=pchMed, col=colMed)
      }
    }
    
    
  }
  
  
  invisible (list( upper=upper, lower=lower, median=med, q1=q1, q3=q3))
}


