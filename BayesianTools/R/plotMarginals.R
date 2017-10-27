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
#' @param priorTop Logical, determining whether the prior should be plotted top (TRUE) or bottom (FALSE)
#' @param nDrawsPrior Integer, number of draws from the prior, when plotPrior is active 
#' @param breaks Integer, number of histogram breaks if histogram is set to TRUE
#' @param res resolution parameter for violinPlot, determining how many descrete points should be used for the density kernel.
#' @param ... additional parameters to pass on to the \code{\link{getSample}}
#' @export
#' @references 
#'          \code{\link{tracePlot}} \cr
#'          \code{\link{correlationPlot}}
#' @example /inst/examples/plotMarginals.R
marginalPlot <- function(mat, thin = "auto", scale = NULL, best = NULL, histogram = FALSE, plotPrior = FALSE, priorTop = FALSE, nDrawsPrior = 1000, breaks=15, res=500,...){
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
    if(scale==TRUE) scale <- t(apply(mat, 2, range))
  }
  
  numPars = ncol(mat)
  names = colnames(mat)

  # TODO this is a hack to make the is.numeric(scale) below work for data.frame inputs, e.g. marginalPlot(out, scale = refPars[parSel, 2:3], best = refPars[parSel,1], start = 5000) . In general, the type conversion in this function should be cleaned up.
  if(is.data.frame(scale)) scale = as.matrix(scale)  
  
  if(is.numeric(scale)) {
    min = scale[,1]
    max = scale[,2]
    mat = scaleMatrix(mat, min, max)
    if (plotPrior) priorMat <- scaleMatrix(priorMat, min, max)
    if(is.numeric(best)) best <- scaleMatrix(best, min, max)
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
      if (histogram == TRUE) {
        histMarginal(posterior = mat[,i], prior = priorMat[,i], at = i, col = c("orangered", "#4682B4A0"), breaks = breaks)
      } else {
        priorPos <- posteriorPos <- NULL
        if (priorTop == TRUE) {
          priorPos <- c("above", "top")
          posteriorPos <- c("below", "bottom")
        } else {
          priorPos <- c("below", "bottom")
          posteriorPos <- c("above", "top")
        }
        
        violinPlot(mat[,i], at = i, .range = 0.475, add = T, col = "orangered", relToAt = posteriorPos[1], which = posteriorPos[2], res=res)
        violinPlot(priorMat[,i], at = i, .range = 0.475, add = T, col = "#4682B4A0", relToAt = priorPos[1], which = priorPos[2], res=res)
      }
    } else {
      if (histogram == TRUE) {
        histMarginal(posterior = mat[,i], prior = NULL, at = i, col = "orangered", breaks = breaks)
      } else {
        violinPlot(mat[,i], at = i, .range = 0.95, add = T, col = "orangered", relToAt = "centered", which = "both", res=res)
      }
    }
  
    axis(side = 2,at=i,labels = names[i],las=1)
  }
  if(is.numeric(best)) points(best,1:length(best), cex = 3, pch = 4, lwd = 2)
  
  if(plotPrior && histogram) legend("bottomright", c("posterior", "prior"), col = c("orangered", "#4682B4A0"), pch = 15, cex = 0.8, bty = "n")
  if(plotPrior && !histogram) legend("topright", c("posterior", "prior"), col = c("orangered", "#4682B4A0"), pch = 15, cex = 0.8, bty = "n")

}

#' @author Tankred Ott
#' @title Violin Plot
#' @description Function to plot classic violin plots, as well as "half violin plots" (density plots).
#' @param x vector of values to plot
#' @param at position of the plot when add is TRUE
#' @param .range maximum height if horizontal or width if vertical of the plot when add is set to FALSE
#' @param add logical, determining whether the plot should be added to an existing plot window
#' @param horizontal logical, determining whether the plot should be horizontal, if FALSE the plot will be vertical
#' @param which string, either "both" for a classic violing plot, or "top" or "bottom" to plot only the half violin
#' @param relToAt string, one of "centered", "above", "below", "left", "right". Determining the relative position to at.
#' @param plotQBox logical, determining whether the quantile box should be plotted
#' @param plotMed logical, determining whether the median should be plotted
#' @param col color of the violin
#' @param border color of the border of the violin
#' @param colQBox color of quantile box
#' @param borderQBox color of the border of the quantile box
#' @param colMed color of the median point
#' @param pchMed pch for median point
#' @param res "resolution" of the violin. Determining how many descrete points should be used to calculate the density kernel.
violinPlot <- function (x, at, .range = 1, add = FALSE, horizontal = TRUE, which = "both", relToAt = "above", plotQBox = TRUE, plotMed = TRUE,
                        col = "orangered", border = "black", colQBox = "black", borderQBox = "black", colMed = "white", pchMed = 19, res = 500) {
  
  q1 <- quantile(x, probs = 0.25)
  q3 <- quantile(x, probs = 0.75)
  med <- median(x)
  medX <- NULL
  medY <- NULL
  
  minX <- min(x)
  maxX <- max(x)
  dens <- density(x, from = minX, to = maxX, n = res)
  
  xVals <- NULL
  yVals <- NULL
  qBoxHeight <- NULL
  qBoxPoints <- vector(mode = "numeric", length = 4)
  qBoxMod <- 1 / 50 # maybe let the user decide how large the qBox should be?
  
  # if add is FALSE create empty plot
  if (add == FALSE) {
    if (horizontal == TRUE) plot(NULL, xlim = c(minX, maxX), ylim = c(0,1), xlab = "", ylab = "")
    else plot(NULL, xlim = c(0,1), ylim = c(minX, maxX), xlab = "", ylab = "")
    
    if (which == "both") {
      at <- 0.5
      .range <- 1
    } else if (which == "top") {
      at <- 0
      .range <- 1
    } else if (which == "bottom") {
      at <- 1
      .range <- 1
    }
  }
  
  # assign points to be plotted
  if (which == "both") {
    .at <- NULL
    if (relToAt == "above" || relToAt == "left") .at <- at + .range / 2
    else if (relToAt == "centered") .at <- at
    else if (relToAt == "below" || relToAt == "right") .at <- at - .range / 2
    qBoxHeight <- .range * qBoxMod
    
    yVals <- rescale(dens$y, c(0, max(dens$y)), c(0, .range / 2))
    yVals <- c(.at + yVals, .at - rev(yVals))
    xVals <- c(dens$x, rev(dens$x))
    
    qBoxPoints <- c(q1, .at - qBoxHeight, q3, .at + qBoxHeight)
    
    medX <- med
    medY <- .at
    
  } else if (which == "top") {
    .at <- NULL
    if (relToAt == "above" || relToAt == "left") .at <- at
    else if (relToAt == "centered") .at <- at - .range / 2
    else if (relToAt == "below" || relToAt == "right") .at <- at - .range
    
    qBoxHeight <- .range * 2 * qBoxMod
    
    yVals <- c(.at, .at + rescale(dens$y, c(0, max(dens$y)), c(0, .range)), .at)
    xVals <- c(minX, dens$x, maxX)
    
    qBoxPoints <- c(q1, .at, q3, .at + qBoxHeight)
    
    medX <- med
    medY <- .at + qBoxHeight / 2
    
  } else if (which == "bottom") {
    .at <- NULL
    if (relToAt == "above" || relToAt == "left") .at <- at + .range
    else if (relToAt == "centered") .at <- at + .range / 2
    else if (relToAt == "below" || relToAt == "right") .at <- at
    
    qBoxHeight <- .range * 2 * qBoxMod
    
    yVals <- c(.at, .at - rescale(dens$y, c(0, max(dens$y)), c(0, .range)), .at)
    xVals <- c(minX, dens$x, maxX)
    
    qBoxPoints <- c(q1, .at  - qBoxHeight, q3, .at)
    
    medX <- med
    medY <- .at - qBoxHeight / 2
    
  } else stop(paste("Did not recognize argument", paste("'", which, "'", sep = ""), "for parameter 'which'"))
  
  # if horizontal not true, turn the plot vertically
  if (!horizontal==TRUE) {
    tmp <- xVals
    xVals <- at + at - yVals
    yVals <- tmp
    
    qBoxPoints <- qBoxPoints[c(2,1,4,3)]
    qBoxPoints[1] <- at + at - qBoxPoints[1]
    qBoxPoints[3] <- at + at - qBoxPoints[3]
    
    tmp <- medX
    medX <- at + at - medY
    medY <- tmp
  }
  
  polygon(xVals, yVals, col = col, border = border)
  
  if (plotQBox == TRUE) rect(xleft = qBoxPoints[1], ybottom = qBoxPoints[2], xright = qBoxPoints[3], ytop = qBoxPoints[4],
                             col = colQBox, border = borderQBox)
  if (plotMed == TRUE) points(medX, medY, col=colMed, pch=pchMed)
}


#' @title histogram for marginalPlot
#' @description Internal function to plot histograms for the marginalPlot function
#' @param posterior vector of posterior values
#' @param prior vector of prior values or NULL
#' @param breaks integer, determining the number of breaks
#' @param at y position of the histograms
#' @param col vector determining posterior and prior color
#' @author Tankred Ott
histMarginal <- function (posterior, prior=NULL, breaks = 15, at=1, col=NULL) {
  
  maxHeight <- 0.90
  minHeight <- 0
  
  matPosterior <- createBreakMat(posterior, breaks)
  matPosterior[,3] <- BayesianTools:::rescale(matPosterior[,3], c(0, sum(matPosterior[,3])), c(minHeight, maxHeight))
  z <- maxHeight / max(matPosterior[,3])
  
  if (!is.null(prior)) {
    matPrior <- createBreakMat(prior, breaks)
    matPrior[,3] <- BayesianTools:::rescale(matPrior[,3], c(0, sum(matPrior[,3])), c(minHeight, maxHeight))
    
    z <- maxHeight / max(max(matPrior[,3]), max(matPosterior[,3]))
    matPrior[,3] <- matPrior[,3] * z
  }
  
  matPosterior[,3] <- matPosterior[,3] * z
  
  plotHist(matPosterior, at, col[1])
  if (!is.null(prior)) plotHist(matPrior, at, col[2])
}

#' @author Tankred Ott
plotHist <- function (x, at, col) {
  for (i in 1:nrow(x)) {
    xl <- x[i,1]
    xr <- x[i,2]
    yb <- 0
    yt <- x[i,3]
    rect(xleft = xl, xright = xr,
         ybottom = yb + at, ytop = yt + at,
         col = col)
  }
}

#' @author Tankred Ott
createBreakMat <- function (x, breaks=15) {
  cut_x <- cut(x, breaks = breaks)
  lvls <- levels(cut_x)
  ranges <- t(sapply(lvls, FUN = function (x) as.numeric(unlist(strsplit(substr(x, 2, nchar(x)-1), split = ",")))))
  frequencies <- as.vector(table(cut_x))
  mat <- cbind(ranges, frequencies)
  return(mat)
}


