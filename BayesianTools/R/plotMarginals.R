#' @export
marginalPlot <- function(x, ...) UseMethod("marginalPlot")

#' Plot MCMC marginals
#' @author Florian Hartig, Tankred Ott
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
#' @param singlePanel Logical, determining whether all histograms/violins should be plotted in a single plot panel or in separate panels.
#' @param dens Logical, determining wheter an density overlay should be plotted when 'histogram' is TRUE
#' @param col vector of colors for posterior and prior
#' @param lwd line width of the violin plots
#' @param ... additional parameters to pass on to the \code{\link{getSample}}
#' @export
#' @references 
#'          \code{\link{tracePlot}} \cr
#'          \code{\link{correlationPlot}}
#' @example /inst/examples/plotMarginals.R
marginalPlot <- function(mat, thin = "auto", scale = NULL, best = NULL, histogram = TRUE, plotPrior = TRUE, priorTop = FALSE,
                         nDrawsPrior = 1000, breaks=15, res=500, singlePanel=FALSE, dens=TRUE, col=c("#FF5000D0","#4682B4A0"),
                         lwd = par("lwd"), ...){
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
  
  
  if (singlePanel == TRUE) {
    plot(NULL, ylim = c(0,numPars +1), type="n", yaxt="n", xlab="", ylab="", xlim = xlim, main = main)
  } else {
    .main <- gsub("/n", "", main)
    panels <- getPanels(numPars)
    oldPar <- par(mfrow=panels, oma=c(4.1,0,2.7,0))
  }
  
  for (i in 1:numPars){
    if (singlePanel == TRUE) {
      main <- ""
      add <- TRUE
      .at <- i
    } else {
      main <- names[i]
      add <- FALSE
    }
    
    if (plotPrior) {
      if (histogram == TRUE) {
        # TODO: add overlay here
        histMarginal(x = list(mat[,i], priorMat[,i]), at = i, .range = 0.95, col = col, breaks = breaks, add = add,
                     main = main, dens = dens, res = res)
      } else {
        priorPos <- posteriorPos <- NULL
        if (priorTop == TRUE) {
          priorPos <- c("above", "top")
          posteriorPos <- c("below", "bottom")
        } else {
          priorPos <- c("below", "bottom")
          posteriorPos <- c("above", "top")
        }
        if (singlePanel == FALSE) {
          plot(NULL, xlim = range(mat[,i], priorMat[,i]), ylim = c(0,1), main=main, ylab = "frequencies", xlab = "values")
          .at <- 0.5
        }
        
        violinPlot(mat[,i], at = .at, .range = 0.475, add = T, col = col[1], relToAt = posteriorPos[1], which = posteriorPos[2], res=res, lwd = lwd)
        violinPlot(priorMat[,i], at = .at, .range = 0.475, add = T, col = col[2], relToAt = priorPos[1], which = priorPos[2], res=res, lwd = lwd)
      }
    } else {
      if (histogram == TRUE) {
        histMarginal(x = list(mat[,i]), at = i, .range = 0.95, col = col[1], breaks = breaks, add = add,
                     main = main, dens = dens, res = res)
      } else {
        violinPlot(mat[,i], at = i, .range = 0.95, add = add, col = col[1], relToAt = "centered", which = "both", res=res, main = main, lwd)
      }
    }
  
    if (singlePanel) axis(side = 2,at=i,labels = names[i],las=1)
  }
  if(is.numeric(best)) points(best,1:length(best), cex = 3, pch = 4, lwd = 2)
  
  if(singlePanel == TRUE) {
    if(plotPrior && histogram) legend("bottomright", c("posterior", "prior"), col = col, pch = 15, cex = 0.8, bty = "n")
    if(plotPrior && !histogram) legend("topright", c("posterior", "prior"), col = col, pch = 15, cex = 0.8, bty = "n")
  } else {
    mtext(text = .main, side = 3, line = -1.25, outer = T, font = 2)
    
    # overlay plot with empty plot to be able to place the legends freely
    oldPar2 <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    
    if(plotPrior && histogram) legend("bottom", c("posterior", "prior"), xpd = TRUE, horiz = TRUE, inset = c(0, 0),
                                      bty = "n", pch = 15, col = col, cex = 2)
    if(plotPrior && !histogram) legend("bottom", c("posterior", "prior"), xpd = TRUE, horiz = TRUE, inset = c(0, 0),
                                       bty = "n", pch = 15, col = col, cex = 2)
    par(oldPar2)
    par(oldPar)
  }
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
#' @param lwd line width if the border of the violin
#' @param colQBox color of quantile box
#' @param borderQBox color of the border of the quantile box
#' @param colMed color of the median point
#' @param pchMed pch for median point
#' @param res "resolution" of the violin. Determining how many descrete points should be used to calculate the density kernel.
#' @param main header text, only applicable, if 'add' is FALSE
violinPlot <- function (x, at, .range = 1, add = FALSE, horizontal = TRUE, which = "both", relToAt = "above", plotQBox = TRUE, plotMed = TRUE,
                        col = "orangered", border = "black", lwd = par("lwd"), colQBox = "black", borderQBox = "black", colMed = "white",
                        pchMed = 19, res = 500, main = "") {
  
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
    if (horizontal == TRUE) plot(NULL, xlim = c(minX, maxX), ylim = c(0,1), xlab = "", ylab = "", main = main)
    else plot(NULL, xlim = c(0,1), ylim = c(minX, maxX), xlab = "", ylab = "", main = main)
    
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
  
  # draw the plot
  polygon(xVals, yVals, col = col, border = border, lwd = lwd)
  
  if (plotQBox == TRUE) rect(xleft = qBoxPoints[1], ybottom = qBoxPoints[2], xright = qBoxPoints[3], ytop = qBoxPoints[4],
                             col = colQBox, border = borderQBox)
  if (plotMed == TRUE) points(medX, medY, col = colMed, pch = pchMed)
}


#' @title histogram for marginalPlot
#' @description Internal function to plot histograms for the marginalPlot function
#' @param x list of vectors. Each vector will be plotted as a separate histogram.
#' @param at y position of the histograms
#' @param .range maximum height of the histogram. If NULL, will be determined automatically
#' @param breaks integer, determining the number of breaks
#' @param col vector of colors determining histogram colors
#' @param border vector of colors determining histogram border colors
#' @param add Logical, determining whether the histogram should be added to and existing plot window.
#' @param main Character, determining the title of the plot
#' @param dens Logical, determining whether a density plot should be plotted additionally to the histogram
#' @param densCol vector of colors for the density plots
#' @param densLwd line width of the density plot
#' @param densLty line type of the density plot
#' @param res resolution of the density overlay
#' @author Tankred Ott
# TODO: let the function automatically chose colors
histMarginal <- function (x, at = 0, .range = NULL, breaks = 15, col=c("#FF5000C0", "#4682B4A0"),
                          border = c("black", "black"),  main="", dens=FALSE, densCol = c("black", "black"),
                          densLwd = c(2, 2), densLty = c(2, 2), add=TRUE, res=500) {
  
  mnX <- Inf
  mxX <- -Inf
  
  mnY <- 0
  mxY <- -Inf
  
  densities <- NULL
  
  breakMats <- rep(list(NA), length(x))
  if (dens) densities <- rep(list(NA), length(x))
  
  for (i in 1:length(x)) {
    mnX <- min(mnX, min(x[[i]]))
    mxX <- max(mxX, max(x[[i]]))
    
    if (dens == TRUE) {
      densities[[i]] <- density(x[[i]], n = res, from = mnX, to = mxX)
      mxY <- max(mxY, max(densities[[i]]$y))
    }
  }
  
  brkRange <- round((mxX - mnX) / breaks, digits = 5)
  brks <- mnX + cumsum(c(0, rep(brkRange, breaks)))
  
  for (i in 1:length(x)) {
    breakMats[[i]] <- createBreakMat(x[[i]], brks, TRUE)
    mxY <- max(mxY, max(breakMats[[i]][,3]))
  }
  
  if (!is.null(.range)) {
    for (i in 1:length(breakMats)) {
      breakMats[[i]][,3] <- rescale(breakMats[[i]][,3], from = c(0, mxY), to = c(0, .range))
      if (dens == TRUE) densities[[i]]$y <- rescale(densities[[i]]$y, from = c(0, mxY), to = c(0, .range))
    }
    mxY <- .range
  }
  
  if (add == FALSE) {
    at <- 0
    plot(NULL, xlim = c(mnX, mxX), ylim = c(mnY, mxY), main = main, xlab = "value", ylab = "frequency")
  }
  
  for (i in 1:length(breakMats)) {
    plotHist(breakMats[[i]], at = at, col = col[i], border = border[i])
    if (dens == TRUE) lines(densities[[i]]$x, densities[[i]]$y + at, col = densCol[i], lty = densLty[i], lwd = densLwd[i])
  }
  
}

#' @author Tankred Ott
#' @title plot histogram
#' @description A simple custom histogram plotting function 
#' @param x matrix (e.g. constructed with \code{\link{createBreakMat}}) or dataframe containing the lower bounds, the upper bounds, and the frequencies of the breaks in the columns and the individual breaks in the rows. 
#' @param at y position of the histogram. If NULL a new plot will be generated.
#' @param col color of the histogram
#' @param border border color
plotHist <- function (x, at=NULL, col="orangered", border="black") {
  if (is.null(at)) {
    at <- 0
    xlim <- c(min(x[,1:2]), max(x[,1:2]))
    ylim <- c(0, max(x[,3]))
    plot(x = NULL, xlim = xlim, ylim = ylim, xlab = "values", ylab = "frequencies")
  }
  
  for (i in 1:nrow(x)) {
    xl <- x[i,1]
    xr <- x[i,2]
    yb <- 0
    yt <- x[i,3]
    rect(xleft = xl, xright = xr, ybottom = yb + at, ytop = yt + at,
         col = col, border = border)
  }
}

#' @author Tankred Ott
#' @title create break matrix
#' @description A function for use with plotHist. Creates a matrix representing breaks of a histogram. The matrix will contain the upper bounds, the lower bounds and the frequencies of the breaks in the columns, and the individual breaks in the rows.
#' @param x vector of values
#' @param breaks number of breaks
#' @param scale logical, if TRUE, the area within the rectangle will be scaled to one (density)
createBreakMat <- function (x, breaks=15, scale=FALSE) {
  cut_x <- cut(x, breaks = breaks)
  lvls <- levels(cut_x)
  ranges <- t(sapply(lvls, FUN = function (x) as.numeric(unlist(strsplit(substr(x, 2, nchar(x)-1), split = ",")))))
  frequencies <- as.vector(table(cut_x))
  mat <- cbind(ranges, frequencies)
  if (scale == TRUE) {
    areas <- apply(mat, 1, function(x) return ((x[2] - x[1]) * x[3]))
    mat[,3] <- rescale(mat[,3], from=c(0, sum(areas)), to=c(0,1))
  }
  colnames(mat) <- c("lower", "upper", "frequencies")
  return(mat)
}


