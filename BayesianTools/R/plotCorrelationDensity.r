#' Flexible function to create correlation density plots
#' @author Florian Hartig
#' @param mat object of class "bayesianOutput" or a matrix or data frame of variables 
#' @param density type of plot to do. Either "smooth" (default), "corellipseCor", or "ellipse"
#' @param thin thinning of the matrix to make things faster. Default is to thin to 5000 
#' @param method method for calculating correlations. Possible choices are "pearson" (default), "kendall" and "spearman"
#' @param whichParameters indices of parameters that should be plotted
#' @param ... additional parameters to pass on to the \code{\link{getSample}}, for example parametersOnly =F, or start = 1000
#' @references The code for the correlation density plot originates from Hartig, F.; Dislich, C.; Wiegand, T. & Huth, A. (2014) Technical Note: Approximate Bayesian parameterization of a process-based tropical forest model. Biogeosciences, 11, 1261-1272.
#' @export
#' @seealso \code{\link{marginalPlot}} \cr
#'          \code{\link{plotTimeSeries}} \cr
#'          \code{\link{tracePlot}} \cr
#' @example /inst/examples/plotCorrelationDensity.R

correlationPlot<- function(mat, density = "smooth", thin = "auto", method = "pearson", whichParameters = NULL, ...){
  
  mat = getSample(mat, thin = thin, whichParameters = whichParameters, ...)
  
  numPars = ncol(mat)
  names = colnames(mat)
  
  panel.hist.dens <- function(x, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="blue4", ...)
  }
  
  # replaced by spearman 
  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = "complete.obs", method = method)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * abs(r))
  }
  
  plotEllipse <- function(x,y){ 
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    cor <- cor(x,y) 
    el = ellipse::ellipse(cor) 
    polygon(el[,1] + mean(x), el[,2] + mean(y), col = "red")
  }
  
  
  correlationEllipse <- function(x){
    cor = cor(x)
    ToRGB <- function(x){grDevices::rgb(x[1]/255, x[2]/255, x[3]/255)}
    C1 <- ToRGB(c(178, 24, 43))
    C2 <- ToRGB(c(214, 96, 77))
    C3 <- ToRGB(c(244, 165, 130))
    C4 <- ToRGB(c(253, 219, 199))
    C5 <- ToRGB(c(247, 247, 247))
    C6 <- ToRGB(c(209, 229, 240))
    C7 <- ToRGB(c(146, 197, 222))
    C8 <- ToRGB(c(67, 147, 195))
    C9 <- ToRGB(c(33, 102, 172))
    CustomPalette <- grDevices::colorRampPalette(rev(c(C1, C2, C3, C4, C5, C6, C7, C8, C9)))
    ord <- order(cor[1, ])
    xc <- cor[ord, ord]
    colors <- unlist(CustomPalette(100))
    ellipse::plotcorr(xc, col=colors[xc * 50 + 50])
  }
  
  if (density == "smooth"){ 
    return(pairs(mat, lower.panel=function(...) {par(new=TRUE);IDPmisc::ipanel.smooth(...)}, diag.panel=panel.hist.dens, upper.panel=panel.cor))
  }else if (density == "corellipseCor"){
      return(pairs(mat, lower.panel=plotEllipse, diag.panel=panel.hist.dens, upper.panel=panel.cor))  
  }else if (density == "ellipse"){
    correlationEllipse(mat)   
  }else if (density == F){
    return(pairs(mat, lower.panel=panel.cor, diag.panel=panel.hist.dens, upper.panel=panel.cor))      
  }else stop("wrong sensity argument")
  
}

