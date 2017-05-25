#' Plots a time series, with the option to include confidence and prediction band
#' @author Florian Hartig
#' @param x optional values for x axis (time)
#' @param observed observed values
#' @param predicted predicted values
#' @param confidenceBand matrix with confidenceBand
#' @param predictionBand matrix with predictionBand
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param ... further arguments passed to \code{\link[graphics]{plot}}
#' @seealso \code{\link{plotTimeSeriesResults}}  \cr
#'          \code{\link{marginalPlot}} \cr
#'          \code{\link{tracePlot}} \cr
#'          \code{\link{correlationPlot}}
#' @example /inst/examples/plotTimeSeriesHelp.R
#' @export
plotTimeSeries <- function(observed = NULL, predicted = NULL, x = NULL, confidenceBand = NULL, predictionBand = NULL, xlab = "Time", ylab = "Observed / predicted values", ...){
  
  ylim = range(observed, predicted, confidenceBand, predictionBand,na.rm=T)
  
  if (is.null(x)){
    if(!is.null(observed)) x = 1:length(observed)
    else if(!is.null(predicted)) x = 1:length(predicted)
    else stop("either observed or predicted must be supplied")
  }
  
  len = length(x)
  
  plot(x, ylim = ylim, type = "n", xlab = xlab, ylab = ylab, ...)
  
  if(!is.null(predictionBand)) polygon(c(1:len,len:1),c(predictionBand[1,],predictionBand[2,len:1]),col="moccasin",border=NA)
  
  if(!is.null(confidenceBand)) polygon(c(1:len,len:1),c(confidenceBand[1,],confidenceBand[2,len:1]),col="#99333380",border=NA)    
    
  if(!is.null(predicted)) lines(predicted, col = "red")
  if(!is.null(observed)) points(observed, col = "black", pch = 3, cex = 0.6)
  
}



#' Plots residuals of a time series
#' @author Florian Hartig
#' @param residuals x
#' @param x optional values for x axis (time)
#' @param main title of the plot
#' @export
plotTimeSeriesResiduals <- function(residuals, x = NULL, main = "residuals"){
  
  ylim = range(residuals)
  
  if (is.null(x)){
    x  = 1:length(residuals)
  }
  barplot(residuals)
}


#' Creates a time series plot typical for an MCMC / SMC fit
#' @author Florian Hartig
#' @param sampler Either a) a matrix b) an MCMC object (list or not), or c) an SMC object
#' @param model function that calculates model predictions for a given parameter vector
#' @param observed observed values
#' @param error function with signature f(mean, par) that generates error expectations from mean model predictions. Par is a vector from the matrix with the parameter samples (full length). f needs to know which of these parameters are parameters of the error function
#' @param start numeric start value for the plot (see \code{\link{getSample}})
#' @param plotResiduals logical determining whether residuals should be plotted
#' @param prior if a prior sampler is implemented, setting this parameter to T will draw model parameters from the prior instead of the posterior distribution
#' @export
plotTimeSeriesResults <- function(sampler, model, observed, error = NULL, plotResiduals = T, start = 1, prior = F){
  
  if(prior == F){
    if(inherits(sampler,"bayesianOutput")) parMatrix = getSample(sampler, start = start)
    else if (class(sampler) == "matrix") parMatrix = sampler
    else stop("wrong type give to variable sampler")    
  }else if (prior == T){
    if(class(sampler)[1] == "mcmcSamplerList") parMatrix = sampler[[1]]$setup$prior$sampler(1000)
    else parMatrix = sampler$setup$prior$sampler(1000)
  }else stop("BayesianTools::plotTimeSeriesResults - wrong argument to prior")

  thin = min(1000, nrow(parMatrix))
  
  pred <- getPredictiveIntervals(parMatrix = parMatrix, model = model, thin = thin, quantiles = c(0.025, 0.5, 0.975), error = error)
  
  if(!is.null(error)) plotTimeSeries(observed = observed, predicted = pred[2,], confidenceBand = pred[c(1,3),], predictionBand = pred[c(4,6),] )
  else plotTimeSeries(observed = observed, predicted = pred[2,], confidenceBand = pred[c(1,3),])
  
  
  # TODO - plotResiduals needs to be implemented
  
}









