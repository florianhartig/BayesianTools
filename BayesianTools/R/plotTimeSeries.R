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
  
  ylim = range(observed, predicted, confidenceBand, predictionBand,na.rm=TRUE)
  
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
#' @param error function with signature f(mean, par) that generates observations with error (error = stochasticity according to what is assumed in the likelihood) from mean model predictions. Par is a vector from the matrix with the parameter samples (full length). f needs to know which of these parameters are parameters of the error function. See example in \code{\link{VSEM}}
#' @param start numeric start value for the plot (see \code{\link{getSample}})
#' @param plotResiduals logical determining whether residuals should be plotted
#' @param prior if a prior sampler is implemented, setting this parameter to TRUE will draw model parameters from the prior instead of the posterior distribution
#' @param ... further arguments passed to \code{\link[graphics]{plot}}
#' @export
plotTimeSeriesResults <- function(sampler, model, observed, error = NULL, plotResiduals = TRUE, start = 1, prior = FALSE, ...){
  oldPar = par(no.readonly = TRUE)
  
  if (plotResiduals == TRUE && is.null(error)) {
    warning("Can not plot residuals without an error function.")
  }
  
  if (plotResiduals == TRUE && !is.null(error)) {
    layout(matrix(c(1, 1, 1, 2, 3, 4), 2, 3, byrow = TRUE))
    par(mar = c(3, 3, 3, 3), oma = c(2, 2, 2, 2))
  }
  
  # ... can we pass on to both getSample and plot? 
  
  if(prior == FALSE){
    if(inherits(sampler,"bayesianOutput")) parMatrix = getSample(sampler, start = start)
    else if (class(sampler) == "matrix") parMatrix = sampler
    else if ("mcmc.list" %in% class(sampler) || "mcmc" %in% class(sampler)) parMatrix = getSample(sampler, start = start)
    else stop("wrong type given to variable sampler")    
  }else if (prior == TRUE){
    if(inherits(sampler,"bayesianOutput")) {
      if(class(sampler)[1] == "mcmcSamplerList") parMatrix = sampler[[1]]$setup$prior$sampler(1000)
      else parMatrix <- sampler$setup$prior$sampler(1000)
    } else {
      stop("prior==TRUE is only available for sampler of type bayesianOutput") 
    }
  }else stop("BayesianTools::plotTimeSeriesResults - wrong argument to prior")

  numSamples = min(1000, nrow(parMatrix))
  
  pred <- getPredictiveIntervals(parMatrix = parMatrix,
                                 model = model,
                                 numSamples = numSamples,
                                 quantiles = c(0.025, 0.5, 0.975),
                                 error = error)
  
  if(!is.null(error)) plotTimeSeries(observed = observed,
                                     predicted = pred$posteriorPredictivePredictionInterval[2,],
                                     confidenceBand = pred$posteriorPredictiveCredibleInterval[c(1,3),],
                                     predictionBand = pred$posteriorPredictivePredictionInterval[c(1,3),],
                                     ...)
  else plotTimeSeries(observed = observed, predicted = pred$posteriorPredictiveSimulations,
                      confidenceBand = pred$posteriorPredictiveSimulations[c(1,3),],
                      ...)
  
  if (plotResiduals && !is.null(error)) {
    dh = getDharmaResiduals(model = model,
                            parMatrix = parMatrix,
                            numSamples = numSamples,
                            observed = observed,
                            error = error,
                            plot = FALSE)
    # qq-plot
    gap::qqunif(dh$scaledResiduals, pch=2, bty="n", logscale = F, col = "black", cex = 0.6, main = "QQ plot residuals", cex.main = 1)
    
    # residuals vs fitted
    DHARMa::plotResiduals(dh$fittedPredictedResponse, dh$scaledResiduals, main = "Residual vs. predicted\n quantile lines should be\n horizontal lines at 0.25, 0.5, 0.75", cex.main = 1, xlab = "Predicted value", ylab = "Standardized residual")
    
    # residuals vs time
    t <- 1:length(dh$fittedPredictedResponse)
    DHARMa::plotResiduals(t, dh$scaledResiduals, xlab = "Time", ylab = "Standardized residual", main = "Residual vs. time\n quantile lines should be\n horizontal lines at 0.25, 0.5, 0.75", cex.main = 1)
    
    message("DHARMa::plotTimeSeriesResults called with posterior predictive (residual) diagnostics. Type vignette(\"DHARMa\", package=\"DHARMa\") for a guide on how to interpret these plots")

  }
  
  
  par(oldPar)
}

#' Creates a DHARMa object
#' @author Tankred Ott
#' @param model function that calculates model predictions for a given parameter vector
#' @param parMatrix a parameter matrix from which the simulations will be generated
#' @param numSamples the number of samples
#' @param observed a vector of observed values
#' @param error function with signature f(mean, par) that generates error expectations from mean model predictions. Par is a vector from the matrix with the parameter samples (full length). f needs to know which of these parameters are parameters of the error function
#' @param plot logical, determining whether the simulated residuals should be plotted
# #' @export
getDharmaResiduals <- function(model, parMatrix, numSamples, observed, error, plot = TRUE){

  predDistr <- getPredictiveDistribution(parMatrix = parMatrix,
                                         model = model,
                                         numSamples = numSamples)
  # apply error to predictions
  for (i in 1:nrow(predDistr)){
    predDistr[i,] = error(mean = predDistr[i,], par = parMatrix[i,]) 
  }
  
  fittedPars = apply(parMatrix, 2, median)
  fittedPredictedResponse = model(fittedPars)
  
  dh = DHARMa::createDHARMa(simulatedResponse = t(predDistr),
                            observedResponse = observed,
                            fittedPredictedResponse = fittedPredictedResponse)
  if (plot == TRUE) {
    DHARMa::plotSimulatedResiduals(dh)
  }

  return(dh)
}













