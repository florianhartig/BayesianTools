#' Trace plot for MCMC class
#' @param sampler an object of class MCMC sampler
#' @param thin determines the thinning interval of the chain
#' @param ... additional parameters to pass on to the \code{\link{getSample}}, for example parametersOnly = F, or start = 1000
#' @export
#' @seealso \code{\link{marginalPlot}} \cr
#'          \code{\link{plotTimeSeries}} \cr
#'          \code{\link{correlationPlot}}
#' @example /inst/examples/tracePlotHelp.R
tracePlot <- function(sampler, thin = "auto", ...){
  codaChain = getSample(sampler, coda = T, thin = thin, ...)
  plot(codaChain)
}
