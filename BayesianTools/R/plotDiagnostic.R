#' @author Maximilian Pichler
#' @title Diagnostic Plot
#' @description  This function plots the DIC, WAIC, mPSRF, PSRF(with upper C.I.) and traces of the parameters in dependence of iterations. DIC, WAIC are plotted separately for the chains and the trace plots also for the internal chains.
#' @param out object of class "bayesianOutput"
#' @param start start value for calculating DIC, WAIC, mPSRF and PSRF, default = 50
#' @param numSamples for calculating WAIC, default = 10 because of high computational costs
#' @param window plot range to show, vector of percents or only one value as start value for the window
#' @param plotWAIC whether to calculate WAIC or not, default = T
#' @param plotPSRF calculate and plot mPSRF/PSRF or not, default = T 
#' @param plotDIC calculate and plot DICor not, default = T 
#' @param plotTrace show trace plots or not, default = T
#' @param graphicParameters graphic parameters as list for plot function
#' @param ... parameters to give to getSample
#' @example /inst/examples/plotDiagnosticHelp.R
#' @export



plotDiagnostic <- function(out, start = 50, numSamples = 100, window = 0.2, plotWAIC = F, plotPSRF = T, plotDIC = T, plotTrace = T, graphicParameters = NULL, ...){
  
  oldpar = NULL
  on.exit(par(oldpar))
  
  
  
  if(!"bayesianOutput" %in% class(out)) stop("Wrong input, object of class bayesianOutput required. see runMCMC()")
  
  calcWAIC <- TRUE
  
  if("mcmcSamplerList" %in% class(out) && out[[1]]$setup$pwLikelihood) calcWAIC <- FALSE
  
  if("mcmcSampler" %in% class(out) && out$setup$pwLikelihood)     calcWAIC <- FALSE
  
  if(!plotWAIC) calcWAIC <- FALSE
  
  
  defaultGraphicParameters <- graphicParameters
  
  
  # calculate DIC and WAIC, minimum range: start - start+1
  if("mcmcSamplerList" %in% class(out)){

    if(is.matrix(out[[1]]$chain)) len <- out[[1]]$settings$iterations
    else len <- round(out[[1]]$settings$iterations / length(out[[1]]$chain))
    
    iter = out[[1]]$settings$iterations
    
    internal = length(out[[1]]$chain)
    
    start = start + 1
    
    lenW <- length(seq(start , by = 10, to = len))

    DICResult <- matrix(NA, nrow = length(out), ncol = len - start)

    WAICResult<- matrix(NA, nrow = length(out), ncol = length(seq(start , by = 10, to = len)))

    numPars   <- out[[1]]$setup$numPars
    
    Wseq <- seq(start , by = 10, to = len)

    for(i in 1:length(out)) {
      if(plotDIC) DICResult[i,] <- sapply(start:len, FUN = function(x){return(DIC(out[[i]], start = start - 1 , end = x, ...)$DIC)})
      if(calcWAIC) WAICResult[i,] <- sapply(seq(start , by = 10, to = len), FUN = function(x){return(WAIC(out[[i]], start = start - 1 ,end = x, numSamples = numSamples, ...)$WAIC1)})
    }
    
  } else {
    if(is.matrix(out$chain)) len <- out$settings$iterations
    
    else len <- round(out$settings$iterations / length(out$chain))
    
    internal = length(out$chain)
    
    iter = out$settings$iterations
    
    start = start + 1
    
    lenW<- length(seq(start, by = 10, to = len))
    
    Wseq <- seq(start , by = 10, to = len)

    if(plotDIC) DICResult <- sapply(start:len, FUN = function(x){return(DIC(out, start = start - 1, end = x, ...)$DIC)})

    if(calcWAIC) WAICResult<- sapply(seq(start, by = 10, to = len), FUN = function(x){return(WAIC(out, end = x, start = start - 1, numSamples = numSamples, ...)$WAIC1)})

    numPars   <- out$setup$numPars
  }
  
  
  # RB: missing: check if sampler with multiple chains
  # should user call method with plotPSFR=F for one-chain-sampler?
  
  # calc mPSRF, first checking which low values we could calculate
  if(plotPSRF){
    
    seq <- vector()
    for(i in start:len){
      success <- try(coda::gelman.diag(getSample(out, start = start - 1, parametersOnly = T, coda = T,  end = i, ...))$mpsrf, silent = T)
      if(!"try-error" %in% class(success)){
        # break
        seq[i] <- i
      }
    }
    seq <- seq[complete.cases(seq)]
    
    # calculate the actual PSRF values
    if(numPars > 1) PSRF <- matrix(0, nrow = length(seq), ncol = numPars*2 + 1)
    else PSRF <- matrix(0, nrow = length(seq), ncol = numPars*2 )
  
    for(i in 1:length(seq)){ 
      res <- coda::gelman.diag(getSample(out, start = start - 1, parametersOnly = T, coda = T, end = seq[i], ...))
      if(numPars > 1)PSRF[i,] <- c(as.vector(res$psrf), res$mpsrf)
      else PSRF[i,] <- c(as.vector(res$psrf))
    }
  }
  
  
  # Get number of plots
  nrPlots <- 2
  if(calcWAIC) nrPlots <- nrPlots + 1
  if(plotDIC)  nrPlots <- nrPlots + 1
  if(plotPSRF) nrPlots <- nrPlots + 2
  if(plotTrace) nrPlots<- numPars*2 + nrPlots
  par(mfrow = getPanels(nrPlots))
  
  
  
  
  # set graphicParameters
  if(is.null(graphicParameters)){
    graphicParameters = list(lty = 1, lwd = 1, type = "l", xlab = "Iterations", ylab = "", col = 1:6)
  } else {
    if(is.null(graphicParameters$lty)) graphicParameters$lty = 1
    if(is.null(graphicParameters$lwd)) graphicParameters$lwd = 1
    if(is.null(graphicParameters$type)) graphicParameters$type = "l"
    if(is.null(graphicParameters$xlab)) graphicParameters$xlab = "Iterations"
    if(is.null(graphicParameters$ylab)) graphicParameters$ylab = ""
    if(is.null(graphicParameters$col)) graphicParameters$col = 1:6
  }
  

  
  # plot DIC
  if(plotDIC){
    
    
    if(is.matrix(DICResult)){
      # col <- 1:ncol(DICResult)
      if(is.na(window[2])) endV <- nrow(DICResult)
      else endV <- window[2]*nrow(DICResult)
      startV <- window[1]*nrow(DICResult)
      x = nrow(DICResult)
      ylim = c(min(DICResult[startV:endV,])*0.99, max(DICResult[startV:endV,])*1.01)
    } else {
      if(is.na(window[2])) endV <- length(DICResult)
      else endV <- window[2]*length(DICResult)
      startV <- window[1]*length(DICResult)
      x = length(DICResult)
      ylim = c(min(DICResult[startV:endV])*0.99, max(DICResult[startV:endV])*1.01)
    }
      graphicParameters$y = DICResult
      graphicParameters$x = 1:x
      graphicParameters$main = "DIC"
      graphicParameters$xlim = c(startV, endV)
      graphicParameters$ylim = ylim
      if(is.null(graphicParameters$xaxt)) graphicParameters$xaxt = "n" 
      do.call(matplot, graphicParameters)
      if(graphicParameters$xaxt == "n" ){
        axis(1, at = seq(startV, by = 100, to = endV), labels = seq(startV, by = 100, to = endV)*internal)
        graphicParameters$xaxt <- NULL
      }
  }
  
  
  # plot WAIC
  if(calcWAIC){
    if(is.matrix(WAICResult)){
      # col <- 1:ncol(DICResult)
      if(is.na(window[2])) endV <- nrow(WAICResult)
      else endV <- window[2]*nrow(WAICResult)
      startV <- window[1]*nrow(WAICResult)
      x = nrow(WAICResult)
      ylim = c(min(WAICResult[startV:endV,])*0.99, max(WAICResult[startV:endV,])*1.01)
    } else {
      if(is.na(window[2])) endV <- length(WAICResult)
      else endV <- window[2]*length(WAICResult)
      startV <- window[1]*length(WAICResult)
      x = length(WAICResult)
      ylim = c(min(WAICResult[startV:endV])*0.99, max(WAICResult[startV:endV])*1.01)
    }
    graphicParameters$y = WAICResult
    graphicParameters$x = 1:x
    graphicParameters$main = "WAIC"
    graphicParameters$xlim = c(startV, endV)
    graphicParameters$ylim = ylim
    if(is.null(graphicParameters$xaxt)) graphicParameters$xaxt = "n" 
    do.call(matplot, graphicParameters)
    if(graphicParameters$xaxt == "n" ){
      axis(1, at = seq(startV, by = 10, to = endV), labels = seq(startV, by = 10, to = endV)*10*internal)
      graphicParameters$xaxt <- NULL
    }
    
  }
  

  if(plotPSRF){
    if(is.na(window[2])) endV <- nrow(PSRF)
    else endV <- window[2]*nrow(PSRF)
    startV <- window[1]*nrow(PSRF)
    graphicParameters$xlim = c(startV, endV)
    graphicParameters$x = 1:nrow(PSRF)
    # plot mPSRF
    if(numPars > 1){
      if(!typeof(seq) == "logical" ) {
        
        graphicParameters$ylim = c(min(PSRF[startV:endV,ncol(PSRF)])*0.99, max(PSRF[startV:endV,ncol(PSRF)])*1.01)
        graphicParameters$y = PSRF[,ncol(PSRF)]
        graphicParameters$main = "mPSRF"
        do.call(plot, graphicParameters)
      }
    }
    
    graphicParameters$ylim = c(min(PSRF[startV:endV,-ncol(PSRF)])*0.99, max(PSRF[startV:endV,-ncol(PSRF)])*1.01)
    graphicParameters$y = PSRF[,-ncol(PSRF)]
    graphicParameters$main = "PSRF"
    
    lty = NULL
    for(i in 1:numPars)lty <- c(lty, c(1,2))
    graphicParameters$lty <- lty
    
    col = NULL
    for(i in 1:6)col <- c(col, c(i,i))
    graphicParameters$col <- col
    
    do.call(matplot, graphicParameters)
   
  }
    # plot parameter traces
  if(plotTrace){
    # if(is.null(defaultGraphicParameters)) defaultGraphicParameters <- list()
    # if(is.na(window[2])) endV <- len
    # else endV <- window[2]*len
    # defaultGraphicParameters$xlim <- c(len*window[1], endV)
    # defaultGraphicParameters$ask = F
    # defaultGraphicParameters$auto.layout = F
    # defaultGraphicParameters$x = getSample(out, start = start, coda = T, parametersOnly = T,...)
    # do.call(coda::cumuplot, defaultGraphicParameters)
    
    coda::cumuplot(getSample(out, start = start, coda = T, parametersOnly = T, ...), ask = F, auto.layout = F)
  }
}

