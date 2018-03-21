#' @author Maximilian Pichler
#' @title Diagnostic Plot
#' @describeIn This function plots the DIC, WAIC, mPSRF, PSRF(with upper C.I.) and traces of the parameters in dependence of iterations. DIC, WAIC are plotted separately for the chains and the trace plots also for the internal chains.
#' @param out object of class "bayesianOutput"
#' @param start start value for showing DIC, WAIC, mPSRF and PSRF, default = 50
#' @param plotWAIC whether to calculate WAIC or not, default = T
#' @param plotPSRF calculate and plot mPSRF/PSRF or not, default = T 
#' @param plotDIC calculate and plot DICor not, default = T 
#' @param plotTrace show trace plots or not, default = T
#' @param lty line typ
#' @param lwd line width
#' @param numSamples for calculating WAIC, default = 10 because of high computational cost
#' @example /inst/examples/plotDiagnosticHelp.R
#' @export


plotDiagnostic <- function(out, lty = 1, lwd = 1, start = 50, numSamples = 10, plotWAIC = T, plotPSRF = T, plotDIC = T, plotTrace = T){
  
  oldpar = NULL
  on.exit(par(oldpar))
  
  
  
  if(!"bayesianOutput" %in% class(out)) stop("Wrong input, object of class bayesianOutput required. see runMCMC()")
  
  calcWAIC <- TRUE
  
  if("mcmcSamplerList" %in% class(out) && out[[1]]$setup$pwLikelihood) calcWaic <- FALSE
  
  if("mcmcSampler" %in% class(out) && out$setup$pwLikelihood)     calcWaic <- FALSE
  
  if(!plotWAIC) calcWaic <- FALSE
  
  
  
  
  
  # calculate DIC and WAIC, minimum range: start - start+1
  if("mcmcSamplerList" %in% class(out)){

    if(is.matrix(out[[1]]$chain)) len <- out[[1]]$settings$iterations
    else len <- round(out[[1]]$settings$iterations / length(out[[1]]$chain))

    lenW <- length(seq(2 , by = 10, to = len))

    DICResult <- matrix(NA, nrow = length(out), ncol = len - 2)

    WAICResult<- matrix(NA, nrow = length(out), ncol = length(seq(2 , by = 10, to = len)))

    numPars   <- out[[1]]$setup$numPars
    
    Wseq <- seq(2 , by = 10, to = len)

    for(i in 1:length(out)) {
      if(plotDIC) DICResult[i,] <- sapply(2:len, FUN = function(x){return(DIC(out[[i]], end = x)$DIC)})
      if(calcWAIC) WAICResult[i,] <- sapply(seq(2 , by = 10, to = len), FUN = function(x){return(WAIC(out[[i]], end = x, numSamples = numSamples )$WAIC1)})
    }
    
  } else {
    if(is.matrix(out$chain)) len <- out$settings$iterations
    
    else len <- round(out$settings$iterations / length(out$chain))

    lenW<- length(seq(2, by = 10, to = len))
    
    Wseq <- seq(2 , by = 10, to = len)

    if(plotDIC) DICResult <- sapply(2:len, FUN = function(x){return(DIC(out,  end = x)$DIC)})

    if(calcWAIC) WAICResult<- sapply(seq(2, by = 10, to = len), FUN = function(x){return(WAIC(out, end = x, numSamples = numSamples )$WAIC1)})

    numPars   <- out$setup$numPars
  }
  
  
  
  
  # calc mPSRF, first checking which low values we could calculate
  if(plotPSRF){
    seq <- vector()
    for(i in 2:len){
      success <- try(coda::gelman.diag(getSample(out, parametersOnly = T, coda = T,  end = i))$mpsrf, silent = T)
      if(!"try-error" %in% class(success)){
        # break
        seq[i] <- i
      }
    }
    seq <- seq[complete.cases(seq)]
    
    # calculate the actual PSRF values
    PSRF <- matrix(0, nrow = length(seq), ncol = numPars*2 + 1)
  
    for(i in 1:length(seq)){ 
      res <- coda::gelman.diag(getSample(out, parametersOnly = T, coda = T,  end = seq[i]))
      PSRF[i,] <- c(res$psrf[1,], res$psrf[2,], res$mpsrf)
    }
  }
  
  
  # Get number of plots
  nrPlots <- 1
  if(calcWAIC) nrPlots <- nrPlots + 1
  if(plotDIC)  nrPlots <- nrPlots + 1
  if(plotPSRF) nrPlots <- nrPlots + 2
  if(plotTrace) nrPlots<- numPars*2 + nrPlots
  par(mfrow = getPanels(nrPlots))
  
  
  
  
  # plot DIC
  if(plotDIC){
    
    if(is.matrix(DICResult)){
      col <- 1:ncol(DICResult)
      matplot(y = DICResult[start:nrow(DICResult),], x = start:nrow(?plot?ploDICResult),   lty =  lty, lwd = lwd, type = "l", main = "DIC", xlab = "Iterations", ylab = "", col = col)
      
    }else{
      col <- "red"
      
      plot(y = DICResult[start:length(DICResult)], x = start:length(DICResult),  lty =  lty, lwd = lwd, type = "l", main = "DIC", xlab = "Iterations", ylab = "", col = "red")
    }
  }
  
  
  # plot WAIC
  if(calcWAIC){
    
    plotSeq <- which(Wseq > start, arr.ind = T)
    
    if(is.matrix(WAICResult)){
      col <- 1:ncol(WAICResult)
      matplot(y = WAICResult[plotSeq,],x = plotSeq,lty =  lty, lwd = lwd, type = "l", main = "WAIC", xlab = "Iterations", ylab = "", col = col, xaxt = "n")
      axis(1, at = plotSeq[seq(1, by = 10, to = length(plotSeq))], labels = plotSeq[seq(1, by = 10, to = length(plotSeq))]*10)
    } else {
      plot(y = WAICResult[plotSeq],x = plotSeq,lty =  lty, lwd = lwd, type = "l", main = "WAIC", xlab = "Iterations", ylab = "", col = "red", xaxt = "n")
      axis(1, at = plotSeq[seq(1, by = 10, to = length(plotSeq))], labels = plotSeq[seq(1, by = 10, to = length(plotSeq))]*10)
    }
    
    
    
  }
  
  # plot mPSRF
  if(plotPSRF){
    if(!typeof(seq) == "logical" ) plot(y = PSRF[start:nrow(PSRF),ncol(PSRF)], x = start:nrow(PSRF),lty =  lty, lwd = lwd, type = "l", main = "mPSRF", xlab = "Iterations", ylab = "" )
    
    # plot PSRF
    plotSeqP <- which(seq > start, arr.ind = T)
    plot(ylim = c(min(PSRF[plotSeqP, 1:(ncol(PSRF) - 1)]), max(PSRF[plotSeqP, c(seq(1, numPars*2,2))])*1.1), x = 0, xlim = c(min(plotSeqP), max(plotSeqP)),  lty = lty, lwd = lwd, type = "l", main = "PSRF", xlab = "Iterations", ylab = "")
    for(i in 1:numPars){
      lines(y = PSRF[plotSeqP, i*2 - 1], x = plotSeqP, lty = lty, lwd = lwd, col = i) # PSRF
      lines(y = PSRF[plotSeqP, i*2],  x = plotSeqP, lty = 3, lwd = 0.7, col = i) # C.I.
    }
  }
    # plot parameter traces
  if(plotTrace) coda::cumuplot(getSample(out, coda = T, parametersOnly = T), ask = F, auto.layout = F)
}

