#' @author Maximilian Pichler
#' @title Diagnostic Plot
#' @describeIn This function plots the DIC, WAIC, mPSRF, PSRF(with upper C.I.) and traces of the parameters in dependence of iterations. DIC, WAIC are plotted separately for the chains and the trace plots also for the internal chains.
#' @param out object of class "bayesianOutput" 
#' @param lty line typ
#' @param lwd line width
#' @param numSamples for calculating WAIC, default = 10 because of high computational cost
#' @example /inst/examples/plotDiagnosticHelp.R
#' @export


plotDiagnostic <- function(out, lty = 1, lwd = 1, numSamples = 10){
  
  oldpar = NULL
  on.exit(par(oldpar))
  
  if(!"bayesianOutput" %in% class(out)) stop("Wrong input, object of class bayesianOutput required. see runMCMC()")
  calcWAIC <- TRUE
  
  if("mcmcSamplerList" %in% class(out) && out[[1]]$setup$pwLikelihood) calcWaic <- FALSE
  
  if("mcmcSampler" %in% class(out) && out$setup$pwLikelihood)     calcWaic <- FALSE
  
  # calculate DIC and WAIC
  if("mcmcSamplerList" %in% class(out)){
    
    if(is.matrix(out[[1]]$chain)) len <- out[[1]]$settings$iterations
    else len <- round(out[[1]]$settings$iterations / length(out[[1]]$chain))
    
    lenW <- length(seq(2 , by = 10, to = len))
    
    DICResult <- matrix(NA, nrow = length(out), ncol = len - 1)
    
    WAICResult<- matrix(NA, nrow = length(out), ncol = length(seq(2 , by = 10, to = len)))
    
    numPars   <- out[[1]]$setup$numPars
    
    for(i in 1:length(out)) {
      DICResult[i,] <- sapply(2 : len, FUN = function(x){return(DIC(out[[i]], end = x)$DIC)}) 
      if(calcWAIC) WAICResult[i,] <- sapply(seq(2 , by = 10, to = len), FUN = function(x){return(WAIC(out[[i]], end = x, numSamples = numSamples )$WAIC1)}) 
    }
    
  } else {
    if(is.matrix(out$chain)) len <- out$settings$iterations
    else len <- round(out$settings$iterations / length(out$chain))
    
    lenW<- length(seq(2 , by = 10, to = len))
    
    DICResult <- sapply(2 : len, FUN = function(x){return(DIC(out, end = x)$DIC)}) 
    
    if(calcWAIC) WAICResult<- sapply(seq(2 , by = 10, to = len), FUN = function(x){return(WAIC(out, end = x, numSamples = numSamples )$WAIC1)})
    
    numPars   <- out$setup$numPars
  }
  
  # calc mPSRF
  seq <- vector()
  for(i in 1:len){
    success <- try(coda::gelman.diag(getSample(out, parametersOnly = T, coda = T, end = i))$mpsrf, silent = T)
    if(!"try-error" %in% class(success)){
      # break
      seq[i] <- i
    }
  }
  seq <- seq[complete.cases(seq)]
  
  PSRF <- matrix(0, nrow = length(seq), ncol = numPars*2 + 1)

  for(i in 1:length(seq)){ 
    res <- coda::gelman.diag(getSample(out, parametersOnly = T, coda = T, end = seq[i]))
    PSRF[i,] <- c(res$psrf[1,], res$psrf[2,], res$mpsrf)
  }
  
  if(calcWAIC) par(mfrow = n2mfrow(length(out) + numPars + 1))
  else par(mfrow = n2mfrow(length(out) + numPars))
  
  # plot DIC
  plot(ylim = c(min(DICResult), max(DICResult)), x = 0, xlim = c(1,lenW), lty =  lty, lwd = lwd, type = "l", main = "DIC", xlab = "Iterations", ylab = "")
  if(is.matrix(DICResult)){
    for(i in 1:nrow(DICResult)){
      lines(y = DICResult[i,], x = 1:ncol(DICResult), lty = lty, lwd = lwd, col = rainbow(nrow(DICResult))[i])
    }
  } else {
    lines(y = DICResult, x = 1:length(DICResult), lty = lty, lwd = lwd, col = "red")
  }
  
  # plot WAIC
  if(calcWAIC){
    plot(ylim = c(min(WAICResult), max(WAICResult)), x = 0, xlim = c(1,lenW),lty =  lty, lwd = lwd, type = "l", main = "WAIC", xlab = "Iterations", ylab = "")
    if(is.matrix(WAICResult)){
      for(i in 1:nrow(WAICResult)){
        lines(y = WAICResult[i,], x = 1:ncol(WAICResult), lty = lty, lwd = lwd, col = rainbow(nrow(WAICResult))[i])
      }
    } else {
      lines(y = WAICResult, x = 1:length(WAICResult), lty = lty, lwd = lwd, col = "blue")
    }
  }
  
  # plot mPSRF
  if(!typeof(seq) == "logical") plot(y = PSRF[,ncol(PSRF)], x = seq, lty =  lty, lwd = lwd, type = "l", main = "mPSRF", xlab = "Iterations", ylab = "" )
  
  # plot PSRF
  
  plot(ylim = c(min(PSRF[, 1:(ncol(PSRF) - 1)]), max(PSRF[, c(seq(1, numPars*2,2))])), x = 0 , xlim = c(1, nrow(PSRF)), lty = lty, lwd = lwd, type = "l", main = "PSRF", xlab = "Iterations", ylab = "")
  for(i in 1:numPars){
    lines(y = PSRF[,i*2 - 1], x = seq, lty = lty, lwd = lwd, col = rainbow(numPars)[i])
    lines(y = PSRF[,i*2], x = seq, lty = 3, lwd = 0.7, col = rainbow(numPars)[i])
  }
  # plot parameter traces
  coda::traceplot(getSample(out, coda = T, parametersOnly = T))
}

