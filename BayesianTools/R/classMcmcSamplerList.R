#' Convenience function to create an object of class mcmcSamplerList from a list of mcmc samplers
#' @author Florian Hartig
#' @param mcmcList a list with each object being an mcmcSampler
#' @return Object of class "mcmcSamplerList"
#' @export
createMcmcSamplerList <- function(mcmcList){
  # mcmcList <- list(mcmcList) -> This line didn't make any sense at all. Better would be to allow the user to simply provide several inputs without a list, but I guess the list option should be maintained, as this is convenient when scripting.
  for (i in 1:length(mcmcList)){
    if (! ("mcmcSampler" %in% class(mcmcList[[i]])) ) stop("list objects are not of class mcmcSampler")
  }
  class(mcmcList) = c("mcmcSamplerList", "bayesianOutput")
  return(mcmcList)
}

#' @author Stefan Paul
#' @method summary mcmcSamplerList
#' @export
summary.mcmcSamplerList <- function(object, ...){
  #codaChain = getSample(sampler, parametersOnly = parametersOnly, coda = T, ...)
  #summary(codaChain)
  #rejectionRate(sampler$codaChain)
  #effectiveSize(sampler$codaChain)
  #DIC(sampler)
  #max()

  sampler <- object

  DInf <- DIC(sampler)
  MAPvals <- round(MAP(sampler)$parametersMAP,3)

  gelDiag <- gelmanDiagnostics(sampler)
  psf <- round(gelDiag$psrf[,1], 3)
  
  mcmcsampler <- sampler[[1]]$settings$sampler
  
  runtime <- 0
  for(i in 1:length(sampler)) runtime <- runtime+sampler[[i]]$settings$runtime[3]

  correlations <- round(cor(getSample(sampler)),3)

  
  sampler <- getSample(sampler, parametersOnly = T, coda = T, ...)
  if("mcmc.list" %in% class(sampler)){
    nrChain <- length(sampler)
    nrIter <- nrow(sampler[[1]])
    conv <- round(gelDiag$mpsrf,3)
    npar <- ncol(sampler[[1]])
    lowerq <- upperq <- numeric(npar)
    medi <- numeric(npar)
    parnames <- colnames(sampler[[1]])
    for(i in 1:npar){
      tmp <- unlist(sampler[,i])
      tmp <- quantile(tmp, probs = c(0.025, 0.5, 0.975))
      lowerq[i] <- round(tmp[1],3)
      medi[i] <- round(tmp[2],3)
      upperq[i] <- round(tmp[3],3)
    }

  }else{
    nrChain <- 1
    nrIter <- nrow(sampler)
    npar <- ncol(sampler)
    conv <- "Only one chain; convergence cannot be determined!"
    medi <- numeric(npar)
    lowerq <- upperq <- numeric(npar)
    parnames <- colnames(sampler)
    for(i in 1:npar){
      tmp <- quantile(sampler[,i], probs = c(0.025, 0.5, 0.975))
      lowerq[i] <- round(tmp[1],3)
      medi[i] <- round(tmp[2],3)
      upperq[i] <- round(tmp[3],3)
    }

  }
  
  # output for parameter metrics
  parOutDF <- cbind(psf, MAPvals, lowerq, medi, upperq)
  colnames(parOutDF) <- c("psf", "MAP", "2.5%", "median", "97.5%")
  row.names(parOutDF) <- parnames

  
  cat(rep("#", 25), "\n")
  cat("## MCMC chain summary ##","\n")
  cat(rep("#", 25), "\n", "\n")
  cat("# MCMC sampler: ",mcmcsampler, "\n")
  cat("# Nr. Chains: ", nrChain, "\n")
  cat("# Iterations per chain: ", nrIter, "\n")
  cat("# Rejection rate: ", ifelse(object[[1]]$setup$numPars == 1, # this is a hack because coda::rejectionRate does not work for 1-d MCMC lists
                                   round(mean(sapply(sampler, coda::rejectionRate)),3), 
                                   round(mean(coda::rejectionRate(sampler)),3) ), "\n")
  cat("# Effective sample size: ", round(mean(coda::effectiveSize(sampler)),0), "\n")
  cat("# Runtime: ", runtime, " sec.","\n", "\n")
  cat("# Parameters\n")
  print(parOutDF)
  cat("\n")
  cat("## DIC: ", round(DInf$DIC,3), "\n")
  cat("## Convergence" ,"\n", "Gelman Rubin multivariate psrf: ", conv, "\n","\n")
  cat("## Correlations", "\n")
  print(correlations)
  
}

#' @author Florian Hartig
#' @method print mcmcSamplerList
#' @export
print.mcmcSamplerList <- function(x, ...){
  print("mcmcSamplerList - you can use the following methods to summarize, plot or reduce this class:")
  print(methods(class ="mcmcSamplerList"))
  #codaChain = getSample(sampler, coda = T, ...)
  #rejectionRate(sampler$codaChain)
  #effectiveSize(sampler$codaChain)
}

#' @method plot mcmcSamplerList
#' @export
plot.mcmcSamplerList <- function(x, ...){
  tracePlot(x, ...)
}

#' @author Florian Hartig
#' @export
getSample.mcmcSamplerList <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = 1, numSamples = NULL, whichParameters = NULL, includesProbabilities = F, reportDiagnostics, ...){

  
  if(!is.null(numSamples)) nS = ceiling(numSamples/length(sampler))
  
  
  # check here if due to number of chains numSamples has to be adjustet, print warning and mute internal warnings
  muteInternalGetSample = TRUE
  if(nS*length(sampler) > numSamples) {
    
    internalChains <- sampler[[1]]$chain
    if (class(internalChains)[1] == "mcmc.list") {
      nSamplesPerInternalChain <- ceiling(nS/length(internalChains))
      if (nSamplesPerInternalChain*length(internalChains) > nS) {
        warning("Due to number of external and internal chains, numSamples was rounded to the next number divisble by the number of chains.", call. = FALSE)
      }
    } else {
      warning("Due to number of external chains, numSamples was rounded to the next number divisble by the number of chains.", call. = FALSE)
    }
  }  
  
    
    
    
  
  numSamples = nS
  
    # out = NULL
    out <- list()
    for (i in 1:length(sampler)){
      # out = rbind(out, getSample(sampler[[i]], parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, numSamples = numSamples, whichParameters = whichParameters, reportDiagnostics= F))
      out[[i]] <- getSample(sampler[[i]], parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, numSamples = numSamples, whichParameters = whichParameters, reportDiagnostics= F, muteInternalGetSample=mutedInternalGetSample)
    }
    out <- combineChains(out)
  
  if(coda == T) out = makeObjectClassCodaMCMC(out, start = start, end = end, thin = thin)

  return(out)
}


