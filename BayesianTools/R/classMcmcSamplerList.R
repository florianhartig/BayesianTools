#' Convenience function to create an object of class mcmcSamplerList from a list of mcmc samplers
#' @param mcmcList a list with each object being an mcmcSampler
#' @return Object of class "mcmcSamplerList"
#' @export

createMcmcSamplerList <- function(mcmcList){
  mcmcList <- list(mcmcList)
  for (i in 1:length(mcmcList)){
    if (! ("mcmcSampler" %in% class(mcmcList[[i]])) ) stop("list objects are not of class mcmcSampler")
  }
  class(mcmcList) = c("mcmcSamplerList", "bayesianOutput")
  return(mcmcList)
}

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
  MAPvals <- format(round(MAP(sampler)$parametersMAP,3), nsmall = 3)
  mcmcsampler <- sampler[[1]]$settings$sampler
  runtime <- 0
  for(i in 1:length(sampler)) runtime <- runtime+sampler[[i]]$settings$runtime[3]

  correlations <- round(cor(getSample(sampler)),3)
  
  
  sampler <- getSample(sampler, parametersOnly = T, coda = T)
  if("mcmc.list" %in% class(sampler)){
    nrChain <- length(sampler)
    nrIter <- nrow(sampler[[1]])
    conv <- round(coda::gelman.diag(sampler)$mpsrf,3)
    npar <- ncol(sampler[[1]])
    lowerq <- upperq <- numeric(npar)
    medi <- numeric(npar)
    parnames <- colnames(sampler[[1]])
    for(i in 1:npar){
      tmp <- unlist(sampler[,i])
      tmp <- quantile(tmp, probs = c(0.025, 0.5, 0.975))
      lowerq[i] <- format(round(tmp[1],3), nsmall = 3)
      medi[i] <- format(round(tmp[2],3), nsmall = 3)
      upperq[i] <- format(round(tmp[3],3), nsmall = 3)
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
      lowerq[i] <- format(round(tmp[1],3), nsmall = 3)
      medi[i] <- format(round(tmp[2],3), nsmall = 3)
      upperq[i] <- format(round(tmp[3],3), nsmall = 3)
    }
    
  }
  
  cat(rep("#", 25), "\n")
  cat("## MCMC chain summary ##","\n")
  cat(rep("#", 25), "\n", "\n")
  cat("# MCMC sampler: ",mcmcsampler, "\n")
  cat("# Nr. Chains: ", nrChain, "\n")
  cat("# Iterations per chain: ", nrIter, "\n")
  cat("# Rejection rate: ", round(mean(coda::rejectionRate(sampler)),3), "\n")
  cat("# Effective sample size: ", round(mean(coda::effectiveSize(sampler)),0), "\n")
  cat("# Runtime: ", runtime, " sec.","\n", "\n")
  cat("# Parameter", rep("", 6), "MAP", "    ", "2.5%", "  ", "median", " ", "97.5%", "\n")
  for(i in 1:npar){
    cat("# ", parnames[i],": ",rep("", max(1,10-nchar(parnames[i]))),
        MAPvals[i],rep("",max(0,8-nchar(lowerq[i]))), lowerq[i],rep("",max(0,8-nchar(medi[i]))), medi[i],
        rep("",max(0,8-nchar(upperq[i]))), upperq[i], "\n")
  }
  cat("\n")
  cat("## DIC: ", round(DInf$DIC,3), "\n")
  cat("## Convergence" ,"\n", "Gelman Rubin multivariate psrf: ", conv, "\n","\n")
  cat("## Correlations", "\n")
  print(correlations)
  
  
}

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
  tracePlot(x, parametersOnly = T, ...)
}

#' @export
getSample.mcmcSamplerList <- function(sampler, parametersOnly = T, coda = F, start = 1, 
                                      end = NULL, thin = 1, numSamples, whichParameters = NULL, reportDiagnostics,
                                      ...){
  
  if(coda == F){
    out = NULL
    for (i in 1:length(sampler)){
      out = rbind(out, getSample(sampler[[i]], parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, whichParameters = whichParameters, reportDiagnostics= F))
    }
    #out = BayesianTools:::combineChains(out)
  }  
  if(coda == T){
    
    out = list()
    
    for (i in 1:length(sampler)){
   
      out[[i]] = getSample(sampler[[i]], parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, whichParameters = whichParameters, reportDiagnostics= F)
    }
    
    if(class(out[[1]]) == "mcmc.list") out = unlist(out, recursive = F)
    class(out) = "mcmc.list" 
    out = out
  } 
    
  return(out)
}


