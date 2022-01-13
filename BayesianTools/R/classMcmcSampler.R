# Functions for class mcmcSamper

#' @author Florian Hartig
#' @export
getSample.mcmcSampler <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = 1, numSamples = NULL, whichParameters = NULL, includesProbabilities = F, reportDiagnostics= F, ...){
  
  if (class(sampler$chain)[1] == "matrix"){
    
    out <- getmcmcSamplerSampleHelper(sampler, parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, numSamples = numSamples, whichParameters = whichParameters, includesProbabilities = includesProbabilities, reportDiagnostics= reportDiagnostics, ...)
    
    # Sample size
    if(thin == 1 && !is.null(numSamples)){
      out <- sampleEquallySpaced(out, numSamples) # RB: here colnames get lost
    }

    # TODO - see matrix, need to check if both thing and numSamples is set 
    
    #############

  } 
  else if (class(sampler$chain) == "mcmc.list"){
    
    out = list()
    
    
    for (i in 1:length(sampler$chain)){
      temp = sampler
      temp$chain = sampler$chain[[i]]
      if(is.null(end)) end = nrow(temp$chain)
      
      temp <- getmcmcSamplerSampleHelper(temp, parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, numSamples = numSamples, whichParameters = whichParameters, includesProbabilities = includesProbabilities, reportDiagnostics= reportDiagnostics, ...)
      # Sample size
      if(thin == 1 && !is.null(numSamples)){
        nSamplesPerChain <- ceiling(numSamples/length(sampler$chain))
        
        if(i == 1 && nSamplesPerChain*length(sampler$chain) > numSamples){
          ## if getSample.mcmcSampler was called from getSample.mcmcSamplerList this warning is muted
          if (!hasArg(muteInternalGetSample)) {
            warning("Due to internal chains, numSamples was rounded to the next number divisble by the number of chains.", call. = FALSE)
          }
        }
        
        temp <- sampleEquallySpaced(temp, nSamplesPerChain)
      }
      
      out[[i]] = temp
      
    }
    class(out) = "mcmc.list" 
    
    #trueNumSamples <- sum(unlist(lapply(out, FUN = nrow)))
    #if (!is.null(numSamples) && trueNumSamples > numSamples) warning(paste(c("Could not draw ", numSamples, " samples due to rounding errors. Instead ", trueNumSamples," were drawn.")))
    out = combineChains(out)
  }else stop("sampler appears not to be of class mcmcSampler")
  
  
  if (!is.null(whichParameters)) out = out[,whichParameters, drop = F]
  if(coda == T) out = makeObjectClassCodaMCMC(out, start = start, end = end, thin = thin)
  
  if(reportDiagnostics == T){
    return(list(chain = out, start = start, end = end, thin = thin))
  } else return(out)
}

#' @author Robert Bosek
#' @export
getmcmcSamplerSampleHelper <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = 1, numSamples = NULL, whichParameters = NULL, includesProbabilities = F, reportDiagnostics= F, ...){
  
  if(is.null(end)) end = nrow(sampler$chain)
  
  if(parametersOnly == T) {
    out = sampler$chain[start:end,1:sampler$setup$numPars,drop=F] # RB: drop=F delete next lines
    #if(class(out)[1] == "numeric") out = as.matrix(out) # case 1 parameter
    if(!is.null(sampler$setup$names)) colnames(out) = sampler$setup$names
  }
  else {
    out = sampler$chain[start:end,,drop=F] 
    if(!is.null(sampler$setup$names)) colnames(out) = c(sampler$setup$names, "Lposterior", "Llikelihood", "Lprior")
  }
  
  ########################
  # THINNING
  
  if (thin == "auto"){
    thin = max(floor(nrow(out) / 5000),1)
  }
  if(is.null(thin) || thin == F || thin < 1 || is.nan(thin)) thin = 1
  if(thin > nrow(out)) warning("thin is greater than the total number of samples!")
  if (! thin == 1){
    sel = seq(1,dim(out)[1], by = thin )
    out = out[sel, ,drop=F]
  }
  return(out)
}



#' @method summary mcmcSampler
#' @author Stefan Paul
#' @export
summary.mcmcSampler <- function(object, ...){
  #codaChain = getSample(sampler, parametersOnly = parametersOnly, coda = T, ...)
  #summary(codaChain)
  #rejectionRate(sampler$codaChain)
  #effectiveSize(sampler$codaChain)
  #DIC(sampler)
  #max()
  #}
  
  sampler <- object
  
  try(DInf <- DIC(sampler), silent = TRUE)
  MAPvals <- round(MAP(sampler)$parametersMAP,3)
  psf <- FALSE
  
  mcmcsampler <- sampler$settings$sampler
  runtime <- sampler$settings$runtime[3]
  correlations <- round(cor(getSample(sampler)),3)
  
  chain <- getSample(sampler, parametersOnly = T, coda = T, ...)
  # chain <- getSample(sampler, parametersOnly = T, coda = T)
  if("mcmc.list" %in% class(chain)){
    psf <- TRUE
    nrChain <- length(chain)
    nrIter <- nrow(chain[[1]])
    conv <- ifelse(chain$setup$numPars > 1, round(coda::gelman.diag(chain)$mpsrf,3), round(coda::gelman.diag(chain)$mpsrf,3)$psrf[1])
    npar <- sampler$setup$numPars
    lowerq <- upperq <- numeric(npar)
    medi <- numeric(npar)
    parnames <- colnames(chain[[1]])
    
    # Shorthen parameter names
    for (i in 1:npar) {
      if (nchar(parnames[i]) > 8)
        parnames[i] <- paste(substring(parnames[i], 1, 6), "...", sep = "")
    }
    
    for (i in 1:npar) {
      tmp <- unlist(chain[, i])
      tmp <- quantile(tmp, probs = c(0.025, 0.5, 0.975))
      lowerq[i] <- round(tmp[1], 3)
      medi[i] <- round(tmp[2], 3)
      upperq[i] <- round(tmp[3], 3)
    }
    
  } else{
    nrChain <- 1
    nrIter <- nrow(chain)
    npar <- sampler$setup$numPars
    conv <- "Only one chain; convergence cannot be determined!"
    medi <- numeric(npar)
    lowerq <- upperq <- numeric(npar)
    parnames <- colnames(chain)
    for (i in 1:npar) {
      tmp <- quantile(chain[, i], probs = c(0.025, 0.5, 0.975))
      lowerq[i] <- round(tmp[1], 3)
      medi[i] <- round(tmp[2], 3)
      upperq[i] <- round(tmp[3], 3)
    }
    
  }
  
  parOutDF <- cbind(MAPvals, lowerq, medi, upperq)
  colnames(parOutDF) <- c("MAP", "2.5%", "median", "97.5%")
  if (psf == TRUE) {
    psf <- round(gelmanDiagnostics(sampler)$psrf[,1], 3)
    parOutDF <- cbind(psf, parOutDF)
  }
  row.names(parOutDF) <- parnames
  
  cat(rep("#", 25), "\n")
  cat("## MCMC chain summary ##","\n")
  cat(rep("#", 25), "\n", "\n")
  cat("# MCMC sampler: ",mcmcsampler, "\n")
  cat("# Nr. Chains: ", nrChain, "\n")
  cat("# Iterations per chain: ", nrIter, "\n")
  cat("# Rejection rate: ", ifelse(object$setup$numPars == 1 & class(chain) == "mcmc.list", # this is a hack because coda::rejectionRate does not work for 1-d MCMC lists
                                   round(mean(sapply(chain, coda::rejectionRate)),3), 
                                   round(mean(coda::rejectionRate(chain)),3) ), "\n")
  cat("# Effective sample size: ", ifelse(sampler$setup$numPars == 1, round(coda::effectiveSize(chain),0), round(mean(coda::effectiveSize(chain)),0) ) , "\n")
  cat("# Runtime: ", runtime, " sec.","\n", "\n")
  cat("# Parameters\n")
  print(parOutDF)
  cat("\n")
  
  try(cat("## DIC: ", round(DInf$DIC,3), "\n"), silent = TRUE)
  cat("## Convergence" ,"\n", "Gelman Rubin multivariate psrf: ", conv, "\n","\n")
  cat("## Correlations", "\n")
  print(correlations)
}

#' @author Florian Hartig
#' @method print mcmcSampler
#' @export
print.mcmcSampler <- function(x, ...){
  print("mcmcSampler - you can use the following methods to summarize, plot or reduce this class:")
  print(methods(class ="mcmcSampler"))
  #codaChain = getSample(sampler, coda = T, ...)
  #rejectionRate(sampler$codaChain)
  #effectiveSize(sampler$codaChain)
}

#' @author Florian Hartig
#' @method plot mcmcSampler
#' @export
plot.mcmcSampler <- function(x, ...){
  tracePlot(x, ...)
}




