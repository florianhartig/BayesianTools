# Functions for class mcmcSamper

#' @export
getSample.mcmcSampler <- function(sampler, parametersOnly = T, coda = F, 
                                  start = 1, end = NULL, thin = 1, numSamples = NULL, whichParameters = NULL, 
                                  reportDiagnostics= F){
  
  if (class(sampler$chain) == "matrix"){
    
    if(is.null(end)) end = nrow(sampler$chain)
    
    if(parametersOnly == T) {
      out = sampler$chain[start:end,1:sampler$setup$numPars] 
      if(class(out) == "numeric") out = as.matrix(out) # case 1 parameter
      if(!is.null(sampler$setup$names)) colnames(out) = sampler$setup$names
    }
    else {
      out = sampler$chain[start:end,] 
      if(!is.null(sampler$setup$names)) colnames(out) = c(sampler$setup$names, "Lposterior", "Llikelihood", "Lprior")
    }
    
    ########################
    # THINNING
    if (thin == "auto"){
      thin = max(floor(nrow(out) / 5000),1)
    }
    if(is.null(thin) | thin == F) thin = 1
    if (! thin == 1){
      sel = seq(1,dim(out)[1], by = thin )
      out = out[sel,]
    }
    # Sample size
    if(thin == 1 && !is.null(numSamples)){
      sel <- seq(1,dim(out)[1], len = numSamples)
      out <- out[sel,] 
    }
    
    #############
    
    if (!is.null(whichParameters)) out = out[,whichParameters]
    if(coda == T) out = makeObjectClassCodaMCMC(out, start = start, end = end, thin = thin)
  } 
  else if (class(sampler$chain) == "mcmc.list"){
    
    out = list()
    
    for (i in 1:length(sampler$chain)){
      
      if(is.null(end)) end = nrow(sampler$chain[[1]])
      
      temp = sampler$chain[[i]][start:end,]
      
      if(parametersOnly == T) {
        temp = temp[,1:sampler$setup$numPars] 
        if(class(temp) == "numeric") temp = as.matrix(temp) # case 1 parameter
        if(!is.null(sampler$setup$names)) colnames(temp) = sampler$setup$names
      }
      else {
        if(!is.null(sampler$setup$names)) colnames(temp) = c(sampler$setup$names, "Lposterior", "Llikelihod", "Lprior")
      }
      
      ########################
      # THINNING
      if (thin == "auto"){
        thin = max(floor(nrow(temp) / 5000),1)
      }
      if(is.null(thin) | thin == F) thin = 1
      if (! thin == 1){
        sel = seq(1,dim(temp)[1], by = thin )
        temp = temp[sel,]
      }
      
      # Sample size
      if(thin == 1 && !is.null(numSamples)){
        sel <- seq(1,dim(temp)[1], len = (floor(numSamples/length(sampler$chain))))
        temp <- temp[sel,] 
      }
      
      
      #############
      
      if (!is.null(whichParameters)) temp = temp[,whichParameters]
      
      out[[i]] = coda::mcmc(temp)
    }
    class(out) = "mcmc.list" 
    
    if(coda == F){
      out = combineChains(out)
    }  
    if(coda == T){
      out = out
    } 
  }else stop("sampler appears not to be of class mcmcSampler")
  
  if(reportDiagnostics == T){
    return(list(chain = out, start = start, end = end, thin = thin))
  } else return(out)
}



#' @method summary mcmcSampler
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
  MAPvals <- format(round(MAP(sampler)$parametersMAP,3), nsmall = 3)
  mcmcsampler <- sampler$settings$sampler
  runtime <- sampler$settings$runtime[3]
  correlations <- round(cor(getSample(sampler)),3)
  
  
  chain <- getSample(sampler, parametersOnly = T, coda = T, ...)
  # chain <- getSample(sampler, parametersOnly = T, coda = T)
  if("mcmc.list" %in% class(chain)){
   nrChain <- length(chain)
   nrIter <- nrow(chain[[1]])
   conv <- ifelse(chain$setup$numPars > 1, round(coda::gelman.diag(chain)$mpsrf,3), round(coda::gelman.diag(chain)$mpsrf,3)$psrf[1])
   npar <- sampler$setup$numPars
   lowerq <- upperq <- numeric(npar)
   medi <- numeric(npar)
   parnames <- colnames(chain[[1]])
   
   # Shorthen parameter names
   for(i in 1:npar){
     if(nchar(parnames[i])>8) parnames[i] <- paste(substring(parnames[i], 1,6), "...", sep = "")
   }
   
     for(i in 1:npar){
      tmp <- unlist(chain[,i])
      tmp <- quantile(tmp, probs = c(0.025, 0.5, 0.975))
      lowerq[i] <- format(round(tmp[1],3), nsmall = 3)
      medi[i] <- format(round(tmp[2],3), nsmall = 3)
      upperq[i] <- format(round(tmp[3],3), nsmall = 3)
     }
  
  }else{
   nrChain <- 1
   nrIter <- nrow(chain)
   npar <- sampler$setup$numPars
   conv <- "Only one chain; convergence cannot be determined!"
   medi <- numeric(npar)
   lowerq <- upperq <- numeric(npar)
   parnames <- colnames(chain)
   for(i in 1:npar){
     tmp <- quantile(chain[,i], probs = c(0.025, 0.5, 0.975))
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
  cat("# Rejection rate: ", ifelse(sampler$setup$numPars == 1, round(mean(sapply(chain, coda::rejectionRate)),3), round(mean(coda::rejectionRate(chain)),3) ) , "\n")
  cat("# Effective sample size: ", ifelse(sampler$setup$numPars == 1, round(coda::effectiveSize(chain),0), round(mean(coda::effectiveSize(chain)),0) ) , "\n")
  cat("# Runtime: ", runtime, " sec.","\n", "\n")
  cat("# Parameter", rep("", 6), "MAP", "    ", "2.5%", "  ", "median", " ", "97.5%", "\n")
  
  space <- numeric(4)
  for(i in 1:npar){
    
    space[1] <- try(max(1,10-nchar(parnames[i])), silent = T)
    space[2] <- try(8-nchar(lowerq[i]), silent = T)
    space[3] <- try(8-nchar(medi[i]), silent = T)
    space[4] <- try(8-nchar(upperq[i]), silent = T)
    
    for(k in 1:4){
      if(class(space[k]) == "character" | space[k] < 1) space[k] <- 2
    }
    
    cat("# ", parnames[i],": ",rep("", space[1]),
        MAPvals[i],rep("",space[2]), lowerq[i],rep("", space[3]), medi[i],
        rep("", space[4]), upperq[i], "\n")
  }
  cat("\n")
  try(cat("## DIC: ", round(DInf$DIC,3), "\n"), silent = TRUE)
  cat("## Convergence" ,"\n", "Gelman Rubin multivariate psrf: ", conv, "\n","\n")
  cat("## Correlations", "\n")
  print(correlations)
}


#' @method print mcmcSampler
#' @export
print.mcmcSampler <- function(x, ...){
  print("mcmcSampler - you can use the following methods to summarize, plot or reduce this class:")
  print(methods(class ="mcmcSampler"))
  #codaChain = getSample(sampler, coda = T, ...)
  #rejectionRate(sampler$codaChain)
  #effectiveSize(sampler$codaChain)
}

#' @method plot mcmcSampler
#' @export
plot.mcmcSampler <- function(x, ...){
  tracePlot(x, ...)
}




