# S3 Functions for class 'smcSampler'



#' @rdname getSample
#' @author Florian Hartig
#' @export
getSample.smcSampler <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = 1, numSamples = NULL, whichParameters = NULL, reportDiagnostics = FALSE, ...){
  
  if(is.null(end)) end = nrow(sampler$particles)
  
  if(parametersOnly == T) {
    out = sampler$particles[start:end,] 
    if(!is.null(sampler$setup$names)) colnames(out) = sampler$setup$names
  }
  else {
    out = cbind(sampler$particles[start:end,] , sampler$posterior[start:end,] )
    if(!is.null(sampler$setup$names)) colnames(out) = c(sampler$setup$names, "Lposterior", "Llikelihood", "Lprior")
  }
  
  ########################
  # THINNING
  if (thin == "auto"){
    thin = max(floor(nrow(out) / 5000),1)
  }
  if(is.null(thin) || thin == F || thin < 1) thin = 1
  if (! thin == 1){
    sel = seq(1,dim(out)[1], by = thin )
    out = out[sel,]
  }
  # Sample size
  if(thin == 1 && !is.null(numSamples)){
    if (numSamples > nrow(out)) {
      numSamples = nrow(out)
      warning("numSamples is greater than the total number of samples! All samples were selected.")
    }
    if (numSamples < 1) numSamples = 1;
    sel <- seq(1,dim(out)[1], len = numSamples)
    out <- out[sel,] 
  }
  
  #############
  
  if (!is.null(whichParameters)) out = out[,whichParameters]
  
  if(reportDiagnostics == T){
    return(list(chain = out, start = start, end = end, thin = thin))
  } else return(out)
}

#' Summary for class 'smcSampler'
#' @description
#' Creates a summary table of a 'smcSampler' output
#' 
#' @author Florian Hartig
#' @method summary smcSampler
#' @describeIn summary.mcmcSampler Summary for smcSampler objects
#' @export
summary.smcSampler<- function(object, ...){
  sampler <- object
  print("SMC sampler output")
  summary(getSample(sampler, ...))
}

#' Plots of smcSampler output
#' @description
#' Plots smcSampler output
#' 
#' @method plot smcSampler
#' @export
plot.smcSampler<- function(x, ...){
  marginalPlot(x, ...)
}

#' Print of smcSampler output
#' @description
#' Print smcSampler output
#' @author Florian Hartig
#' @method print smcSampler
#' @export
print.smcSampler <- function(x, ...){
  print("smcSampler - you can use the following methods to summarize, plot or reduce this class:")
  print(methods(class ="smcSampler"))
}


