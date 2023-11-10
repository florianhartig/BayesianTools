# NOTE: The functions in this class are just templates that are to be implemented for all subclasses of BayesianOutput. They are not functional. 



#' Extracts the sample from a bayesianOutput
#' @author Florian Hartig
#' @param sampler an object of class mcmcSampler, mcmcSamplerList, smcSampler, smcSamplerList, mcmc, mcmc.list, double, numeric
#' @param parametersOnly for a BT output, if F, likelihood, posterior and prior values are also provided in the output
#' @param coda works only for mcmc classes - returns output as a coda object. Note: if mcmcSamplerList contains mcmc samplers such as DE that have several chains, the internal chains will be collapsed. This may not be desired for all applications.
#' @param start for mcmc samplers, start value in the chain. For SMC samplers, start particle
#' @param end for mcmc samplers end value in the chain. For SMC samplers, end particle
#' @param thin thinning parameter. Either an integer determining the thinning interval (default is 1) or "auto" for automatic thinning.
#' @param numSamples sample size (only used if thin = 1). If you want to use numSamples, set thin to 1.
#' @param whichParameters possibility to select parameters by index
#' @param reportDiagnostics logical, determines whether settings should be included in the output
#' @param ... further arguments
#' @example /inst/examples/getSampleHelp.R
#' @details If thin is greater than the total number of samples in the sampler object, the first and the last element (of each chain if a sampler with multiples chains is used) are sampled. If numSamples is greater than the total number of samples all samples are selected. A warning will be displayed in both cases.
#' @details If both thin and numSamples are provided, the function will use thin only if it is valid and greater than 1; otherwise, numSamples will be used.
#' @export
getSample <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = 1, numSamples = NULL, whichParameters = NULL, reportDiagnostics = FALSE, ...) UseMethod("getSample")



#' @rdname getSample
#' @author Florian Hartig
#' @export
getSample.matrix <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", numSamples = NULL, whichParameters = NULL, reportDiagnostics = F, ...){

    if(is.null(end)) end = nrow(sampler)

    out = sampler[start:end,, drop=F]

    ########################
    # THINNING
    nTotalSamples <- nrow(out)
    thin <- correctThin(nTotalSamples, thin = thin)

    if (thin == 1 && !is.null(numSamples)) {
      out <- sampleEquallySpaced(out, numSamples)
    } else {
      sel = seq(1, nTotalSamples, by = thin)
      out = out[sel,, drop=F]
    }

    if (!is.null(whichParameters)) out = out[,whichParameters, drop = FALSE]
    if(coda == T) out = makeObjectClassCodaMCMC(out, start = start, end = end, thin = thin)

    if(reportDiagnostics == T){
      return(list(chain = out, start = start, end = end, thin = thin))
    } else return(out)
}


#' @rdname getSample
#' @author Tankred Ott
#' @export
# TODO: This is right now only a helper function for getSample.mcmc. It is needed to return a vector istead of a matrix, if
#       the mcmc object passed to getSample.mcmc contains a vector.
getSample.double <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", numSamples = NULL, whichParameters = NULL, reportDiagnostics = F, ...){
  if(is.null(end)) end = length(sampler)
  out <- sampler[start:end]

  nTotalSamples <- length(out)

  thin = correctThin(nTotalSamples, thin)

  if (thin == 1 && !is.null(numSamples)) {
    out <- sampleEquallySpaced(out, numSamples)
  } else {
    sel = seq(1, nTotalSamples, by = thin)
    out = out[sel]
  }

  return(out)
}


#' @rdname getSample
#' @author Tankred Ott
#' @export
# TODO: This is right now only a helper function for getSample.mcmc. It is needed to return a vector instead of a matrix, if
#       the mcmc object passed to getSample.mcmc contains a vector.
getSample.integer <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", numSamples = NULL, whichParameters = NULL, reportDiagnostics = F, ...){
  if(is.null(end)) end = length(sampler)
  out <- sampler[start:end]

  nTotalSamples <- length(out)

  thin = correctThin(nTotalSamples, thin)

  if (thin == 1 && !is.null(numSamples)) {
    out <- sampleEquallySpaced(out, numSamples)
  } else {
    sel = seq(1, nTotalSamples, by = thin)
    out = out[sel]
  }

  return(out)
}

#' @rdname getSample
#' @author Tankred Ott
#' @export
getSample.data.frame <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", numSamples = NULL, whichParameters = NULL, reportDiagnostics = F, ...){
  getSample(as.matrix(sampler), parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, whichParameters = whichParameters, reportDiagnostics = reportDiagnostics)
}

#' @rdname getSample
#' @author Tankred Ott
#' @export
getSample.list <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", numSamples = NULL, whichParameters = NULL, reportDiagnostics = F, ...){

  if(!is.null(numSamples)) numSamples = ceiling(numSamples/length(sampler))

  if(coda == F){
    # out = NULL
    out <- rep(list(NA), length(sampler))
    for (i in 1:length(sampler)){
      # out = rbind(out, getSample(sampler[[i]], parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, numSamples = numSamples, whichParameters = whichParameters, reportDiagnostics= F))
      out[[i]] <- getSample(sampler[[i]], parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, numSamples = numSamples, whichParameters = whichParameters, reportDiagnostics= F)
    }
    out <- combineChains(out)
  }

  if(coda == T){

    out = list()

    for (i in 1:length(sampler)){

      out[[i]] = getSample(sampler[[i]], parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, numSamples = numSamples, whichParameters = whichParameters, reportDiagnostics= F)
    }

    if(inherits(out[[1]], "mcmc.list")) out = unlist(out, recursive = F)
    class(out) = "mcmc.list"
    out = out
  }

  return(out)
}

# The following two S3 implementations make getSample compatible with coda::mcmc and coda::mcmc.list

#' @rdname getSample
#' @author Tankred Ott
#' @export
getSample.mcmc <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", numSamples = NULL, whichParameters = NULL, reportDiagnostics = F, ...){

  if(coda == T){
    # mcmc objects can contain matrices or vectors
    if (is.matrix(sampler)) {
      nTotalSamples <- nrow(sampler)
    } else {
      nTotalSamples <- length(sampler)
    }

    if (is.null(end)) end = nTotalSamples

    # check/correct thin
    thin <- correctThin(nTotalSamples, thin)

    # see http://svitsrv25.epfl.ch/R-doc/library/coda/html/window.mcmc.html
    # for coda's window implementation
    return(window(sampler, start = start, end = end, thin = thin))

  } else if(coda == F){
    # mcmc objects can contain matrices or vectors
    # TODO:  do vector case as 1-d matrix?
    if (is.matrix(sampler)) {
      out <- getSample(as.matrix(sampler), parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, numSamples = numSamples, whichParameters = whichParameters, reportDiagnostics = reportDiagnostics)
    } else {
      out <- getSample(as.vector(sampler), parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, numSamples = numSamples, whichParameters = whichParameters, reportDiagnostics = reportDiagnostics)
    }
    return(out)
  }
}


#' @author Tankred Ott
#' @rdname getSample
#' @export
getSample.mcmc.list <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", numSamples = NULL, whichParameters = NULL, reportDiagnostics = F, ...){

  # TODO: implement handling of wrong inputs?

  if(coda == T){

    if (is.matrix(sampler[[1]])) {
      nTotalSamples <- nrow(sampler[[1]])
    } else {
      nTotalSamples <- length(sampler[[1]])
    }

    if (is.null(end)) end = nTotalSamples

    # check/correct thin
    thin <- correctThin(nTotalSamples, thin)

    # see http://svitsrv25.epfl.ch/R-doc/library/coda/html/window.mcmc.html
    # for coda's window implementation
    return(window(sampler, start = start, end = end, thin = thin))

  } else if(coda == F){
    if(is.matrix(sampler[[1]])) {
      return(getSample(combineChains(sampler), parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, numSamples = numSamples, whichParameters = whichParameters, reportDiagnostics = reportDiagnostics))
    } else {
      return(as.vector(getSample(combineChains(sampler), parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, numSamples = numSamples, whichParameters = whichParameters, reportDiagnostics = reportDiagnostics)))
    }
  }
}


# getSample implementation for nimble objects of class MCMC

#' @rdname getSample
#' @author Tankred Ott
#' @export
getSample.MCMC <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", numSamples = NULL, whichParameters = NULL, reportDiagnostics = F, ...){
  return(getSample(as.matrix(sampler$mvSamples), parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, numSamples = numSamples, whichParameters = whichParameters, reportDiagnostics = reportDiagnostics))
}

#' @rdname getSample
#' @author Tankred Ott
#' @export
getSample.MCMC_refClass <- function(sampler, parametersOnly = T, coda = F, start = 1, end = NULL, thin = "auto", numSamples = NULL, whichParameters = NULL, reportDiagnostics = F, ...){
  return(getSample(as.matrix(sampler$mvSamples), parametersOnly = parametersOnly, coda = coda, start = start, end = end, thin = thin, numSamples = numSamples, whichParameters = whichParameters, reportDiagnostics = reportDiagnostics))
}


#' Merge Chains
#'
#' Merge a list of outputs from MCMC / SMC samplers
#'
#' The function merges a list of outputs from MCMC / SMC samplers into a single matrix. Requirement is that the list contains classes for which the \code{\link{getSample}} function works
#'
#' @param l a list with objects that can be accessed with \code{\link{getSample}}
#' @param ... arguments to be passed on to \code{\link{getSample}}
#'
#' @return a matrix
#'
#' @author Florian Hartig
#'
#' @export
mergeChains <- function(l, ...){

  x = getSample(l[[1]], ...)

  for(i in 2:length(l)){
    x = rbind(x, getSample(l[[i]], ...))
  }
  return(x)
}
