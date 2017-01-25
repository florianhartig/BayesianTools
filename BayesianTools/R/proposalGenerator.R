#' Factory that creates a proposal generator
#' @author Florian Hartig
#' @param covariance covariance matrix. Can also be vector of the sqrt of diagonal elements --> standard deviation
#' @param gibbsProbabilities optional probabilities for the number of parameters to vary in a Metropolis within gibbs style - for 4 parameters, c(1,1,0.5,0) means that at most 3 parameters will be varied, and it is double as likely to vary one or two than varying 3 
#' @param gibbsWeights optional probabilities for parameters to be varied in a Metropolis within gibbs style - default ist equal weight for all parameters - for 4 parameters, c(1,1,1,100) would mean that if 2 parameters would be selected, parameter 4 would be 100 times more likely to be picked than the others. If 4 is selected, the remaining parameters have equal probability.
#' @param otherDistribution optional additinal distribution to be mixed with the default multivariate normal. The distribution needs to accept a parameter vector (to allow for the option of making the distribution dependend on the parameter values), but it is still assumed that the change from the current values is returned, not the new absolute values. 
#' @param otherDistributionLocation a vector with 0 and 1, denoting which parameters are modified by the otherDistribution
#' @param otherDistributionScaled should the other distribution be scaled if gibbs updates are calculated?
#' @param message print out parameter settings
#' @param method method for covariance decomposition
#' @param scalingFactor scaling factor for the proposals
#' @seealso \code{\link{updateProposalGenerator}}
#' @export 
#' @example /inst/examples/proposalGeneratorHelp.R
#'
 
#@param covarianceDecomp composed covariance matrix. If provided, faster TODO?

createProposalGenerator <- function(
  covariance, # covariance matrix for the multivariate proposal
  gibbsProbabilities = NULL, #  changes
  gibbsWeights = NULL,
  otherDistribution = NULL,
  otherDistributionLocation = NULL, 
  otherDistributionScaled = F,
  message = F,
  method = "chol",
  scalingFactor = 2.38
  ) {
  
  # To provide the option of defining via sd of individual normal
  if (is.vector(covariance)) {
    covariance = diag(covariance^2) 
  }
  
  if(ncol(covariance) == 0 && nrow(covariance) == 0) covariance = 1
  
  if(is.null(otherDistribution)) numberOfParameters = max(1,nrow(covariance)) 
  else numberOfParameters = length(otherDistributionLocation)
  
  if(is.null(method) | numberOfParameters < 2){
    covarianceDecomp = NULL
    if(numberOfParameters > 1) samplingFunction = function() as.vector(mvtnorm::rmvnorm(n = 1, sigma = covariance))
    else samplingFunction = function() rnorm(n = 1, sd = sqrt(covariance)) 
  } else {
    covarianceDecomp = factorMatrice(covariance, method = method)
    samplingFunction = function() as.vector(getRmvnorm(n = 1, R = covarianceDecomp))
  }
  
  ## Assertions

  if(!is.null(otherDistribution)){
    stopifnot(class(otherDistribution) == "function")
    stopifnot(!is.null(otherDistributionLocation))
    if(is.numeric(otherDistributionLocation)) otherDistributionLocation = as.logical(otherDistributionLocation) 
    stopifnot(is.logical(otherDistributionLocation))      
    stopifnot(is.logical(otherDistributionScaled))  
  } 
  
  #scalingFactor = 2.38/sqrt(numberOfParameters) # CHECK ???
  #scalingFactorN = (2.38^2)/numberOfParameters 
  # note - scaling is 2.38 * sqrt because it is applied on the change, not directly on the sigma

  ##########################
  # Definition of proposal function
  
  returnProposal <- function(x, scale = 1){     
    
    # Possibility to mix with other distribution
    if(!is.null(otherDistribution)){
      move = rep(NA, numberOfParameters)
      move[otherDistributionLocation] = otherDistribution(x[otherDistributionLocation])
      move[!otherDistributionLocation] = samplingFunction()
    }else{
      move = samplingFunction()
    }
    
    ## Gibbs updates
    if (!is.null(gibbsProbabilities)) {

      nGibbs <- sample.int(length(x), size = 1, replace = F, prob = gibbsProbabilities)
      whichParametersLoc <- sample.int(length(x), nGibbs, replace = F, prob = gibbsWeights)
      move[! (1:numberOfParameters %in% whichParametersLoc)] = 0
    } else {
      nGibbs = numberOfParameters
    }
    
    ### 
    
    if(!is.null(otherDistribution) & otherDistributionScaled == F){
      nGibbs = nrow(covariance) 
      move[!otherDistributionLocation] = move[!otherDistributionLocation] * scalingFactor / sqrt(nGibbs)
    }else{
      move = move * scalingFactor / sqrt(nGibbs)
    }
      
    newParam = x + move * scale
    
    return(newParam)
  }
  
  returnProposalMatrix <- function(x, scale = 1){
    
    numPar <- ncol(x)
    
    if (numPar == 1){
      out = matrix(apply(x, 1, returnProposal, scale = scale), ncol = 1) 
    } else {
      out = t(apply(x, 1, returnProposal, scale = scale)) 
    }
    return(out)
  }
  
  
  returnDensity <- function(x, y, scale = 1){
    if (!is.null(gibbsProbabilities) & !(is.null(otherDistribution)))stop("proposal density not implemented if Gibbs or other distribution is activated in the proposal. This error may appear if you have chosen both gibbs and delayed rejection in an MCMC algorith. This option is currently not implemented")
    
    sigmaDensity = scalingFactor^2 / numberOfParameters * covariance * scalingFactor^2
    if(length(sigmaDensity) > 1) dens = mvtnorm::dmvnorm(x, mean = y, sigma = sigmaDensity, log = T)
    else dens = dnorm(x, mean = y, sd = sqrt(sigmaDensity), log = T) 
    return(dens)
    
  }
  
 
  ##########################
  # Wrap up class fields
  
  classFields = list(
    covariance = covariance, 
    covarianceDecomp = covarianceDecomp,
    gibbsProbabilities = gibbsProbabilities, 
    gibbsWeights = gibbsWeights,
    otherDistribution = otherDistribution, 
    otherDistributionLocation = otherDistributionLocation, 
    otherDistributionScaled = otherDistributionScaled,
    returnProposal = returnProposal, 
    returnProposalMatrix = returnProposalMatrix, 
    returnDensity = returnDensity,
    updateProposalGenerator = updateProposalGenerator ,
    samplingFunction = samplingFunction
  )
  class(classFields) <- c("proposalGenerator")
  
  if(message == T){
    cat("Proposalgenerator created")
    print(classFields)
  }
  
  return(classFields)
}

#' @method print proposalGenerator
#' @export
print.proposalGenerator <- function(x, ...){
  
  names = names(x)
  
  for(i in 1:6){
    cat(names[i], "set to:\n ")
    print(x[[i]])
  } 
 
}


#' To update settings of an existing proposal genenerator
#' @param proposal an object of class proposalGenerator
#' @param chain a chain to create the covariance matrix from (optional)
#' @param message whether to print an updating message
#' @param eps numeric tolerance for covariance
#' @param manualScaleAdjustment optional adjustment for the covariance scale (multiplicative)
#' @details The this function can be applied in 2 ways 1) update the covariance given an MCMC chain, and 2) update the proposal generator after parameters have been changed
#' @export
updateProposalGenerator <- function(proposal,chain = NULL,  message = F, eps = 1e-10, manualScaleAdjustment = 1){
  
  if(!is.null(chain)){
    npar = ncol(chain)
    if(is.null(npar)) npar = 1
    
    if (npar > 1){
      covar = cov(chain) * manualScaleAdjustment
      covar = as.matrix(Matrix::nearPD(covar + diag(eps, npar))$mat)    
    }else{
      covar = var(chain) * manualScaleAdjustment
    }
    if(!any(is.na(covar))) proposal$covariance = covar
  }
  
  out <- createProposalGenerator(
    covariance = proposal$covariance, 
    gibbsProbabilities = proposal$gibbsProbabilities, 
    gibbsWeights = proposal$gibbsWeights,
    otherDistribution = proposal$otherDistribution, 
    otherDistributionLocation = proposal$otherDistributionLocation, 
    otherDistributionScaled = proposal$otherDistributionScaled
  )
  
  if(message == T){
    cat("Proposalgenerator settings changed")
    print(out)
  }
  return(out)
}



#' Produce multivariate normal proposal
#' @param n n
#' @param R R
#' @return X

getRmvnorm <- function(n=1, R){
  X <- matrix(rnorm(n * ncol(R)), nrow=n )%*%  R
  return(X)
}

#' factorMatrice
#' @param sigma sigma
#' @param method either "eigen", "svd" or "chol"
factorMatrice <- function(sigma, method){
  if(method == "eigen") {
    ev <- eigen(sigma, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))){
      warning("sigma is numerically not positive definite")
    }
    ## ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
    ## faster for large  nrow(sigma):
    t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
  }
  else if(method == "svd"){
    s. <- svd(sigma)
    if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))){
      warning("sigma is numerically not positive definite")
    }
    t(s.$v %*% (t(s.$u) * sqrt(s.$d)))
  }
  else if(method == "chol"){
    R <- chol(sigma, pivot = TRUE)
    R[, order(attr(R, "pivot"))]
  }
}

# adapt
#proposalGenerator$covariance = factorMatrice(proposalGenerator$covariance, method)







