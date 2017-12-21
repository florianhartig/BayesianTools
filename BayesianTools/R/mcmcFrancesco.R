#' The Metropolis Algorithm
#' @author Francesco Minunno
#' @description The Metropolis Algorithm (Metropolis et al. 1953)
#' @param startValue vector with the start values for the algorithm. Can be NULL if FUN is of class BayesianSetup. In this case startValues are sampled from the prior.
#' @param iterations iterations to run
#' @param nBI number of burnin
#' @param parmin minimum values for the parameter vector or NULL if FUN is of class BayesianSetup
#' @param parmax maximum values for the parameter vector or NULL if FUN is of class BayesianSetup
#' @param f scaling factor
#' @param FUN function to be sampled from or object of class bayesianSetup
#' @param consoleUpdates interger, determines the frequency with which sampler progress is printed to the console
#' @references Metropolis, Nicholas, et al. "Equation of state calculations by fast computing machines." The journal of chemical physics 21.6 (1953): 1087-1092.
#' @keywords internal
# #' @export
M <- function(startValue = NULL, iterations  = 10000, nBI = 0 , parmin = NULL, parmax= NULL, f = 1, FUN, consoleUpdates=1000) {
  
  
  if(class(FUN) == "BayesianSetup"){
     if(FUN$numPars==1) stop("Sampler cannot be started for 1 parameter")
      
    if(is.null(startValue)){
      startValue <- FUN$prior$sampler()
    }
    parmin <- FUN$prior$lower
    parmax <- FUN$prior$upper
    FUN <- FUN$posterior$density
  }
  
  
  pValues = startValue
  lChain = iterations
  
  npar <- length(pValues)
  logMAP <- -Inf
  pChain <- matrix(NA_real_, nrow = lChain - nBI, ncol = npar+3)

#********************************************************************************

# First call to the model. Calculate likelihood and prior
  postL0 <- FUN(pValues, returnAll = T)
  accept.prob <- 0

#********************************************************************************

# Define Variance-covariance matrix (vcovProp) for proposal generation an

  scalProp <- f * 2.4^2/npar # This f is the scaling factor tuned manually
  covPar <- scalProp * diag((0.01 * (parmax - parmin))^2)

#********************************************************************************
# Build up the chain. Candidates for the parameter values (candidatepValues)
# are assumed to stem from a multivariate normal distribution (mvrnorm) with mean
# at the current state and covariance given by scalProp*covPar.
#-----

  for (j in 1:lChain) {
    if (j%%consoleUpdates == 0) print(c(j,postL1[1]))
    candidatepValues <- mvtnorm::rmvnorm(1, pValues, covPar)

    # Call the model and calculate the likelihood
    postL1 <- FUN(candidatepValues, returnAll = T)

    # Check whether the candidates are accepted.
    alpha <- min(exp(postL1[1] - postL0[1]), 1)
    accept <- 0
    if (runif(1) < alpha) {
      postL0 <- postL1
      pValues <- candidatepValues
		  accept <-  1
		  if (postL0[1] > logMAP)
		  {
		    logMAP <- postL0[1]
		    psetMAP <- pValues
		  }
		  
    }
    if (j > nBI) {
      pChain[j-nBI,] <- c(pValues,postL0)
      accept.prob <- accept.prob + accept
    }
  }
  accept.prob <- accept.prob/(lChain-nBI)
  list(Draws = pChain, accept.prob = accept.prob,psetMAP=psetMAP)
}


#' The Adaptive Metropolis Algorithm
#' @author Francesco Minunno
#' @description The Adaptive Metropolis Algorithm (Haario et al. 2001)
#' @param startValue vector with the start values for the algorithm. Can be NULL if FUN is of class BayesianSetup. In this case startValues are sampled from the prior.
#' @param iterations iterations to run
#' @param nBI number of burnin
#' @param parmin minimum values for the parameter vector or NULL if FUN is of class BayesianSetup
#' @param parmax maximum values for the parameter vector or NULL if FUN is of class BayesianSetup
#' @param f scaling factor
#' @param FUN function to be sampled from or object of class bayesianSetup
#' @param eps small number to avoid singularity
#' @references  Haario, Heikki, Eero Saksman, and Johanna Tamminen. "An adaptive Metropolis algorithm." Bernoulli (2001): 223-242.
#' @keywords internal
# #' @export
AM <- function(startValue = NULL, iterations = 10000, nBI = 0, parmin = NULL, parmax = NULL, FUN, f = 1, eps = 0) {
  
  if(class(FUN) == "BayesianSetup"){
    if(FUN$numPars==1) stop("Sampler cannot be started for 1 parameter")
    if(is.null(startValue)){
      startValue <- FUN$prior$sampler()
    }
    parmin <- FUN$prior$lower
    parmax <- FUN$prior$upper
    FUN <- FUN$posterior$density
  }
  
  
  
  pValues = startValue
  lChain = iterations
  
  noAdapt <- 1000
  n.iter <- lChain + noAdapt
  npar = length(pValues)
  pChain <- matrix(NA_real_, nrow = n.iter - nBI, ncol = npar+3)
  
  #********************************************************************************
  
  # First call to the model. Calculate likelihood and prior
  postL0 <- FUN(pValues, returnAll = T)  
  accept.prob <- 0
  
  epsDiag <- eps * diag(npar)
  scalProp <- f * (2.4^2/npar)
  covPar <- scalProp * diag((0.01*(parmax - parmin))^2)
  
  for (j in 1:n.iter) {
    candidatepValues <- as.vector(mvtnorm::rmvnorm(1, pValues, covPar))
    
    postL1 <- FUN(candidatepValues, returnAll = T)

    alpha <- min(exp(postL1[1] - postL0[1]), 1)
    accept <- 0
    if (runif(1) < alpha) {
      postL0 <- postL1
      pValues <- candidatepValues
      accept <- 1
    }
    
    if (j > nBI) {
      pChain[j-nBI,] <- c(pValues, postL0)
    }

    if (j == (nBI + noAdapt)) {
      avePar <- apply(pChain[1:noAdapt,1:npar], 2, mean)
      covPar <- scalProp * (cov(pChain[1:noAdapt,1:npar], pChain[1:noAdapt,1:npar]) + epsDiag)
    }
    if (j > (nBI + noAdapt)) {
      accept.prob <- accept.prob + accept
      t <- j - nBI
      avePar_new <- as.vector(((t-1) * avePar + pValues) / t)
      covPar_new <- ((t-2) * covPar + scalProp * ((t-1) * (avePar %o% avePar) - t * (avePar_new %o% avePar_new) + (pValues %o% pValues)) + epsDiag) / (t-1)
      avePar <- avePar_new
      covPar <- covPar_new
    }
  }
  accept.prob = accept.prob/(lChain-nBI)
  list(Draws = pChain[(noAdapt+1):(n.iter-nBI),], accept.prob = accept.prob)
}
                     

#' The Delayed Rejection Algorithm 
#' @author Francesco Minunno
#' @description The Delayed Rejection Algorithm (Tierney and Mira, 1999)
#' @param startValue vector with the start values for the algorithm. Can be NULL if FUN is of class BayesianSetup. In this case startValues are sampled from the prior.
#' @param iterations iterations to run
#' @param nBI number of burnin
#' @param parmin minimum values for the parameter vector or NULL if FUN is of class BayesianSetup
#' @param parmax maximum values for the parameter vector or NULL if FUN is of class BayesianSetup
#' @param f1 scaling factor for first proposal
#' @param f2 scaling factor for second proposal
#' @param FUN function to be sampled from or object of class bayesianSetup
#' @references Tierney, Luke, and Antonietta Mira. "Some adaptive Monte Carlo methods for Bayesian inference." Statistics in medicine 18.1718 (1999): 2507-2515.
#' @keywords internal
# #' @export
DR <- function(startValue = NULL, iterations = 10000, nBI=0, parmin = NULL, parmax =NULL, f1 = 1, f2= 0.5, FUN) {
  
  if(class(FUN) == "BayesianSetup"){
    if(FUN$numPars==1) stop("Sampler cannot be started for 1 parameter")
    if(is.null(startValue)){
      startValue <- FUN$prior$sampler()
    }
    parmin <- FUN$prior$lower
    parmax <- FUN$prior$upper
    FUN <- FUN$posterior$density
  }
  
  pValues = startValue
  lChain = iterations
  
  npar = length(pValues)
  pChain <- matrix(NA_real_, nrow = lChain - nBI, ncol = npar+3)
  
  #********************************************************************************
  
  # First call to the model. Calculate likelihood and prior
  postL0 <- FUN(pValues, returnAll = T)
  
  #********************************************************************************
  # Define Variance-covariance matrix (vcovProp) for proposal generation an  
  covPar <- diag((0.01 * (parmax - parmin))^2)
  sP <- (2.4^2/npar) * c(f1, f2)
  accept.prob <- 0
  
  for (j in 1:lChain) {
    candidatepValues <- mvtnorm::rmvnorm(1, pValues, sP[1] * covPar)
    
    # Call the model and calculate the likelihood
    postL1 <- FUN(candidatepValues, returnAll = T)
    
    # Check whether the candidates are accepted. If yes and if burn-in has been completed,
    alpha1 <- min(exp(postL1[1]-postL0[1]), 1.0)
    accept <- 0
    if (runif(1) < alpha1) {
      pValues <- candidatepValues
      postL0 = postL1
      accept <- 1
    } else {        
      candidatepValues2 <- mvtnorm::rmvnorm(1, pValues, sP[2] * covPar)
        
      # Call the model and calculate the likelihood
      postL2 <- FUN(candidatepValues2, returnAll = T)
      
      # Check whether the candidates are accepted. 

      alpha2 <- min(exp(postL1[1]-postL2[1]), 1.0)
      temp <- mvtnorm::dmvnorm(candidatepValues, candidatepValues2, sP[1] * covPar) / mvtnorm::dmvnorm(candidatepValues, pValues, sP[1] * covPar)
      alpha <- min(exp(postL2[1]-postL0[1]) * temp * ((1.0-alpha2)/(1.0-alpha1)), 1.0)
      if(is.nan(alpha)) {
        alpha <- -1
      } 
      if (runif(1) < alpha) {
        pValues <- candidatepValues2
        postL0 <- postL2                 
        accept <-  1
      }
    }
    if (j > nBI) {
      pChain[j-nBI,] <- c(pValues, postL0)
      accept.prob <- accept.prob + accept
    }
  }
  accept.prob = accept.prob/(lChain-nBI)
  list(Draws = pChain, accept.prob = accept.prob)
}

        
#' The Delayed Rejection Adaptive Metropolis Algorithm 
#' @author Francesco Minunno
#' @description The Delayed Rejection Adaptive Metropolis Algorithm (Haario et al. 2001)
#' @param startValue vector with the start values for the algorithm. Can be NULL if FUN is of class BayesianSetup. In this case startValues are sampled from the prior.
#' @param iterations iterations to run
#' @param nBI number of burnin
#' @param parmin minimum values for the parameter vector or NULL if FUN is of class BayesianSetup
#' @param parmax maximum values for the parameter vector or NULL if FUN is of class BayesianSetup
#' @param f scaling factor
#' @param FUN function to be sampled from
#' @param eps small number to avoid singularity or object of class bayesianSetup
#' @references  Haario, Heikki, Eero Saksman, and Johanna Tamminen. "An adaptive Metropolis algorithm." Bernoulli (2001): 223-242.
#' @keywords internal
# #' @export
DRAM <- function(startValue = NULL, iterations = 10000, nBI = 0, parmin = NULL, parmax = NULL, FUN, f = 1, eps = 0) {
  
  if(class(FUN) == "BayesianSetup"){
    if(FUN$numPars==1) stop("Sampler cannot be started for 1 parameter")
    if(is.null(startValue)){
      startValue <- FUN$prior$sampler()
    }
    parmin <- FUN$prior$lower
    parmax <- FUN$prior$upper
    FUN <- FUN$posterior$density
  }
  
  pValues = startValue
  lChain = iterations
  
  noAdapt <- 1000
  n.iter <- lChain + noAdapt
  npar = length(pValues)
  pChain <- matrix(NA_real_, nrow = n.iter - nBI, ncol = npar+3)
  
  #********************************************************************************
  
  # First call to the model. Calculate likelihood and prior
  postL0 <- FUN(pValues, returnAll = T)  
  accept.prob <- 0
  
  epsDiag <- eps * diag(npar)
  scalProp <- f * (2.4^2/npar)
  covPar <- scalProp * diag((0.01*(parmax - parmin))^2)
  
  for (j in 1:n.iter) {
    candidatepValues <- as.vector(mvtnorm::rmvnorm(1, pValues, covPar))
    
    postL1 <- FUN(candidatepValues, returnAll = T)
    
    alpha1 <- min(exp(postL1[1] - postL0[1]), 1)
    accept <- 0
    if (runif(1) < alpha1) {
      postL0 <- postL1
      pValues <- candidatepValues
      accept <-  1
    } else {        
      candidatepValues2 <- as.vector(mvtnorm::rmvnorm(1, pValues, 0.5 * covPar))
      
      # Call the model and calculate the likelihood
      postL2 <- FUN(candidatepValues2, returnAll = T)
      
      # Check whether the candidates are accepted. 
      
      alpha2 <- min(exp(postL1[1]-postL2[1]), 1.0)
      temp <- mvtnorm::dmvnorm(candidatepValues, candidatepValues2, covPar) / mvtnorm::dmvnorm(candidatepValues, pValues, covPar)
      alpha <- min(exp(postL2[1]-postL0[1]) * temp * ((1.0-alpha2)/(1.0-alpha1)), 1.0)
      if(is.nan(alpha)) {
        alpha <- -1
      } 
      if (runif(1) < alpha) {
        pValues <- candidatepValues2
        postL0 <- postL2                 
        accept <-  1
      }
    }
    
    if (j > nBI) {
      pChain[j-nBI,] <- c(pValues, postL0)
    }
    
    if (j == (nBI + noAdapt)) {
      avePar <- apply(pChain[1:noAdapt,1:npar], 2, mean)
      covPar <- scalProp * (cov(pChain[1:noAdapt,1:npar], pChain[1:noAdapt,1:npar]) + epsDiag)
    }
    if (j > (nBI + noAdapt)) {
      accept.prob <- accept.prob + accept
      t <- j - nBI
      avePar_new <- as.vector(((t-1) * avePar + pValues) / t)
      covPar_new <- ((t-2) * covPar + scalProp * ((t-1) * (avePar %o% avePar) - t * (avePar_new %o% avePar_new) + (pValues %o% pValues)) + epsDiag) / (t-1)
      avePar <- avePar_new
      covPar <- covPar_new
    }
  }
  accept.prob = accept.prob/(lChain-nBI)
  list(Draws = pChain[(noAdapt+1):(n.iter-nBI),], accept.prob = accept.prob)
}
