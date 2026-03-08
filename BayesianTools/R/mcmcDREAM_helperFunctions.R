##' Generates matrix of CR values based on pCR
##' @param pCR vector of crossover probabilities. Needs to be of length nCR.
##' @param settings list of settings
##' @param Npop number of chains
##' @return Matrix with CR values
#' @keywords internal
generateCRvalues <- function(pCR,settings, Npop){
  
  # Random vector, add zero to get first position 
  RandomVec <- c(0,cumsum(as.numeric(rmultinom(1, size = Npop*settings$updateInterval, prob = pCR))))
  
  # get candidate points
  cand <- sample(Npop*settings$updateInterval)
  CR <- rep(NA, Npop*settings$updateInterval)
  
  ## Now loop over chains to generate CR values
  for(i in 1:settings$nCR){
    #Start and End
    Start <- RandomVec[i]+1
    End <- RandomVec[i+1]
    
    # get candidates
    candx <- cand[Start:End]
    
    # Assign these indices settings$CR
    CR[candx] <- i/settings$nCR
  }
  ## Reshape CR
  CR <- matrix(CR,Npop,settings$updateInterval)
  
  return(CR)
} 


#' Adapts pCR values
#' @param CR vector of crossover probabilities. Needs to be of length nCR.
#' @param settings list of settings
#' @param delta vector with differences
#' @param lCR values to weight delta
#' @param Npop number of chains.
#' @return Matrix with CR values
#' @keywords internal
AdaptpCR <- function(CR, delta ,lCR, settings, Npop){
  if(any(delta >0)){  ## Adaptions can only be made if there are changes in X
    
    # Change CR to vector
    CR <- c(CR)
    
    # Store old lCR values
    lCROld <- lCR
    ## Determine lCR
    lCR <- rep(NA,settings$nCR)
    
    for (k in 1:settings$nCR){
      
      ##  how many times a CR value is used. This is used to weight delta
      CR_counter <- length(which(CR==k/settings$nCR))
      lCR[k] <- lCROld[k]+ CR_counter
    }                                     
    
    ## Adapt pCR 
    pCR <- Npop * (delta / lCR) / sum(delta)
    
    pCR[which(is.nan(pCR))] <- 1/settings$nCR # catch possible error if delta and lCR = 0
    
    ## Normalize values
    pCR <- pCR/sum(pCR)
    
  }
  return(list(pCR=pCR,lCR=lCR))
} ##AdaptpCR
