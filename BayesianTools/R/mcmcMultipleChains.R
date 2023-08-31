#' Run multiple chains
#' @param bayesianSetup object of class "BayesianSetup"
#' @param settings list with settings for sampler
#' @param sampler character, either "Metropolis" or "DE"
#' @return list containing the single runs ($sampler) and the chains in a coda::mcmc.list ($mcmc.list)
#' @keywords internal
mcmcMultipleChains <- function(bayesianSetup, settings, sampler) {
  # Get number of chains
  nrChains <- settings$nrChains
  
  # Set settings$nrChains to one to avoid infinite loop
  settings$nrChains <- 1
  
  # Initialize output
  out <- list()
  out$sampler <- list()
  
  # Run sampler
  for (i in 1:nrChains) {
    out$sampler[[i]] <-
      runMCMC(bayesianSetup, sampler = sampler, settings = settings)
  }
  
  
  # Make coda::mcmc.list object
  for (i in 1:nrChains) {
    txtemp <- paste("coda::mcmc(out$sampler[[", i, "]]$chain)", sep = "")
    if (i == 1)
      tx = txtemp
    else
      tx <- paste(tx, txtemp, sep = ", ")
  }
  
  tx <- paste("coda::mcmc.list(", tx, ")", sep = "")
  out$mcmc.list <- eval(parse(text = tx))
  
  
  return(out)
}
