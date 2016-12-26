# The aim of this test is to check whether the diagnostic functions
# can be applied to the output of all samplers.
# It throws an error if an output of a diagnostoc is changed
# or the diagnotic function cannot be applied to an 
# object of class bayesianOutput.
# However, to shorten the test, only the structures of the 
# outputs are tested. For the validity of the mcmc sampling
# algorithms see the respective tests.



context("Test workflow of package")



set.seed(1)
library(BayesianTools)

iter = 1000
iterSMC = 2 # TODO implement tests for SMC 


runEverything <- function(bayesianOutput){
  
  out <- bayesianOutput
  if("mcmcSamplerList" %in% class(out)) nPars <- out[[1]]$setup$numPars
  else nPars <- out$setup$numPars
  
  # Test sample
  sam <- getSample(out, numSamples = 60)
  
  # Test diagnostics
  test_that("Diagnostic functions work for all outputs",{
    tmp <- DIC(out)
    expect_equal(names(tmp), c("DIC", "IC", "pD", "pV", "Dbar", "Dhat"))
    
    tmp <- MAP(out)
    expect_equal(names(tmp), c("parametersMAP", "valuesMAP"))
    
    tmp <- WAIC(out)
    expect_equal(names(tmp), c("WAIC1", "WAIC2", "lppd", "pWAIC1", "pWAIC2"))
    

    if(nPars > 1){
      suppressWarnings(tmp <- marginalLikelihood(out))
      expect_equal(names(tmp), c("marginalLikelihod", "ln.lik.star", "ln.pi.star", "ln.pi.hat", "method"))
      
    }
  })
  
} # runEverything




##############################################################


test_that("1-d samplers work", {
  
  
  
  ll = function(x, sum = TRUE){
    if(sum) sum(dnorm(x, log = T))
    else dnorm(x, log = T)
    
  } 
  
  setup = createBayesianSetup(ll, lower = c(-10), upper = c(10))
  
  samp = getPossibleSamplerTypes()
  
  for(i in 1:length(samp$BTname)){
    
    
    if(samp$univariatePossible[i] == T){
      if (samp$restartable[i] == T){
        settings = list(iterations = iter/2, adaptive = F, consoleUpdates = 1e+8)
        invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = setup, sampler = samp$BTname[i], settings = settings))))
        invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = out, sampler = samp$BTname[i], settings = settings))))
      }else{
        settings = list(adaptive = F, iterations = iter, consoleUpdates = 1e+8)
        invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = setup, sampler = samp$BTname[i], settings = settings))))
      } 
    }
    
    # run tests
    runEverything(out)
    
  }
})



test_that("2-D samplers with restarting works", {
  
  ll = function(x, sum = TRUE){
    if(sum) sum(dnorm(x, log = T))
    else dnorm(x, log = T)
  } 
  
  setup = createBayesianSetup(ll, lower = c(-10,-10), upper = c(10,10))
  
  
  samp = getPossibleSamplerTypes()
  
  for(i in 1:(length(samp$BTname)-1)){
    
    if (samp$restartable[i] == T){
      settings = list(iterations = iter/2, consoleUpdates = 1e+8)
      invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = setup, sampler = samp$BTname[i], settings = settings))))
      invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = out, sampler = samp$BTname[i], settings = settings))))
    }else{
      settings = list(iterations = iter, consoleUpdates = 1e+8)
      if(samp$BTname[i] == "SMC") settings = list(iterations = iterSMC, 
                                                  initialParticles = 10000, consoleUpdates = 1e+8)
      invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = setup, sampler = samp$BTname[i], settings = settings))))
    } 
    
    # run tests
    runEverything(out)
    
  }
  
})



test_that("2-d with parallel chains works", {
  
  ll = function(x, sum = TRUE) {
    if(sum) sum(dnorm(x, log = T))
    else dnorm(x, log = T)
  }
  setup = createBayesianSetup(ll, lower = c(-10,-10), upper = c(10,10), parallel = 2)
  
  samp = getPossibleSamplerTypes()
  
  for(i in 1:(length(samp$BTname)-1)){
    
    settings = list(nrChains = 3, iterations = iter, consoleUpdates = 1e+8)
    if(samp$BTname[i] == "SMC") settings = list(nrChains = 3, iterations = iterSMC, consoleUpdates = 1e+8)
    invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = setup, sampler = samp$BTname[i], settings = settings)))) 
    
    
    # run tests
    runEverything(out)
  }
  # Stop cluster
  rmBayesianSetup(setup)
  
})


test_that("2-d with parallel chains and restart works", {
  
  ll = function(x, sum = TRUE){
    if(sum)  sum(dnorm(x, log = T))
    else dnorm(x, log = T)
  }
  setup = createBayesianSetup(ll, lower = c(-10,-10), upper = c(10,10), parallel = 2)
  
  samp = getPossibleSamplerTypes()
  
  for(i in 1:(length(samp$BTname)-1)){
    
    if (samp$restartable[i] == T){
      settings = list(nrChains = 3, iterations = iter/2, consoleUpdates = 1e+8)
      if(samp$BTname[i] == "SMC") settings = list(nrChains = 3, iterations = iterSMC, consoleUpdates = 1e+8)
      invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = setup, sampler = samp$BTname[i], settings = settings))))
      invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = out, sampler = samp$BTname[i], settings = settings))))
    }else{
      settings = list(nrChains = 3, iterations = iter, consoleUpdates = 1e+8)
      invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = setup, sampler = samp$BTname[i], settings = settings))))
    } 
    
    # run tests
    runEverything(out)
    
  }
  
  # Stop cluster
  rmBayesianSetup(setup)
  
}
)



test_that("2-d with external parallelization and restart works", {
  
  # Define function
  FUN <- function(x, sum = TRUE){
   if(sum) sum(dnorm(x, log = T))
   else dnorm(x, log = T)
  }
  
  # Make cluster
  cl <- parallel::makeCluster(2)
  
  # Likelihood function
  ll <- function(pars,...){
    if(is.matrix(pars)) res = parallel::parApply(cl = cl, pars, 1, FUN, ...)
    else res = FUN(pars, ...)
    return(res)
  }
  
  setup <- createBayesianSetup(ll, lower = c(-10,-10), upper = c(10,10),
                                       parallel = "external")
  
  samp = getPossibleSamplerTypes()
  
  for(i in 1:(length(samp$BTname)-1)){
    
    if (samp$restartable[i] == T){
      settings = list(nrChains = 1, iterations = iter/2, consoleUpdates = 1e+8)
      if(samp$BTname[i] == "SMC") settings = list(nrChains = 3, iterations = iterSMC, consoleUpdates = 1e+8)
      invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = setup, sampler = samp$BTname[i], settings = settings))))
      invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = out, sampler = samp$BTname[i], settings = settings))))
    }else{
      settings = list(nrChains = 1, iterations = iter, consoleUpdates = 1e+8)
      invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = setup, sampler = samp$BTname[i], settings = settings))))
    } 
    
    # run tests
    runEverything(out)
    
  }
  
  # Stop cluster
  parallel::stopCluster(cl)
  
}
)



test_that("2-d with external parallelization, restart and multiple chains works", {
  
  # Define function
  FUN <- function(x, sum = TRUE){
    if(sum) sum(dnorm(x, log = T))
    else dnorm(x, log = T)
  }
  
  # Make cluster
  cl <- parallel::makeCluster(2)
  
  # Likelihood function
  ll <- function(pars,...){
    if(is.matrix(pars)) res = parallel::parApply(cl = cl, pars, 1, FUN, ...)
    else res = FUN(pars, ...)
    return(res)
  }
  
  setup <- createBayesianSetup(ll, lower = c(-10,-10), upper = c(10,10),
                               parallel = "external")
  
  samp = getPossibleSamplerTypes()
  
  for(i in 1:(length(samp$BTname)-1)){
    
    if (samp$restartable[i] == T){
      settings = list(nrChains = 3, iterations = iter/2, consoleUpdates = 1e+8)
      if(samp$BTname[i] == "SMC") settings = list(nrChains = 3, iterations = iterSMC, consoleUpdates = 1e+8)
      invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = setup, sampler = samp$BTname[i], settings = settings))))
      invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = out, sampler = samp$BTname[i], settings = settings))))
    }else{
      settings = list(nrChains = 3, iterations = iter, consoleUpdates = 1e+8)
      invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = setup, sampler = samp$BTname[i], settings = settings))))
    } 
    
    # run tests
    runEverything(out)
    
  }
  # Stop cluster
  parallel::stopCluster(cl)
}
)


