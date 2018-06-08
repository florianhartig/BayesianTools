
context("Test basic functionality of all samplers")


set.seed(1)
library(BayesianTools)


test <- "CRAN" # avoids exact tests
# test <- "exact"
# library(Matching)

if(test == "exact"){
  iter = 500000
  start = 500
  iterSMC = 400
  library(msm) 
  # library(Matching) # Because this is only needed for manual testing I removed it from
  # the package dependencies. If you are running the manual tests please install
  # the package yourself.


test_that("1-d samplers work", {
  
  
  if(test == "exact") skip_on_cran()
  
  ll = function(x) sum(dnorm(x, log = T))
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

        ks <- ks.boot(getSample(out, numSamples = 10000), rnorm(10000))$ks.boot.pvalue  
        
        # Test that distribution is not significally different from gaussian
        expect_true(ks>0.05)
      
      
    }
    
    
  }
  
  
})



test_that("2-D samplers with restarting works", {
  
  ll = function(x) sum(dnorm(x, log = T))
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

      # Get sample 
      x = getSample(out, numSamples  = 10000)
      y <- rnorm(10000)
      for(z in 1:ncol(x)){
        
        #  ks <- ks.test(x[,z], pnorm)$p.value
        
        ks <- ks.boot(x[,z], y)$ks.boot.pvalue 
        
        # Test that distribution is not significally different from gaussian
        expect_true(ks>0.05)
        
      }
    
  }
  
})



test_that("2-d with parallel chains works", {
  
  
  ll = function(x) sum(dnorm(x, log = T))
  setup = createBayesianSetup(ll, lower = c(-10,-10), upper = c(10,10), parallel = 2)
  
  
  samp = getPossibleSamplerTypes()
  
  for(i in 1:(length(samp$BTname)-1)){

    settings = list(nrChains = 3, iterations = iter, consoleUpdates = 1e+8)
    if(samp$BTname[i] == "SMC") settings = list(nrChains = 3, iterations = iterSMC, consoleUpdates = 1e+8)
    invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = setup, sampler = samp$BTname[i], settings = settings)))) 


      # Get sample 
      x = getSample(out, numSamples  = 10000)
      y <- rnorm(10000)
      
      for(z in 1:ncol(x)){
        
        #ks <- ks.test(x[,z],pnorm)$p.value  
        ks <- ks.boot(x[,z], y)$ks.boot.pvalue 
        # Test that distribution is not significally different from gaussian
        expect_true(ks>0.05)
        
      }
    
    
  }
})


test_that("2-d with parallel chains and restart works", {
  
  ll = function(x) sum(dnorm(x, log = T))
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
    
      # Get sample and remove first 1000 iterations as burn-in
      x = getSample(out, numSamples  = 10000)
      y <- rnorm(10000)
      
      for(z in 1:ncol(x)){
        
        # ks <- ks.test(x[,z], pnorm)$p.value
        
        ks <- ks.boot(x[,z], y)$ks.boot.pvalue 
        
        # Test that distribution is not significally different from gaussian
        expect_true(ks>0.05)
        
      }
    
  }
  
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
    
    # Get sample and remove first 1000 iterations as burn-in
    x = getSample(out, numSamples  = 10000)
    y <- rnorm(10000)
    
    for(z in 1:ncol(x)){
      
      # ks <- ks.test(x[,z], pnorm)$p.value
      
      ks <- ks.boot(x[,z], y)$ks.boot.pvalue 
      
      # Test that distribution is not significally different from gaussian
      expect_true(ks>0.05)
      
    }
   
    
  }
  
  # Stop cluster
  parallel::stopCluster(cl)
  
}
)

} # if test == "exact"
