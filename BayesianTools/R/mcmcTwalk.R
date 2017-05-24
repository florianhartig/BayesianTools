#' T-walk MCMC
#' @author Stefan Paul
#' @param bayesianSetup Object of class 'bayesianSetup' or 'bayesianOuput'.
#' @param settings  list with parameter values. 
#' @param iterations Number of model evaluations
#' @param at "traverse" move proposal parameter. Default to 6
#' @param aw "walk" move proposal parameter. Default to 1.5
#' @param pn1 Probability determining the number of parameters that are changed
#' @param Ptrav Move probability of "traverse" moves, default to 0.4918
#' @param Pwalk Move probability of "walk" moves, default to 0.4918
#' @param Pblow Move probability of "traverse" moves, default to 0.0082
#' @param burnin number of iterations treated as burn-in. These iterations are not recorded in the chain.
#' @param thin thinning parameter. Determines the interval in which values are recorded.
#' @param startValue Matrix with start values
#' @param consoleUpdates Intervall in which the sampling progress is printed to the console
#' @param message logical determines whether the sampler's progress should be printed
#' @details  
##' The probability of "hop" moves is 1 minus the sum of all other probabilities.
#' @return Object of class bayesianOutput.
#' @references Christen, J. Andres, and Colin Fox. "A general purpose sampling algorithm for continuous distributions (the t-walk)." Bayesian Analysis 5.2 (2010): 263-281.
#' @export
Twalk <- function (bayesianSetup, settings = list(iterations = 10000, at = 6, aw = 1.5, 
                                                  pn1 = NULL, Ptrav = 0.4918, Pwalk = 0.4918, 
                                                  Pblow = 0.0082, burnin = 0, thin= 1, startValue = NULL, consoleUpdates = 100,
                                                  message = TRUE)) 
{
  if("bayesianOutput" %in% class(bayesianSetup)){
    restart <- TRUE
    setup <- bayesianSetup$setup
  }else{
    restart <- FALSE
    setup <- bayesianSetup
  }
  
  setup <- checkBayesianSetup(setup, parallel = settings$parallel) # calling parallel will check if requested parallelization in settings is provided by the BayesianSetup
  if(is.null(settings$parallel)) settings$parallel = setup$parallel # checking back - if no parallelization is provided, we use the parallelization in the BayesianSetup. We could also set parallel = F, but I feel it makes more sense to use the Bayesiansetup as default
  
  
  aw <- settings$aw
  at <- settings$at
  Npar <- setup$numPars
  iterations <- floor(settings$iterations/2) # Divided by 2 because two chains are run
  if(is.null(settings$pn1)) pn1 <- min(Npar,4)/Npar
  else pn1 <- settings$pn1
  Ptrav <- settings$Ptrav
  if(is.null(settings$Pwalk)) Pwalk <-  0.4918
  else Pwalk <- settings$Pwalk
  if(is.null(settings$Pblow)) Pblow <- 0.0082
  else Pblow <- settings$Pblow
  
  
  # Set burnin and thin
  burnin <- settings$burnin
  thin <- settings$thin
  
  # Set Phop
  Phop <- 1-(Ptrav+Pwalk+Pblow)
  
  # Check for consistency of move probabilities
  if((Pwalk + Ptrav + Pblow) > 1) stop("Move probabilities larger one")
  
  consoleUpdates <- settings$consoleUpdates
  
  
  FUN <- setup$posterior$density
  
  if(!restart){
  # Initialize x and x2
    
  if(is.null(settings$startValue)){
      settings$startValue = setup$prior$sampler(2)
    }
    if(is.function(settings$startValue)){
      settings$startValue = settings$startValue(2)
    }
  x <- settings$startValue[1,]
  x2 <- settings$startValue[2,]
  
  # Evaluate 
  Eval <- FUN(x, returnAll = T) 
  Eval2 <- FUN(x2, returnAll = T)
  }else{
    x <- bayesianSetup$chain[[1]][nrow(bayesianSetup$chain[[1]]), 1:Npar]
    x2 <- bayesianSetup$chain[[2]][nrow(bayesianSetup$chain[[2]]), 1:Npar]
    
    Eval <- bayesianSetup$chain[[1]][nrow(bayesianSetup$chain[[1]]), (Npar+1):(Npar+3)]
    Eval2 <- bayesianSetup$chain[[2]][nrow(bayesianSetup$chain[[2]]), (Npar+1):(Npar+3)]
    
  }
  
  # Initialize chains
  chain <- matrix(NA, nrow = floor((iterations+1-burnin)/thin), ncol = Npar+3)
  chain2 <-  matrix(NA, nrow = floor((iterations+1-burnin)/thin), ncol = Npar+3)
  
  # Fill first values in chain
  chain[1,] <- c(x,Eval)
  chain2[1,] <- c(x2,Eval2)
  
  # Initialize counter for acceptance rate
  acceptance <- 0

  # Initialize counter
  counter <- 0
  
  
  for (i in 1:iterations) {
    
    move <- TwalkMove(Npar = Npar, FUN = FUN, x = x, 
                 Eval = Eval, x2 = x2, Eval2 = Eval2, at = at, aw = aw, pn1 = pn1, Ptrav = Ptrav, 
                 Pwalk = Pwalk, Pblow = Pblow, Phop = Phop)
    if(!is.na(move$alpha)){
    if (runif(1) < move$alpha) {
      x <- move$y
      Eval<- move$val
      x2 <- move$y2
      Eval2 <- move$val2
    }
    }
    
    if((i > burnin) && (i %% thin == 0) ){ # retain sample
    counter <- counter + 1
    chain[counter,] <- c(x, Eval)
    chain2[counter,] <- c(x2, Eval2)
    }
    
    if(settings$message){
    if( (i %% consoleUpdates == 0) | (i == iterations)) {
      cat("\r","Running Twalk-MCMC, chain ", settings$currentChain ,  "iteration" ,(i*2),"of",(iterations*2),
          ". Current logp ", Eval[1], Eval2[1]  ,". Please wait!","\r")
    flush.console()
    }
  }
  }
  colnames(chain) <- c(setup$names,"LP", "LL", "LPr")
  colnames(chain2) <- c(setup$names,"LP", "LL", "LPr")
  
  if(restart){ # Combine chains
    chain <- rbind(bayesianSetup$chain[[1]], chain)
    chain2 <- rbind(bayesianSetup$chain[[2]], chain2)
  }
  
  # Make sure chains have the right size
  chain <- chain[1:counter,] 
  chain2 <- chain2[1:counter,] 
  
  chain <- coda::mcmc.list(coda::mcmc(chain), coda::mcmc(chain2))
  
  out <- list(chain = chain, settings = settings)
  class(out) <- c("mcmcSampler", "bayesianOutput")
  return(out)
}

