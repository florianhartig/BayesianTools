

#' Differential-Evolution MCMC
#' @author Francesco Minunno and Stefan Paul
#' @param bayesianSetup a BayesianSetup with the posterior density function to be sampled from
#' @param settings list with parameter settings
#' @param startValue (optional) eiter a matrix with start population, a number to define the number of chains that are run or a function that samples a starting population.
#' @param iterations number of function evaluations.
#' @param burnin number of iterations treated as burn-in. These iterations are not recorded in the chain.
#' @param thin thinning parameter. Determines the interval in which values are recorded.
#' @param f scaling factor gamma
#' @param eps small number to avoid singularity
#' @param blockUpdate list determining whether parameters should be updated in blocks. For possible settings see Details.
#' @param message logical determines whether the sampler's progress should be printed
#' @references Braak, Cajo JF Ter. "A Markov Chain Monte Carlo version of the genetic algorithm Differential Evolution: easy Bayesian computing for real parameter spaces." Statistics and Computing 16.3 (2006): 239-249.
#' @export
#' @seealso \code{\link{DEzs}}
#' @details For blockUpdate the first element in the list determines the type of blocking.
#' Possible choices are
#' \itemize{
#'  \item{"none"}{ (default), no blocking of parameters}
#'  \item{"correlation"} { blocking based on correlation of parameters. Using h or k (see below)}
#'  \item{"random"} { random blocking. Using k (see below)}
#'  \item{"user"} { user defined groups. Using groups (see below)}
#'  }
#'  Further seven parameters can be specified. "k" determnined the number of groups, "h" the strength
#'  of the correlation used to group parameter and "groups" is used for user defined groups.
#'  "groups" is a vector containing the group number for each parameter. E.g. for three parameters 
#'  with the first two in one group, "groups" would be c(1,1,2).
#'  Further pSel and pGroup can be used to influence the choice of groups. In the sampling process
#'  a number of groups is randomly drawn and updated. pSel is a vector containing relative probabilities
#'  for an update of the respective number of groups. E.g. for always updating only one group pSel = 1.
#'  For updating one or two groups with the same probability pSel = c(1,1). By default all numbers
#'  have the same probability.
#'  The same principle is used in pGroup. Here the user can influence the probability of each group
#'  to be updated. By default all groups have the same probability.
#'  Finally "groupStart" defines the starting point of the groupUpdate and "groupIntervall" the intervall
#'  in which the groups are evaluated.

DE <- function(bayesianSetup, 
                  settings = list(
                    startValue = NULL, 
                    iterations = 10000, 
                    f = -2.38, 
                    burnin = 0, 
                    thin = 1,
                    eps = 0, 
                    consoleUpdates = 100, 
                    blockUpdate = list("none", k = NULL, h = NULL, pSel = NULL, pGroup = NULL, 
                    groupStart = 1000, groupIntervall = 1000), 
                    currentChain = 1,
                    message = TRUE
                  )
               ){

  if("bayesianOutput" %in% class(bayesianSetup)){
    restart <- TRUE
  } else restart <- FALSE
  
  
  if(restart){
    if(is.null(settings)) settings <- bayesianSetup$settings
    else settings <- applySettingsDefault(settings = settings, sampler = "DE")
    
  }else{
    # If nothing provided use default settings
    settings <- applySettingsDefault(settings = settings, sampler = "DE")
  }
  
  if(!restart){ 
    setup <- bayesianSetup
  }else{
    setup <- bayesianSetup$setup
  } 
  
  
  setup <- checkBayesianSetup(setup, parallel = settings$parallel) # calling parallel will check if requested parallelization in settings is provided by the BayesianSetup
  if(is.null(settings$parallel)) settings$parallel = setup$parallel # checking back - if no parallelization is provided, we use the parallelization in the BayesianSetup. We could also set parallel = F, but I feel it makes more sense to use the Bayesiansetup as default
  
  if(!restart){
    if(is.null(settings$startValue)){
      parLen = length(bayesianSetup$prior$sampler(1))
      X = bayesianSetup$prior$sampler(3 * parLen)
    }
    if(is.function(settings$startValue)){
      X = settings$startValue()
    }
    if(class(settings$startValue)[1] == "numeric"){
        X = bayesianSetup$prior$sampler(settings$startValue)
    }
    if(is.matrix(settings$startValue)) X <- settings$startValue
  }else{
    X <- bayesianSetup$X
  }
  
 # X = startValue
 if (!is.matrix(X)) stop("wrong starting values")

  FUN = setup$posterior$density
  
  ## Initialize blockUpdate parameters and settings
  blockdefault <- list("none", k = NULL, h = NULL, pSel = NULL, pGroup = NULL, 
                       groupStart = 1000, groupIntervall = 1000)
  
  if(!is.null(settings$blockUpdate)){
    blockUpdate <- modifyList(blockdefault, settings$blockUpdate)
    blockUpdate[[1]] <- settings$blockUpdate[[1]] # to catch first argument
    if(blockUpdate[[1]] == "none"){
      blockUpdateType <- "none"
      blocks = FALSE
      BlockStart = FALSE
    }else{
      groupStart <- blockUpdate$groupStart
      groupIntervall <- blockUpdate$groupIntervall
      blockUpdateType = blockUpdate[[1]] 
      blocks = TRUE
      ## Initialize BlockStart
      BlockStart = FALSE
      Bcount = 0
    }
  }else{
    blockUpdateType <- "none"
    blocks = FALSE
    BlockStart = FALSE
  }
  
  
 
  Npar <- ncol(X)
  Npop <- nrow(X)
  burnin <- settings$burnin/Npop
  n.iter <- ceiling(settings$iterations/Npop)
  
  if (n.iter < 2) stop ("The total number of iterations must be greater than the number of parameters to fit times 3.")
  
  lChain <- ceiling((n.iter - burnin)/settings$thin)+1
  #pChain <- array(NA, dim=c(n.iter*Npop, Npar+3))
  
  pChain <- array(NA, dim=c(lChain, Npar+3, Npop))
  
  
  colnames(pChain) <- c(setup$names, "LP", "LL", "LPr")
  
  counter <- 1
  iseq <- 1:Npop

  
  F2 = abs(settings$f)/sqrt(2*Npar)
  if (settings$f>0) F1 = F2 else F1 = 0.98
  
  logfitness_X <- FUN(X, returnAll = T)
  
  # Write first values in chain
  pChain[1,,] <- t(cbind(X,logfitness_X))
  
  # Print adjusted iterations
 # cat("Iterations adjusted to", n.iter*Npop,"to fit settings", "\n")

  ####
  eps <- settings$eps
  currentChain <- settings$currentChain
  iterations <- settings$iterations
  
  for (iter in 2:n.iter) {
    
    if (iter%%10) F_cur = F2 else F_cur = F1

    
    if(blocks){
      ### Update the groups. 
      if(iter == groupStart+ Bcount*groupIntervall){
        blockSettings <- updateGroups(chain = pChain[1:counter,, ], blockUpdate)
        BlockStart <- TRUE
        Bcount <- Bcount + 1
      }
    }
   ####
     
    for (i in iseq){
      # select to random different individuals (and different from i) in rr, a 2-vector
      
      rr <- sample(iseq[-i], 2, replace = FALSE)
      x_prop <- X[i,] + F_cur * (X[rr[1],]-X[rr[2],]) + eps * rnorm(Npar,0,1)
      
      if(BlockStart){
        # Get the current group and update the proposal accordingly
        Member <- getBlock(blockSettings)
        x_prop[-Member] <- X[i,-Member]
        ####
      }
      
      logfitness_x_prop <- FUN(x_prop, returnAll = T)
      if(!is.na(logfitness_x_prop[1] - logfitness_X[i,1])){ # To catch possible error
        if ((logfitness_x_prop[1] - logfitness_X[i,1] ) > log(runif(1))){
          X[i,] <- x_prop
          logfitness_X[i,] <- logfitness_x_prop
        }
      }
    } #iseq
    if ((iter > burnin) && (iter %% settings$thin == 0) ) { # retain sample
      counter <- counter+1
      pChain[counter,,] <- t(cbind(X,logfitness_X))
      
    }
    
    if(settings$message){  
      if( (iter %% settings$consoleUpdates == 0) | (iter == n.iter)) cat("\r","Running DE-MCMC, chain ", currentChain, 
                                                                         "iteration" ,iter*Npop,"of",n.iter*Npop,". Current logp ",
                                                                         logfitness_X[,1],
                                                                         "Please wait!","\r")
      flush.console()
    }

  } # n.iter
  iterationsOld <- 0
  
  pChain <- pChain[1:counter,,]
  
  if(restart){ # Combine chains
    newchains <- array(NA, dim = c((counter+nrow(bayesianSetup$chain[[1]])), (Npar+3), Npop))
    
    for(i in 1:Npop){
      for(k in 1:(Npar+3)){
        newchains[,k,i] <- c(bayesianSetup$chain[[i]][,k],pChain[,k,i])
      }
    }
    pChain <- newchains
  }
  


  pChain<- coda::as.mcmc.list(lapply(1:Npop,function(i) coda::as.mcmc(pChain[,1:(Npar+3),i])))
  
  
  list(Draws = pChain, X = as.matrix(X[,1:Npar]))
  }
