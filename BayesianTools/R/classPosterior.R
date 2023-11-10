#' Creates a standardized posterior class
#' @author Florian Hartig
#' @param prior prior class
#' @param likelihood log likelihood density
#' @details Function is internally used in \code{\link{createBayesianSetup}} to create a standardized posterior class.
#' @export
createPosterior <- function(prior, likelihood){

  posterior <- function(x, returnAll = F){
    
    if (is.vector(x)){
      priorResult = prior$density(x) # Checking if outside prior to save calculation time
      if (! (priorResult == -Inf)) ll = likelihood$density(x)
      else ll = -Inf
      if (returnAll == F) return(ll + priorResult)    
      else return(c(ll + priorResult, ll, priorResult)) 
    
    } else if(is.matrix(x)){
      
      priorResult = prior$density(x) # Checking first if outside the prior to save calculation time
      feasible <- (! priorResult == -Inf)
      if (dim(x)[2] == 1) llResult <- likelihood$density(matrix(x[feasible, ], ncol = 1))
      else{
        if(TRUE %in% feasible) llResult <- likelihood$density(x[feasible, ])
        else llResult <- -Inf 
      }
      post = priorResult
      ll = priorResult
      ll[!feasible] = NA
      ll[feasible] = llResult
      post[feasible] = post[feasible] + llResult
      post[!feasible] = -Inf
      if (returnAll == F) return(post)    
      else{
        out <- cbind(post, ll, priorResult)
        colnames(out) = c("posterior", "likelihood", "prior")
        return(out)    
      } 
    }
    else stop("parameter must be vector or matrix")
  }

  out<- list(density = posterior)
  class(out) <- "posterior"
  return(out)
}


# likelihood <- function(x)stop("a")
# prior <- createPrior(function(x) sum(dunif(x, log = T)))
# 
# x = createPosterior(prior, likelihood)
# 
# x$density(c(0.2,0.2))
# prior$density(c(2,2))
# 
# 
# x = c(0.2,0.2)
