# inspired from https://gist.github.com/doobwa/941125
# maybe replace with logSumExp {matrixStats} which might be faster?

#' Funktion to compute log(sum(exp(x))
#' @author Florian Hartig
#' @param x values
#' @param mean logical, determines whether the mean should be used instead of the sum
#' @details This function computes log(sum(exp(x)), using the offset trick to avoid numeric overflow, see, e.g. http://jblevins.org/notes/log-sum-exp. The mean option allows calculating logMeanExp
#' 
#' @keywords internal
logSumExp<- function(x, mean = F) {
  # 
  
  nObs = length(x)   
  
  if(any(x == Inf)) stop("BayesianTools::logSumExp: positive infinity values in log probabilities")
  if(any(x == -Inf )){
    message("BayesianTools::logSumExp: encountered -Inf in logSumExp - value was removed")    
    x = x[x != -Inf] 
  } 
  
  # seems that this created problems in the presence of small values,
  # doesn't seem to be a need to shift towards min
  # if ( max(abs(x)) > max(x) ) offset <- min(x) else offset <- max(x)
  offset <- max(x)
  if (mean == T) out = log(sum(exp(x - offset))/nObs) + offset
  else out = log(sum(exp(x - offset))) + offset
  return(out)
}

# Unit test in test-utils-Numerics