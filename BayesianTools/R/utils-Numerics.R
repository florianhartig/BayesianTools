# inspired from https://gist.github.com/doobwa/941125
# maybe replace with logSumExp {matrixStats} which might be faster?

#' Funktion to compute log(sum(exp(x))
#' @author Florian Hartig
#' @param x values
#' @param mean logical, determines whether the mean should be used
#' @keywords internal
logSumExp<- function(x, mean = F) {
  # Computes log(sum(exp(x))
  # Uses offset trick to avoid numeric overflow: http://jblevins.org/notes/log-sum-exp
  if ( max(abs(x)) > max(x) ) offset <- min(x) else offset <- max(x)
  if (mean == T) out = log(mean(exp(x - offset))) + offset
  else out = log(sum(exp(x - offset))) + offset
  return(out)
}


#x = rep(-100, 100) 

#log(sum(exp(x)))
#logSumExp(x, mean = F)
#log(mean(exp(x)))
#logSumExp(x, mean = T)
