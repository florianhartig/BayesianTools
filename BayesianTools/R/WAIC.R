
# TODO - implement WAIC as AIC, can look at https://github.com/jrnold/mcmcStats/blob/master/R/waic.R, check against http://finzi.psych.upenn.edu/library/blmeco/html/WAIC.html, https://cran.r-project.org/web/packages/loo/index.html, http://stats.stackexchange.com/questions/173128/watanabe-akaike-widely-applicable-information-criterion-waic-using-pymc


#' calculates the WAIC 
#' @author Florian Hartig
#' @param bayesianOutput an object of class BayesianOutput. Must implement a log-likelihood density function that can return point-wise log-likelihood values ("sum" argument).
#' @param numSamples the number of samples to calculate the WAIC
#' @param ... optional values to be passed on the the getSample function 
#' @note The function requires that the likelihood passed on to BayesianSetup contains the option sum = T/F, with defaul F. If set to true, the likelihood for each data point must be returned. 
#' @details
#'
#'
#' The WAIC is constructed as
#' \ifelse{html}{\out{<center><em>WAIC = -2 * (lppd - p<sub>WAIC</sub>)</em></center><br>}}{\deqn{WAIC = -2 * (lppd - p_{WAIC})}}
#'
#' The lppd (log pointwise predictive density), defined in Gelman et al., 2013, eq. 4 as
#' \ifelse{html}{\out{<center><em>lpdd = &Sigma;<sub>i=1</sub><sup>n</sup> log 1/S &Sigma;<sub>s=1</sub><sup>S</sup> p(y<sub>i</sub> | &Theta;<sup>s</sup>)</em></center><br>}}{\deqn{lppd = \sum_{i=1}^n \log \left(\frac{1}{S} \sum_{s=1}^S p(y_i | \theta^s)\right)}}
#' 
#'
#' The value of \ifelse{html}{\out{<em>p<sub>WAIC</sub></em>}}{\eqn{p_WAIC}} can be calculated in two ways, the method used is determined by the
#' \code{method} argument.
#'
#' Method 1 is defined as,
#' \ifelse{html}{\out{<center><em>p<sub>WAIC1</sub> = 2 &Sigma;<sub>i=1</sub><sup>n</sup> \{log[1/S &Sigma;<sub>s=1</sub><sup>S</sup> p(y<sub>i</sub> | &Theta;<sup>s</sup>)] - 1/S &Sigma;<sub>s=1</sub><sup>S</sup> log p(y<sub>i</sub> | &Theta;<sup>s</sup>)\}</em></center><br>}}{\deqn{p_{WAIC1} = 2 \sum_{i=1}^{n} (\log (\frac{1}{S} \sum_{s=1}^{S} p(y_i \ \theta^s)) - \frac{1}{S} \sum_{s = 1}^{S} \log p(y_i | \theta^s))}}
#' Method 2 is defined as,
#' \ifelse{html}{\out{<center><em>p<sub>WAIC2</sub> = 2 &Sigma;<sub>i=1</sub><sup>n</sup> V<sub>s=1</sub><sup>S</sup>[log p(y<sub>i</sub> | &Theta;<sup>s</sup>)]</em></center><br>}}{\deqn{p_{WAIC2} = 2 \sum_{i=1}^{n} V_{s=1}^{S} (\log p(y_i | \theta^s))}}
#' where  \ifelse{html}{\out{<em>V<sub>s=1</sub><sup>S</sup></em>}}{\eqn{V_{s=1}^{S}}} is the sample variance.
#' 
#' @references
#' 
#' Gelman, Andrew and Jessica Hwang and Aki Vehtari (2013), "Understanding Predictive Information Criteria for Bayesian Models," \url{http://www.stat.columbia.edu/~gelman/research/unpublished/waic_understand_final.pdf}.
#'
#' Watanabe, S. (2010). "Asymptotic Equivalence of Bayes Cross Validation and Widely Applicable Information Criterion in Singular Learning Theory", Journal of Machine Learning Research, \url{http://www.jmlr.org/papers/v11/watanabe10a.html}.
#' @example /inst/examples/WAICHelp.R
#' @seealso \code{\link{DIC}}, \code{\link{MAP}}, \code{\link{marginalLikelihood}}
#' @export
WAIC <- function(bayesianOutput, numSamples = 1000, ...){
  
 
  
  x = getSample(bayesianOutput, parametersOnly = F,  ...)
  
  # Catch nPars and ll density for mcmcList
  if("mcmcSamplerList" %in% class(bayesianOutput)){
    if (bayesianOutput[[1]]$setup$pwLikelihood == FALSE)
      stop("WAIC can only be applied if the likelihood density can be returned point-wise ('sum' argument, see examples).")
    nPars = bayesianOutput[[1]]$setup$numPars
    llDensity <- bayesianOutput[[1]]$setup$likelihood$density
  }else{
    if (bayesianOutput$setup$pwLikelihood == FALSE)
      stop("WAIC can only be applied if the likelihood density can be returned point-wise ('sum' argument, see examples).")
    nPars = bayesianOutput$setup$numPars
    llDensity <- bayesianOutput$setup$likelihood$density
  } 

  # x = getSample(bayesianOutput, parametersOnly = F)  
  
  i <- sample.int(nrow(x),numSamples,replace=TRUE) # should replace really be true?

  # Calculating log pointwise posterior predictive density, see Gelman et al., 2013, eq. 4
  # This is simply the mean likelihood of the posterior samples 
  
  pointWiseLikelihood = t(llDensity(x[i,1:nPars], sum = F))
  
  lppd <- sum(apply(pointWiseLikelihood, 2, BayesianTools::logSumExp, mean = T))

  # Calculating pWAIC options

  pWAIC1 <- 2 * sum(apply(pointWiseLikelihood, 2, function(y) logSumExp(y, mean = T) - mean(y)) )

  pWAIC2 <- sum(apply(pointWiseLikelihood, 2, var))
  out = list ( WAIC1 = -2 * (lppd - pWAIC1), WAIC2 = -2 * (lppd - pWAIC2),  lppd = lppd, pWAIC1 = pWAIC1, pWAIC2 = pWAIC2)
  return(out)
}