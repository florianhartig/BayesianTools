p = createPrior(function(x) sum(dnorm(x, 0, 50, TRUE)), function(n=1) rnorm(n, 0, 50))
ll = function(x) sum(dnorm(x, 3, 2, TRUE))
bs = createBayesianSetup(ll, p)

p = createPrior(function(x) sum(mvtnorm::dmvnorm(x, c(0,0,0), diag(3) * c(100, 100, 100), TRUE)),
                function(n=1) mvtnorm::rmvnorm(n, c(0,0,0), diag(3) * c(100, 100, 100)))
ll = function(x) sum(mvtnorm::dmvnorm(x, c(4, 2, 7), diag(3) * c(5, 2, 15), TRUE))
bs = createBayesianSetup(ll, p)

preOptimize = function(bayesianSetup) {
  target <- function(x) {
    out <- bayesianSetup$posterior$density(x)
    if (out == -Inf) out =  -1e20 # rnorm(1, mean = -1e20, sd = 1e-20)
    return(out)
  }
  
  settings = bayesianSetup$settings
  nPar = bayesianSetup$numPars
  
  lower = if (!is.null(bayesianSetup$prior$lower)) bayesianSetup$prior$lower else rep(-1e20, nPar)
  upper = if (!is.null(bayesianSetup$prior$upper)) bayesianSetup$prior$upper else rep(+1e20, nPar)
  startValue = if (!is.null(settings$startValue)) settings$startValue else (lower + upper) / 2

  if(nPar > 1) optresul <- optim(par=startValue, fn = target, method="Nelder-Mead", hessian=F, control=list("fnscale"=-1, "maxit" = 10000))
  else optresul <- optim(par=startValue, fn = target, method="Brent", hessian=F, control=list("fnscale"=-1, "maxit" = 10000), lower = lower, upper = upper)      
  optresul
}

preOpt = preOptimize(bs)

# Metropolis needs one startValue
# and a covariance matrix
# negative inverse of the hessian is (asymptotic) covariance matrix
Matrix::nearPD(-MASS::ginv(numDeriv::hessian(target, preOpt$par)))$mat


# DE needs npar * 3 start values
npar = bs$numPars
x = bs$prior$sampler(npar * 3)

startVals = apply(x, 1, optim, fn = target, hessian = F, control=list("fnscale" = -1, alpha = 2, beta = 0.7, gamma = 3.0), method="Nelder-Mead")
t(sapply(startVals, function(optResult) optResult$par, simplify = "matrix"))

# DEzs needs npar * 3 start values
# and a npar * 10 z matrix







# if(bayesianSetup$numPars > 1) optresul <- optim(par=settings$startValue,fn=target, method="Nelder-Mead", hessian=F, control=list("fnscale"=-1, "maxit" = 10000))
# else optresul <- optim(par=settings$startValue,fn=target, method="Brent", hessian=F, control=list("fnscale"=-1, "maxit" = 10000), lower = bayesianSetup$prior$lower, upper = bayesianSetup$prior$upper)
# 
# 
# 
# if(is.null(settings$message) || settings$message == TRUE){
#   cat("BT runMCMC: trying to find optimal start and covariance values", "\b")
# }
# 
# target <- function(x){
#   out <- bayesianSetup$posterior$density(x)
#   if (out == -Inf) out =  -1e20 # rnorm(1, mean = -1e20, sd = 1e-20)
#   return(out)
# }
# 
# try( {
#   if(bayesianSetup$numPars > 1) optresul <- optim(par=settings$startValue,fn=target, method="Nelder-Mead", hessian=F, control=list("fnscale"=-1, "maxit" = 10000))
#   else optresul <- optim(par=settings$startValue,fn=target, method="Brent", hessian=F, control=list("fnscale"=-1, "maxit" = 10000), lower = bayesianSetup$prior$lower, upper = bayesianSetup$prior$upper)
#   settings$startValue = optresul$par
#   hessian = numDeriv::hessian(target, optresul$par)
# 
# 
#   proposalGenerator$covariance = as.matrix(Matrix::nearPD(MASS::ginv(-hessian))$mat)
#   #proposalGenerator$covariance = MASS::ginv(-optresul$hessian)
# 
#   # Create objects for startValues and covariance to add space between values
#   startV <-covV <- character()
# 
#   for(i in 1:length(settings$startValue)){
#     startV[i] <- paste(settings$startValue[i], "")
#   }
#   for(i in 1:length( proposalGenerator$covariance)){
#     covV[i] <- paste( proposalGenerator$covariance[i], "")
#   }
# 
#   if(is.null(settings$message) || settings$message == TRUE){
#     message("BT runMCMC: Optimization finished, setting startValues to " ,
#             startV, " - Setting covariance to " , covV)
#   }
# 
#   proposalGenerator = updateProposalGenerator(proposalGenerator)
# 
# }
# , silent = FALSE)
# 
