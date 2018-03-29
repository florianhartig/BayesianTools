p1 = createPrior(function(x) sum(dnorm(x, 0, 50, TRUE)), function(n=1) rnorm(n, 0, 50))
ll1 = function(x) sum(dnorm(x, 3, 2, TRUE))
bs1 = createBayesianSetup(ll1, p1)

p2 = createPrior(function(x) sum(mvtnorm::dmvnorm(x, c(0,0,0), diag(3) * c(100, 100, 100), TRUE)),
                function(n=1) mvtnorm::rmvnorm(n, c(0,0,0), diag(3) * c(100, 100, 100)))
ll2 = function(x) sum(mvtnorm::dmvnorm(x, c(4, 2, 7), diag(3) * c(5, 2, 15), TRUE))
bs2 = createBayesianSetup(ll2, p2)


# Metropolis needs one startValue
# and a covariance matrix
# negative inverse of the hessian is (asymptotic) covariance matrix
#  Matrix::nearPD(-MASS::ginv(numDeriv::hessian(target, preOpt$par)))$mat


# DE needs npar * 3 start values
# npar = bs$numPars
#  x = bs$prior$sampler(npar * 3)

#  startVals = apply(x, 1, optim, fn = target, hessian = F, control=list("fnscale" = -1, alpha = 2, beta = 0.7, gamma = 3.0), method="Nelder-Mead")
#  t(sapply(startVals, function(optResult) optResult$par, simplify = "matrix"))

# DEzs needs npar * 3 start values
# and a npar * 10 z matrix
#   draw from a mvnormal distribution around the start values?
#   randomize around startvalues?
#   do multiple optimizations?

# DREAM needs max(4, npar * 2) start values

# DREAMzs needs npar * 3 start values
# and a npar * 10 z matrix


preOptimize = function (bayesianSetup, sampler, ...) {
  target <- function(x) {
    out <- bayesianSetup$posterior$density(x)
    if (out == -Inf) out =  -1e20 # rnorm(1, mean = -1e20, sd = 1e-20)
    return(out)
  }

  nPar = bayesianSetup$numPars

  lower = if (!is.null(bayesianSetup$prior$lower)) bayesianSetup$prior$lower else rep(-1e20, nPar)
  upper = if (!is.null(bayesianSetup$prior$upper)) bayesianSetup$prior$upper else rep(+1e20, nPar)
  # startVal = if (!is.null(settings$startValue)) settings$startValue else (lower + upper) / 2
  
  if (sampler == "Metropolis") preOptimizeMetropolis(bayesianSetup, target, lower, upper, ...)
  else if (sampler == "DE") preOptimizeDE(bayesianSetup, target, lower, upper, ...)
  else stop(paste("Sampler", sampler, "not available."))
}

preOptimizeMetropolis = function(bs, target, lower, upper) {
  nPar = bs$numPars

  startVal = bs$prior$sampler()

  if(nPar > 1) optResul <- optim(par=startVal, fn = target, method="Nelder-Mead", hessian=F, control=list("fnscale"=-1, "maxit" = 10000))
  else optResul <- optim(par=startVal, fn = target, method="Brent", hessian=F, control=list("fnscale"=-1, "maxit" = 10000), lower = lower, upper = upper)
  
  optPar = optResul$par
  covMat = Matrix::nearPD(-MASS::ginv(numDeriv::hessian(target, optPar)))$mat
  list(startValue = optPar, cov = covMat)
}

preOptimizeDE = function(bs, target, lower, upper) {
  nPar = bs$numPars
  
  startVal = bs$prior$sampler(nPar * 3)
  
  if(nPar > 1) f <- function(x) optim(par=x, fn = target, method="Nelder-Mead", hessian=F, control=list("fnscale"=-1, "maxit" = 10000))
  else f <- function(x) optim(par=x, fn = target, method="Brent", hessian=F, control=list("fnscale"=-1, "maxit" = 10000), lower = lower, upper = upper)
  optResul = apply(startVal, 1, f)
  
  optPar = t(sapply(optResul, function(x) x$par))
  # covMat = Matrix::nearPD(-MASS::ginv(numDeriv::hessian(target, preOpt$par)))$mat
  # list(par = optPar, cov = covMat)
  list(startValue = optPar)
}

preOptimize(bs1, "Metropolis")
preOptimize(bs1, "DE")

preOptimize(bs2, "Metropolis")
preOptimize(bs2, "DE")



