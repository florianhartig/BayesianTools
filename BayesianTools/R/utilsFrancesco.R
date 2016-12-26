as_mcmc.matrix = function(x, names, par = NULL) {
  dimnames(x)=list(NULL, names)
  if (is.null(par)) {
    attr(x, 'mcpar') = c(1, nrow(x), 1)
  } else {
    attr(x, 'mcpar') = par
  }
  class(x) = 'mcmc'
  return(x)
}

calc_mpsrf = function(x, end) {
  x = window(x, start(x), end, 1)
  Niter = niter(x)
  Nchain = nchain(x)
  Nvar = nvar(x)
  x = lapply(x, as.matrix)
  S2 = array(sapply(x, var, simplify=TRUE), dim=c(Nvar,Nvar,Nchain))
  W = apply(S2, c(1,2), mean)
  xbar = matrix(sapply(x, apply, 2, mean, simplify=TRUE), nrow=Nvar, ncol=Nchain)
  B = Niter * var(t(xbar))
  CW = chol(W)
  emax = eigen(backsolve(CW, t(backsolve(CW, B, transpose=TRUE)), transpose=TRUE), symmetric=TRUE, only.values=TRUE)$values[1]
  res = sqrt( (1 - 1/Niter) + (1 + 1/Nvar) * emax/Niter )
  return(res)
}

mpsrf = function(x, step = 50, ...) {
  if (nchain(x) < 2 | nvar(x) == 1) 
    stop("You need at least two chains and two parameters.")
  z = seq(start(x)-1+step, end(x), by = step)
  res = sapply(z, function(i) calc_mpsrf(x,i))
  x = data.frame(z,MPSRF = res)
  return(x)
}