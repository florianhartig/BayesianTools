
ll = function(x){
  for(i in 1:10000) mean(rnorm(10))
  sum(dnorm(x, log = T))
} 
setup = createBayesianSetup(ll, lower = c(-10, -10), upper = c(10, 10))

settings = list(iterations = 40000)
out = runMCMC(bayesianSetup = setup, settings = settings)
ks.test(getSample(out, start = 5000)[,1], msm::ptnorm, 0, 1 , -10, 10)
plot(out)

out = smcSampler(bayesianSetup = setup, initialParticles = 2000, iterations = 2, ess.limit = 10000000, diagnostics = function(x)sum(mean(x)))
hist(out$particles)
qqnorm(out$particles)
ks.test(out$particles[,1], msm::ptnorm, 0, 1 , -10, 10)

out$info$diagnostics
out
summary(out)
plot(out)


ptm = proc.time() 

setup = createBayesianSetup(ll, lower = c(-10, -10), upper = c(10, 10))

out = smcSampler(bayesianSetup = setup, initialParticles = 2000, iterations = 2, ess.limit = 10000000, diagnostics = function(x)sum(mean(x)))

proc.time() - ptm


ptm = proc.time() 

setup = createBayesianSetup(ll, lower = c(-10, -10), upper = c(10, 10), parallel = T)

out = smcSampler(bayesianSetup = setup, initialParticles = 2000, iterations = 2, ess.limit = 10000000, diagnostics = function(x)sum(mean(x)))

proc.time() - ptm