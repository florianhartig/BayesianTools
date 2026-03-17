
ll = function(x){
  #for(i in 1:10000) mean(rnorm(10))
  sum(dnorm(x, log = T))
} 
setup = createBayesianSetup(ll, lower = c(-10, -10), upper = c(10, 10))

settings = list(iterations = 40000)
out = runMCMC(bayesianSetup = setup, settings = settings)
ks.test(getSample(out, start = 5000)[,1], msm::ptnorm, 0, 1 , -10, 10)
plot(out)


initialPart = matrix(rnorm(10000), ncol = 2)


out = smcSampler(bayesianSetup = setup, initialParticles = initialPart, iterations = 5, adaptive= F, ess.limit = 20000, mutate = "D")

par(mfrow = c(1,2))
hist(initialPart, xlim  = c(-4,4))
hist(out$particles, xlim  = c(-4,4))

sd(initialPart)
sd(out$particles)
1/sqrt(2)



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