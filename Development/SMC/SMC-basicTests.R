
ll = function(x) sum(dnorm(x, log = T))
setup = createBayesianSetup(ll, lower = c(-10, -10), upper = c(10, 10))

settings = list(iterations = 40000)
out = runMCMC(bayesianSetup = setup, settings = settings)
ks.test(getSample(out, start = 5000)[,1], msm::ptnorm, 0, 1 , -10, 10)


out = smcSampler(bayesianSetup = setup, initialParticles = 2000, iterations = 1, ess.limit = 10000000)
hist(out$particles)
qqnorm(out$particles)
ks.test(out$particles[,1], msm::ptnorm, 0, 1 , -10, 10)

