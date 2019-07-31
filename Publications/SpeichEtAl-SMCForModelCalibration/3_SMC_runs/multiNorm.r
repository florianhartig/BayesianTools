library(BayesianTools)

## Generate a multivariate normal density with strong correlation
#norm.dens <- generateTestDensityMultiNormal(sigma="strongcorrelation")

## Test multivariate normal density with strong correlation
norm.setup <- createBayesianSetup(likelihood = testDensityMultiNormal, lower = c(-5, -5, -5), upper = c(5,5,5))

ref.sampler <- runMCMC(bayesianSetup = norm.setup, sampler = "DEzs", settings = list(iterations = 100000, nrChains = 3))
ref.sample <- getSample(ref.sampler, start = floor(100000/3/2))

# Check convergence
#gelmanDiagnostics(sampler = ref.sampler, start = floor(100000/3/2) )


out.dir <- "/home/fr/fr_fr/fr_ms1719/out/multiNorm/"


# Check how long it takes for MCMC to be as good as the reference

# for(iter in 1:5){
#   out.file <- paste0(out.dir, "mcmcRef", iter, ".txt")
#   cat(c("time", "distance", "mean1", "mean2", "mean3", "sd1", "sd2", "sd3", "converged"), file = out.file, append = FALSE)
#   
#   #time <- distance <- mean1 <- mean2 <- mean3 <- sd1 <- sd2 <- sd3 <- rep(NA, 100)
#   
#   time.count <- system.time(out <- runMCMC(bayesianSetup = norm.setup, sampler = "DEzs", settings = list(iterations = 1000, nrChains = 3)))[3]
#   time <- time.count
#   distance <- getSampleDistance(getSample(ref.sampler, start = floor(100000/3/2)), getSample(out, start = floor(1000/3/2)), type="D")
#   cur.sample <- getSample(out, start = floor((1000)/3/2))
#   mean1 <- mean(cur.sample[,1])
#   mean2 <- mean(cur.sample[,2])
#   mean3 <- mean(cur.sample[,3])
#   
#   sd1 <- sd(cur.sample[,1])
#   sd2 <- sd(cur.sample[,2])
#   sd3 <- sd(cur.sample[,3])
#   
#   for(i in 2:100){
#     cur.time <- system.time(out <- runMCMC(out))[3]
#     time.count <- time.count + cur.time
#     distance <- getSampleDistance(getSample(ref.sampler, start = floor(100000/3/2)), getSample(out, start = floor((i * 1000)/3/2)), type="D")
#     
#     cur.sample <- getSample(out, start = floor((i * 1000)/3/2))
#     mean1 <- mean(cur.sample[,1])
#     mean2 <- mean(cur.sample[,2])
#     mean3 <- mean(cur.sample[,3])
#     
#     sd1 <- sd(cur.sample[,1])
#     sd2 <- sd(cur.sample[,2])
#     sd3 <- sd(cur.sample[,3])
#     
#     gelman <- gelmanDiagnostics(out)
#     if(max(gelman$psrf) <= 1.05 & gelman$mpsrf <=1.2){
#       converged <- TRUE
#     } else {
#       converged <- FALSE
#     }
#     
#     cat("\n", file = out.file, append = TRUE)
#     cat(c(time.count, distance, mean1, mean2, mean3, sd1, sd2, sd3, converged), file = out.file, append = TRUE)
#     
#   }
#   
#   
# }

########
# SMC experiment

settings <- expand.grid(particles = c(50, 100, 1000), ess.limit = c(0.5, 0.75, 0.9), proposalScale = c(0.01, 0.1, 0.333, 0.5), resamplingSteps = c(2,5,10,20, 30))



#norm.setup.par <- createBayesianSetup(likelihood = testDensityMultiNormal, lower = c(-5, -5, -5), upper = c(5,5,5), parallel = TRUE)

out.file <- paste0(out.dir, "smc_out_multiNorm.txt")
cat(c("time", "sd.time", "distance", "sd.distance", "converged"), file = out.file, append = FALSE)


for(i in 1:nrow(settings)){
  smc.list <- list()
  smc.dist.vec <- time.vec <- vector("numeric", 5)
  
  for(j in 1:5){
    cur.time <- system.time(cur.smc <- smcSampler(bayesianSetup = norm.setup, initialParticles = settings$particles[i], ess.limit = (settings$particles[i] * settings$ess.limit[i]), proposalScale=settings$proposalScale[i], mutate = "DE", resamplingSteps = settings$resamplingSteps[i]))[3]
    time.vec[j] <- cur.time
    smc.dist.vec[j] <- getSampleDistance(ref.sample, getSample(cur.smc), type="D")
    smc.list[[j]] <- coda::mcmc(cur.smc$particles)
  }
  
  time <- mean(time.vec)
  sd.time <- sd(time.vec)
  distance <- mean(smc.dist.vec)
  sd.distance <- sd(smc.dist.vec)
  
  gelman <- coda::gelman.diag(smc.list)
  if(max(gelman$psrf) <= 1.05 & gelman$mpsrf <=1.2){
    converged <- TRUE
  } else {
    converged <- FALSE
  }
  
  cat("\n", file = out.file, append = TRUE)
  cat(c(time, sd.time, distance, sd.distance, converged), file = out.file, append = TRUE)
}
