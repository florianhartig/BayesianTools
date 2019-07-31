arg <- commandArgs(TRUE)
number <- eval(parse(text = arg))
#library(Rmpi)
#cat("mpi comm size is ",mpi.comm.size(),"\n")
#cat("mpi universe size is ",mpi.universe.size(),"\n")
#mpi.spawn.Rslaves(nslaves = 199)
#library(Rmpi)
#library(snow)
#cl <- parallel::makeCluster(199, type="MPI")
library(BayesianTools)
library(threePGN)

#################
# Load saved (synthetic) data and reference MCMC output
load("./out/reference_threePGN_list.RData")

# Use only sample from reference; delete reference to free memory
ref.sample <- getSample(save.list[[3]], start = 3333)

ref.ba <- save.list[[1]]
ref.gpp <- save.list[[2]]

rm(save.list)

################# 
# Generate reference data
set.seed(123)

#firstRun <- r3pgn(siteData = data_site[1:3,], climate = data_climate, parameters = data_param[,2], outputs = c(1:5, 10:12, 26:27))


#ref.ba <- firstRun$output[,,1][,3] + rnorm(nrow(firstRun$output[,,1]), sd=1)
#ref.gpp <- firstRun$output[,,1][,9] + rnorm(nrow(firstRun$output[,,1]), sd=1)

#################
# Create Bayesian setup

ref.pars <- c(data_param$mode, 1, 1)
min.pars <- c(data_param$min, 0.2, 0.2)
max.pars <- c(data_param$max, 2, 2)

##### Likelihood function
likelihood <- function(ref.pars){
  
  library(threePGN)
  pgrun <- r3pgn(siteData = data_site[1:3,], climate = data_climate, parameters = ref.pars[1:51], outputs = c(1:5,10:12, 26:27))
  sim.ba  <- pgrun$output[,,1][,3]
  sim.gpp <- pgrun$output[,,1][,9]
  
  like <- sum(dnorm(ref.ba-sim.ba, sd=ref.pars[52], log=TRUE)) + sum(dnorm(ref.gpp-sim.gpp, sd=ref.pars[53], log=TRUE))
  Sys.sleep(0.05)
  return(like)
}

##### Parallel cluster
# library(Rmpi)
# mpi.spawn.Rslaves()
#cl <- parallel::makeCluster(199, type="MPI")
#pLikelihood <- function(ref.pars){
# parallel::parApply(cl=cl, X = ref.pars, MARGIN = 1, FUN = likelihood)
#}
#parallel::clusterExport(cl, varlist = list(likelihood, min.pars, max.pars, ref.pars, ref.ba, ref.gpp))
#print(c("cores", parallel::detectCores(TRUE)))

##### Prior and setup
prior <- createUniformPrior(lower = min.pars, upper = max.pars)
pgSetup <- createBayesianSetup(likelihood = likelihood, prior = prior, parallel = TRUE)
#pgSetup <- createBayesianSetup(likelihood = pLikelihood, prior = prior, parallel = "external")
#pgSetup <- createBayesianSetup(likelihood = likelihood, prior = prior, parallel = 199)

#################
# Create a reference run with MCMC-DEzs

# On first execution: run MCMC and save result
#reference <- runMCMC(bayesianSetup = pgSetup, sampler = "DEzs", settings = list(nrChains=3, iterations = 2000000, thin = 10))
#save(reference, file = "./out/reference_3pg.RData")

# After the script has been run once: load saved result
#load(file = "./out/reference_3pg.RData")

# Use only sample from reference; delete reference to free memory
#ref.sample <- getSample(reference, start = 33333)
#rm(reference)


#################
# SMC runs
diagnostics <- function(particles, ref.sample){
  return(getSampleDistance(particles, ref.sample, type="BHs"))
}

out.dir <- paste0("/home/fr/fr_fr/fr_ms1719/out/threePGN_sleep30_", number, "/")
out.file <- paste0(out.dir, "smc_out.txt")
set.seed(number)

settings <- expand.grid(particles = c(5000, 20000, 50000, 100000), ess.limit = c(0.5, 0.75, 0.9), proposalScale = c(0.01, 0.1, 0.333, 0.5), resamplingSteps = c(2,5,10,20,30))

cat(c("particles", "ess.limit", "proposalScale", "mcmcSteps", "distance.1", "distance.2", "d", "time"), file=out.file, append=FALSE)


for(i in 1:nrow(settings)){
  cur.time <- system.time(smc.3pg <- smcSampler(bayesianSetup = pgSetup, initialParticles = settings$particles[i], ess.limit = (settings$particles[i] * settings$ess.limit[i]), proposalScale=settings$proposalScale[i], mutate = "DE", reference = ref.sample, diagnostics = diagnostics, resamplingSteps = settings$resamplingSteps[i]))
  dist.d <- getSampleDistance(ref.sample, getSample(smc.3pg), type = "D")
  cat("\n", file=out.file, append=TRUE)
  cat(c(settings$particles[i], settings$ess.limit[i], settings$proposalScale[i], settings$resamplingSteps[i], tail(smc.3pg$info$diagnostics, 1)[[1]][1], tail(smc.3pg$info$diagnostics, 1)[[1]][2], dist.d, cur.time[3]), file=out.file, append=TRUE)
  
  save.file <- paste0(out.dir, "p_", settings$particles[i], "l_", settings$ess.limit[i]*100, "s_", settings$proposalScale[i], "r_", settings$resamplingSteps[i], ".RData")
  save(smc.3pg, file = save.file)
}
 