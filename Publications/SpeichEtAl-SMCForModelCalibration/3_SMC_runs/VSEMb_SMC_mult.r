arg <- commandArgs(TRUE)
number <- eval(parse(text = arg))

library(BayesianTools)

#################
# Load saved (synthetic) data and reference MCMC output
load("./out/reference_VSEMb_list.RData")

# Use only sample from reference; delete reference to free memory
ref.sample <- getSample(save.list[[3]], start = 3333)

PAR <- save.list[[1]]
obs <- save.list[[2]]
rm(save.list)

################## 
# Create Bayesian setup

# load reference parameter definition (upper, lower prior)
refPars <- VSEMgetDefaults()
# this adds one additional parameter for the likelihood standard deviation (see below)
refPars[12,] <- c(2, 0.1, 4) 
rownames(refPars)[12] <- "error-sd"

#parSel = c(1:6, 12)
parSel <- c(1:3, 5:7, 12)

# here is the likelihood 
likelihood <- function(par, sum = TRUE){
  # set parameters that are not calibrated on default values 
  x = refPars$best
  x[parSel] = par
  predicted <- VSEM(x[1:11], PAR) # replace here VSEM with your model 
  predicted[,1] = 1000 * predicted[,1] # this is just rescaling
  diff <- c(predicted[,1:4] - obs[,1:4]) # difference betweeno observed and predicted
  # univariate normal likelihood. Note that there is a parameter involved here that is fit
  llValues <- dnorm(diff, sd = x[12], log = TRUE)  
  if (sum == FALSE) return(llValues)
  else return(sum(llValues))
}

# optional, you can also directly provide lower, upper in the createBayesianSetup, see help
prior <- createUniformPrior(lower = refPars$lower[parSel], 
                            upper = refPars$upper[parSel], best = refPars$best[parSel])

bayesianSetup <- createBayesianSetup(likelihood, prior, names = rownames(refPars)[parSel], parallel=TRUE)

####################
# Settings for different SMC runs
out.dir <- paste0("/home/fr/fr_fr/fr_ms1719/out/VSEMb", number, "/")
out.file <- paste0(out.dir, "smc_out.txt")
set.seed(number)

settings <- expand.grid(particles = c(5000, 20000, 50000, 100000), ess.limit = c(0.5, 0.75, 0.9), proposalScale = c(0.01, 0.1, 0.333, 0.5), resamplingSteps = c(2,5,10,20, 30))

cat(c("particles", "ess.limit", "proposalScale", "mcmcSteps", "d", "time"), file=out.file, append=FALSE)

for(i in 1:nrow(settings)){
  cur.time <- system.time(smc.VSEMb <- smcSampler(bayesianSetup = bayesianSetup, initialParticles = settings$particles[i], ess.limit = (settings$particles[i] * settings$ess.limit[i]), proposalScale=settings$proposalScale[i], mutate = "DE", resamplingSteps = settings$resamplingSteps[i]))
  dist.d <- getSampleDistance(ref.sample, getSample(smc.VSEMb), type = "D")
  cat("\n", file=out.file, append=TRUE)
  cat(c(settings$particles[i], settings$ess.limit[i], settings$proposalScale[i], settings$resamplingSteps[i], dist.d, cur.time[3]), file=out.file, append=TRUE)
  
  save.file <- paste0(out.dir, "p_", settings$particles[i], "l_", settings$ess.limit[i]*100, "s_", settings$proposalScale[i], "r_", settings$resamplingSteps[i], ".RData")
  save(smc.VSEMb, file = save.file)
}