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

bayesianSetup <- createBayesianSetup(likelihood, prior, names = rownames(refPars)[parSel])


##### Time MCMC

out.dir <- paste0("/home/fr/fr_fr/fr_ms1719/out/VSEMb", number, "/")
out.file <- paste0(out.dir, "mcmcRef_np.txt")
set.seed(number)

cat(c("time", "distance"), file = out.file, append = FALSE)

time = proc.time()

runtime <- distance <- d1 <- d2 <- vector("numeric", 2000)

settings <- list(iterations = 1000, nrChains = 3, message = TRUE)

time.count <- system.time(out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings))[3]

for(i in 2:2000){
  print(c("i", i))
  time.cur <- system.time(out <- runMCMC(out))[3] 
  time.count <- time.count + time.cur
  runtime[i] <- time.count
  distance[i] <- getSampleDistance(ref.sample,getSample(out, start = floor((i * 1000)/3/2) ),type = "D" )
  cat("\n", file = out.file, append = TRUE)
  cat(c(runtime[i], distance [i]), file = out.file, append = TRUE)
}
