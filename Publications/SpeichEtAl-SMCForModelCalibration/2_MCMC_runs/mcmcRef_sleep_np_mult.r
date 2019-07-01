arg <- commandArgs(TRUE)
number <- eval(parse(text = arg))

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


##### Prior and setup
prior <- createUniformPrior(lower = min.pars, upper = max.pars)
pgSetup <- createBayesianSetup(likelihood = likelihood, prior = prior, parallel = TRUE)

# After the script has been run once: load saved result
#load(file = "./out/reference_3pg.RData")




##### Time MCMC

out.dir <- paste0("/home/fr/fr_fr/fr_ms1719/out/threePGN_sleep", number, "/")
out.file <- paste0(out.dir, "mcmcRef_np.txt")
set.seed(number)

#cat(c("time", "d1", "d2", "distance"), file = "/home/fr/fr_fr/fr_ms1719/out/mcmcRef.txt", append = FALSE)
cat(c("time", "distance"), file = out.file, append = FALSE)

time = proc.time()

#mcmcRes = data.frame(time = rep(NA, 2000), distance = rep(NA, 2000), d1 = rep(NA, 2000), d2 = rep(NA, 2000))
runtime <- distance <- d1 <- d2 <- vector("numeric", 2000)

settings <- list(iterations = 1000, nrChains = 3, message = TRUE)

time.count <- system.time(out <- runMCMC(bayesianSetup = pgSetup, sampler = "DEzs", settings = settings))[3]

for(i in 2:2000){
  print(c("i", i))
  #out <- runMCMC(out)
  
  time.cur <- system.time(out <- runMCMC(out))[3] 
  #print(c("time.cur", time.cur))
  time.count <- time.count + time.cur
  runtime[i] <- time.count
  #runtime[i] <- 666
  #out.sample <- getSample(out)
  #dist.bh <- getSampleDistance(ref.sample,out.sample,type = "BHs" )
  #print(summary(out))
  #d1[i] <- dist.bh[1]
  #d2[i] <- dist.bh[2]
  distance[i] <- getSampleDistance(ref.sample,getSample(out, start = floor((i * 1000)/3/2) ),type = "D" )
  cat("\n", file = out.file, append = TRUE)
  #cat(c(runtime[i], d1[i], d2[i], distance [i]), file = "/home/fr/fr_fr/fr_ms1719/out/mcmcRef_np.txt", append = TRUE)
  cat(c(runtime[i], distance [i]), file = out.file, append = TRUE)
}

