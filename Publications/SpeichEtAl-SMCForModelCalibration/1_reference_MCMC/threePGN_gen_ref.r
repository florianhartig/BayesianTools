library(BayesianTools)
library(threePGN)

################# 
# Generate reference data
set.seed(123)

firstRun <- r3pgn(siteData = data_site[1:3,], climate = data_climate, parameters = data_param[,2], outputs = c(1:5, 10:12, 26:27))


ref.ba <- firstRun$output[,,1][,3] + rnorm(nrow(firstRun$output[,,1]), sd=1)
ref.gpp <- firstRun$output[,,1][,9] + rnorm(nrow(firstRun$output[,,1]), sd=1)

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
  return(like)
} 


##### Prior and setup
prior <- createUniformPrior(lower = min.pars, upper = max.pars)
pgSetup <- createBayesianSetup(likelihood = likelihood, prior = prior, parallel = FALSE)

reference <- runMCMC(pgSetup, sampler = "DEzs", settings = list(nrChains=3, iterations = 2000000, thin=100))
save.list <- list(ref.ba, ref.gpp, reference)
save(save.list, file = "./out/reference_threePGN_list.RData")