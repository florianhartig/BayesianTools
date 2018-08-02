## This example shows how to run and calibrate the VSEM model 

library(BayesianTools)

# Create input data for the model
PAR <- VSEMcreatePAR(1:1000)
plot(PAR, main = "PAR (driving the model)", xlab = "Day")

# load reference parameter definition (upper, lower prior)
refPars <- VSEMgetDefaults()
# this adds one additional parameter for the likelihood standard deviation (see below)
refPars[12,] <- c(2, 0.1, 4) 
rownames(refPars)[12] <- "error-sd"
head(refPars)

# create some simulated test data 
# generally recommended to start with simulated data before moving to real data
referenceData <- VSEM(refPars$best[1:11], PAR) # model predictions with reference parameters  
referenceData[,1] = 1000 * referenceData[,1] 
# this adds the error - needs to conform to the error definition in the likelihood
obs <- referenceData + rnorm(length(referenceData), sd = refPars$best[12])
oldpar <- par(mfrow = c(2,2))
for (i in 1:4) plotTimeSeries(observed = obs[,i], 
                              predicted = referenceData[,i], main = colnames(referenceData)[i])

# Best to program in a way that we can choose easily which parameters to calibrate
parSel = c(1:6, 12)

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

##########################################################################
# Sensitity analysis

plotSensitivity <- function(n, ...){
  lowS = refPars$lower[parSel]
  upS = refPars$upper[parSel]
  parSen <- seq(lowS[n], upS[n], len = 20) 
  refPar <- refPars$best[parSel]
  post <- rep(NA, 20) 
  for (i in 1:20){
    parS <- refPar
    parS[n] = parSen[i]
    post[i] = bayesianSetup$posterior$density(parS)
  }
  plot(x = parSen, y = post, main = rownames(refPars)[parSel][n], type = "l", col =
         "red", ...)#, ylim = c(-7500, -7000) abline(v = refPar[n])
}
opar <- par(mfrow = c(2,3))
for (i in 1:6) plotSensitivity(i)
par(opar)


library(sensitivity)

morrisOut <- morris(model = bayesianSetup$posterior$density, factors = rownames(refPars[parSel, ]), r = 200, design = list(type = "oat", levels = 5, grid.jump = 3), binf = refPars$lower[parSel], bsup = refPars$upper[parSel], scale = TRUE)
par(mfrow=c(1,1))
plot(morrisOut)


# So according to this analysis the most sensitive of the six parameters is Light use efficiency (LUE), and the least sensitive is the residence time of soil organic matter (tauS).

# Now for the unscaled morris, i.e. where we look at the sensitivity for the same amount of numerical change

morrisOutUnscaled <- morris(model = bayesianSetup$posterior$density, factors = rownames(refPars[parSel, ]), r = 200, design = list(type = "oat", levels = 5, grid.jump = 3), binf = refPars$lower[parSel], bsup = refPars$upper[parSel], scale = FALSE)
plot(morrisOutUnscaled)


##########################################################################
# Model fit

# settings for the sampler, iterations should be increased for real applicatoin

settings <- list(iterations = 10000, nrChains = 2)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)

## Not run: 

plot(out)
summary(out)
marginalPlot(out, prior = T)
gelmanDiagnostics(out) # should be below 1.05 for all parameters to demonstrate convergence 

# Posterior predictive simulations

# Create a prediction function
createPredictions <- function(par){
  # set the parameters that are not calibrated on default values 
  x = refPars$best
  x[parSel] = par
  predicted <- VSEM(x[1:11], PAR) # replace here VSEM with your model 
  return(predicted[,1] * 1000)
}

# Create an error function
createError <- function(mean, par){
  return(rnorm(length(mean), mean = mean, sd = par[7]))
}

# plot prior predictive distribution and prior predictive simulations
plotTimeSeriesResults(sampler = out, model = createPredictions, observed = obs[,1],
                      error = createError, prior = TRUE, main = "Prior predictive")

# plot posterior predictive distribution and posterior predictive simulations
plotTimeSeriesResults(sampler = out, model = createPredictions, observed = obs[,1],
                      error = createError, main = "Posterior predictive")


########################################################
# Demonstrating the updating of the prior from old posterior
# Note that it is usually more exact to rerun the MCMC 
# with all (old and new) data, instead of updating the prior
# because likely some information is lost when approximating the
# Posterior by a multivariate normal 

settings <- list(iterations = 5000, nrChains = 2)

out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)

plot(out)
correlationPlot(out, start = 1000)

newPrior = createPriorDensity(out, method = "multivariate",
                              eps = 1e-10,
                              lower = refPars$lower[parSel],
                              upper = refPars$upper[parSel], start= 1000)

bayesianSetup <- createBayesianSetup(likelihood = likelihood,
                                     prior = newPrior,
                                     names = rownames(refPars)[parSel] )

# check boundaries are correct set
bayesianSetup$prior$sampler() < refPars$lower[parSel]
bayesianSetup$prior$sampler() > refPars$upper[parSel]

# check prior looks similar to posterior
x = bayesianSetup$prior$sampler(2000)
correlationPlot(x, thin = F)

out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)

plot(out)
correlationPlot(out)

plotTimeSeriesResults(sampler = out,
                      model = createPredictions,
                      observed = obs[,1],
                      error = createError,
                      prior = F, main = "Posterior predictive")

plotTimeSeriesResults(sampler = out,
                      model = createPredictions,
                      observed = obs[,1],
                      error = createError,
                      prior = T, main = "Prior predictive")





## End(Not run)

par(oldpar)