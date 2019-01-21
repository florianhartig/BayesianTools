library(BayesianTools)

?VSEM

# Create input data for the model
PAR <- VSEMcreatePAR(1:1000)
plot(PAR, main = "PAR (driving the model)", xlab = "Day")

# load reference parameter definition (upper, lower prior)
refPars <- VSEMgetDefaults()
# this adds one additional parameter for the likelihood standard deviation (see below)
refPars[12,] <- c(2, 0.1, 4) 
rownames(refPars)[12] <- "error-sd"
head(refPars)

dat = list()


randomEffect = runif(10, -0.4,0.4)

for(j in 1:10){
  pars = refPars$best[1:11]
  pars[2] = pars[2] + randomEffect[j]
  referenceData <- VSEM(pars, PAR) # model predictions with reference parameters  
  referenceData[,1] = 1000 * referenceData[,1] 
  # this adds the error - needs to conform to the error definition in the likelihood
  obs <- referenceData + rnorm(length(referenceData), sd = refPars$best[12])
  dat[[j]] = obs
}

# Best to program in a way that we can choose easily which parameters to calibrate
parSel = c(1:6, 12)

# here is the likelihood 
singleLikelihood <- function(par, j, sum = TRUE){
  # set parameters that are not calibrated on default values 
  x = refPars$best
  x[parSel] = par[1:7]
  # take x1 out of the calibration later, this just overwrites x1
  x[1] = refPars$best[1]

  randomEffectx = par[9:18]

  # add the random slope to first parameters (which is the light extinction)
  x[2] = x[2] + randomEffectx[j]
  
  # Environment, par[8] ist just an example
  # x[2] = x[2] + par[8] * Env[j]
  
  
  predicted <- VSEM(x[1:11], PAR) # replace here VSEM with your model 
  predicted[,1] = 1000 * predicted[,1] # this is just rescaling
  diff <- c(predicted[,1:4] - dat[[j]][,1:4]) # difference betweeno observed and predicted
  # univariate normal likelihood. Note that there is a parameter involved here that is fit
  llValues <- dnorm(diff, sd = x[12], log = TRUE)  
  if (sum == FALSE) return(llValues)
  else return(sum(llValues))
}

likelihood<-function(par, sum = TRUE){
  count = 0
  
  for(i in 1:10){
   localpar[2] = -1 + 2 * logit( par[2] + par[8] * plotData$Temp[i] + par[9] * plotData$Precip[i] )
    
   count = count + singleLikelihood(localpar, sum = sum)
  }

  rsigmax = par[8]
  randomEffectx = par[9:18]
  
  llRandom = sum(dnorm(randomEffectx, sd = rsigmax, log = T))
  
  # Spatial Random effects, correlationMatrix sigma must be made dependent on the distance, with some decay assumptions, and then we fit the parameter for the decay
  # https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/05-HierarchicalAndSpatialModels/SpatialAutoregressive/gaussianAutoregressiveExample.md
  library(mvtnorm)
  llRandom = sum(dmvnorm(randomEffectx, sigma = invDistance, log = T))
  
  return(count + llRandom)  
}

likelihood(c(refPars$best[parSel], 0.5, rep(0, 10)))


# optional, you can also directly provide lower, upper in the createBayesianSetup, see help
prior <- createUniformPrior(lower = c(refPars$lower[parSel], 0, rep(-0.4, 10)) , 
                            upper = c(refPars$upper[parSel], 2, rep(0.4, 10)), best = c(refPars$best[parSel], 0.5, rep(0, 10)))

bayesianSetup <- createBayesianSetup(likelihood, prior)

# settings for the sampler, iterations should be increased for real applicatoin
settings <- list(iterations = 5000, nrChains = 2)

out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)

plot(out)
marginalPlot(out)
summary(out, start = 500)

test <- getSample(out, start = 500)
test2 = colMeans(test)
plot(test2[9:18], randomEffect)




