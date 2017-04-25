library(BayesianTools)

PAR <- VSEMcreatePAR(1:1000)
plotTimeSeries(observed = PAR)


refPars$upper = refPars$upper * 10

refPars <- VSEMgetDefaults()
refPars[12,] <- c(0.1, 0.001, 0.5)
rownames(refPars)[12] <- "error-sd"
head(refPars)

referenceData <- VSEM(refPars$best[1:11], PAR) 
referenceData[,1] = 1000 * referenceData[,1]
obs <- referenceData + rnorm(length(referenceData), sd = refPars$best[12])
oldpar <- par(mfrow = c(2,2))
for (i in 1:4) plotTimeSeries(observed = obs[,i], predicted = referenceData[,i], main = colnames(referenceData)[i])


parSel = c(1:12)

likelihood <- function(x, sum = TRUE){
  x <- createMixWithDefaults(x, refPars$best, parSel)
  predicted <- VSEM(x[1:11], PAR)
  predicted[,1] = 1000 * predicted[,1]
  diff <- c(predicted[,1:4] - obs[,1:4])
  llValues <- dnorm(diff, sd = x[12], log = T) 
  # if(runif(1) < 0.4) stop()
  if (sum == FALSE) return(llValues)
  else return(sum(llValues))
}


prior <- createUniformPrior(lower = refPars$lower[parSel], upper = refPars$upper[parSel], best = refPars$best[parSel])


bayesianSetup <- createBayesianSetup(likelihood, prior, names = rownames(refPars)[parSel])


settings <- list(iterations = 10000, nrChains = 3)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DREAMzs", settings = settings)
getSample(out, parametersOnly = F)[,13:15]



plot(out)
summary(out)

