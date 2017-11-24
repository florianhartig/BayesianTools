testthat::context("Test plotTimeSeriesResults")

# see VSEMHelp
PAR <- VSEMcreatePAR(1:1000)
refPars <- VSEMgetDefaults()
refPars[12,] <- c(2, 0.1, 4) 
rownames(refPars)[12] <- "error-sd"
referenceData <- VSEM(refPars$best[1:11], PAR) # model predictions with reference parameters  
referenceData[,1] = 1000 * referenceData[,1] 
obs <- referenceData + rnorm(length(referenceData), sd = refPars$best[12])
parSel = c(1:6, 12)

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

prior <- createUniformPrior(lower = refPars$lower[parSel], 
                            upper = refPars$upper[parSel], best = refPars$best[parSel])

createPredictions <- function(par){
  # set the parameters that are not calibrated on default values 
  x = refPars$best
  x[parSel] = par
  predicted <- VSEM(x[1:11], PAR) # replace here VSEM with your model 
  return(predicted[,1] * 1000)
}

createError <- function(mean, par){
  return(rnorm(length(mean), mean = mean, sd = par[7]))
}

bayesianSetup <- createBayesianSetup(likelihood, prior, names = rownames(refPars)[parSel])
settings <- list(iterations = 1000, nrChains = 2)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
out_mat <- getSample(out, numSamples = 2000)
out_coda <- getSample(out, numSamples = 2000)

testthat::test_that("plotTimeSeriesResults works for bayesianOutput, matrix and coda", {
  testthat::expect_error(plotTimeSeriesResults(sampler = out, model = createPredictions, observed = referenceData[,1],
                                               error = createError, prior = TRUE, main = "Prior predictive"),
                         NA)
  testthat::expect_error(plotTimeSeriesResults(sampler = out, model = createPredictions, observed = referenceData[,1],
                                               error = createError, prior = FALSE, main = "Prior predictive"),
                         NA)
  
  testthat::expect_error(plotTimeSeriesResults(sampler = out_mat, model = createPredictions, observed = referenceData[,1],
                                               error = createError, prior = TRUE, main = "Prior predictive"))
  testthat::expect_error(plotTimeSeriesResults(sampler = out, model = createPredictions, observed = referenceData[,1],
                                               error = createError, prior = FALSE, main = "Prior predictive"),
                         NA)
  
  testthat::expect_error(plotTimeSeriesResults(sampler = out_coda, model = createPredictions, observed = referenceData[,1],
                                               error = createError, prior = TRUE, main = "Prior predictive"))
  testthat::expect_error(plotTimeSeriesResults(sampler = out_coda, model = createPredictions, observed = referenceData[,1],
                                               error = createError, prior = FALSE, main = "Prior predictive"),
                         NA)
})
