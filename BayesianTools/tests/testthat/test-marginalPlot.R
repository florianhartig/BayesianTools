testthat::context("Test marginalPlot")

testMarginalPlot <- function (x) {
  testthat::expect_error(marginalPlot(x, singlePanel = T), NA)
  testthat::expect_error(marginalPlot(x, singlePanel = F), NA)
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = T, plotPrior = T, singlePanel = F, densityOnly = F), NA)
  testthat::expect_error(marginalPlot(x, scale = T, histogram = T, plotPrior = T, singlePanel = F, densityOnly = F), NA)
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = T, plotPrior = T, singlePanel = T, densityOnly = F), NA)
  testthat::expect_error(marginalPlot(x, scale = T, histogram = T, plotPrior = T, singlePanel = T, densityOnly = F), NA)
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = T, plotPrior = F, singlePanel = F, densityOnly = F), NA)
  testthat::expect_error(marginalPlot(x, scale = T, histogram = T, plotPrior = F, singlePanel = F, densityOnly = F), NA)
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = T, plotPrior = F, singlePanel = T, densityOnly = F), NA)
  testthat::expect_error(marginalPlot(x, scale = T, histogram = T, plotPrior = F, singlePanel = T, densityOnly = F), NA)
  
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = F, plotPrior = T, singlePanel = F, densityOnly = F), NA)
  testthat::expect_error(marginalPlot(x, scale = T, histogram = F, plotPrior = T, singlePanel = F, densityOnly = F), NA)
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = F, plotPrior = T, singlePanel = T, densityOnly = F), NA)
  testthat::expect_error(marginalPlot(x, scale = T, histogram = F, plotPrior = T, singlePanel = T, densityOnly = F), NA)
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = F, plotPrior = F, singlePanel = F, densityOnly = F), NA)
  testthat::expect_error(marginalPlot(x, scale = T, histogram = F, plotPrior = F, singlePanel = F, densityOnly = F), NA)
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = F, plotPrior = F, singlePanel = T, densityOnly = F), NA)
  testthat::expect_error(marginalPlot(x, scale = T, histogram = F, plotPrior = F, singlePanel = T, densityOnly = F), NA)
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = F, plotPrior = F, singlePanel = T, densityOnly = F, best = T, dens = F), NA)
  testthat::expect_error(marginalPlot(x, scale = F, histogram = F, plotPrior = F, singlePanel = T, densityOnly = F, best = F, dens = F), NA)
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = T, plotPrior = T, singlePanel = F, densityOnly = F, thin = 10), NA)
  testthat::expect_error(marginalPlot(x, scale = F, histogram = T, plotPrior = T, singlePanel = F, densityOnly = F, start = 100), NA)
}


ll <- function(x, sum = TRUE){
  if(sum) sum(dnorm(x, log = T))
  else dnorm(x, log = T)
} 

setup <- createBayesianSetup(ll, lower = c(-10,-10,-10), upper = c(10,10,10))
settings <- list(iterations=10000)

out <- runMCMC(setup)
out_mat <- suppressWarnings(getSample(out, numSamples = 10000))
out_coda <- suppressWarnings(getSample(out, numSamples = 10000, coda = T))

testthat::test_that("marginalPlot runs without throwing an error (bayesianOutput)", testMarginalPlot(out))
testthat::test_that("marginalPlot runs without throwing an error (matrix)", testMarginalPlot(out_mat))
testthat::test_that("marginalPlot runs without throwing an error (coda)", testMarginalPlot(out_coda))



