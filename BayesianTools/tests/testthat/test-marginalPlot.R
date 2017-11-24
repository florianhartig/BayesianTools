testthat::context("Test marginalPlot")

testMarginalPlot <- function (x) {
  testthat::expect_error(marginalPlot(x, scale = F, histogram = T, plotPrior = T, singlePanel = F), NA)
  testthat::expect_error(marginalPlot(x, scale = T, histogram = T, plotPrior = T, singlePanel = F), NA)
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = T, plotPrior = T, singlePanel = T), NA)
  testthat::expect_error(marginalPlot(x, scale = T, histogram = T, plotPrior = T, singlePanel = T), NA)
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = T, plotPrior = F, singlePanel = F), NA)
  testthat::expect_error(marginalPlot(x, scale = T, histogram = T, plotPrior = F, singlePanel = F), NA)
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = T, plotPrior = F, singlePanel = T), NA)
  testthat::expect_error(marginalPlot(x, scale = T, histogram = T, plotPrior = F, singlePanel = T), NA)
  
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = F, plotPrior = T, singlePanel = F), NA)
  testthat::expect_error(marginalPlot(x, scale = T, histogram = F, plotPrior = T, singlePanel = F), NA)
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = F, plotPrior = T, singlePanel = T), NA)
  testthat::expect_error(marginalPlot(x, scale = T, histogram = F, plotPrior = T, singlePanel = T), NA)
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = F, plotPrior = F, singlePanel = F), NA)
  testthat::expect_error(marginalPlot(x, scale = T, histogram = F, plotPrior = F, singlePanel = F), NA)
  
  testthat::expect_error(marginalPlot(x, scale = F, histogram = F, plotPrior = F, singlePanel = T), NA)
  testthat::expect_error(marginalPlot(x, scale = T, histogram = F, plotPrior = F, singlePanel = T), NA)
}


ll <- function(x, sum = TRUE){
  if(sum) sum(dnorm(x, log = T))
  else dnorm(x, log = T)
} 

setup <- createBayesianSetup(ll, lower = c(-10,-10,-10), upper = c(10,10,10))
settings <- list(iterations=10000)

out <- runMCMC(setup)
out_mat <- getSample(out, numSamples = 10000)
out_coda <- getSample(out, numSamples = 10000, coda = T)

testthat::test_that("marginalPlot runs without throwing an error (bayesianOutput)", testMarginalPlot(out))
testthat::test_that("marginalPlot runs without throwing an error (matrix)", testMarginalPlot(out_mat))
testthat::test_that("marginalPlot runs without throwing an error (matrix)", testMarginalPlot(out_coda))



