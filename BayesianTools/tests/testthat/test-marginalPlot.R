testthat::context("Test marginalPlot")

testMarginalPlot <- function (x) {
  if(!'bayesianOutput' %in% class(x)) prior <- FALSE else prior <- TRUE
  
  testthat::expect_error(marginalPlot(x, type = 'd', singlePanel = T, prior = prior), NA)
  testthat::expect_error(marginalPlot(x, type = 'd', singlePanel = F, prior = prior), NA)
  
  testthat::expect_error(marginalPlot(x, type = 'v', singlePanel = T, prior = prior), NA)
  testthat::expect_error(marginalPlot(x, type = 'v', singlePanel = F, prior = prior), NA)
  
  if ('bayesianOutput' %in% class(x)) {
    testthat::expect_error(marginalPlot(x, type = 'v', singlePanel = T, prior = T), NA)
    testthat::expect_error(marginalPlot(x, type = 'v', singlePanel = T, prior = F), NA)
    testthat::expect_error(marginalPlot(x, type = 'd', singlePanel = T, prior = T), NA)
    testthat::expect_error(marginalPlot(x, type = 'd', singlePanel = T, prior = F), NA)
  }
  
  testthat::expect_error(marginalPlot(x, type = 'v', xrange = c(-5, 5), prior = prior), NA)
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



