testthat::context("Test marginalPlot")

testMarginalPlot <- function (x, priorSample) {

  marginalPlot(x, type = 'v', singlePanel = F)
  marginalPlot(x, type = 'v', singlePanel = T)
  marginalPlot(x, type = 'd', singlePanel = F)
  marginalPlot(x, type = 'd', singlePanel = T)
    
  marginalPlot(x, type = 'v', singlePanel = F, prior = F)
  marginalPlot(x, type = 'v', singlePanel = F, prior = T)
  marginalPlot(x, type = 'v', singlePanel = F, prior = priorSample)

  expect_error(marginalPlot(x, type = 'v', singlePanel = F, prior = c(1,2)))
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

priorSample = setup$prior$sampler(10000)

testthat::test_that("marginalPlot runs without throwing an error (bayesianOutput)", testMarginalPlot(out, priorSample))
testthat::test_that("marginalPlot runs without throwing an error (matrix)", testMarginalPlot(out_mat, priorSample))
testthat::test_that("marginalPlot runs without throwing an error (data.frame)", testMarginalPlot(data.frame(out_mat), priorSample))
testthat::test_that("marginalPlot runs without throwing an error (coda)", testMarginalPlot(out_coda, priorSample))

marginalPlot(out,  type = "d", xrange = c(0,1))
marginalPlot(out,  type = "d", prior = F)


