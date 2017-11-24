context("Test marginalPlot")

testthat::test_that("marginalPlot runs without throwing an error",
{
  ll <- function(x, sum = TRUE){
    if(sum) sum(dnorm(x, log = T))
    else dnorm(x, log = T)
  } 
  
  setup <- createBayesianSetup(ll, lower = c(-10,-10,-10), upper = c(10,10,10))
  settings <- list(iterations=10000)
  
  out <- runMCMC(setup)
  
  testthat::expect_error(marginalPlot(out, scale = F, histogram = T, plotPrior = T, singlePanel = F), NA)
  testthat::expect_error(marginalPlot(out, scale = T, histogram = T, plotPrior = T, singlePanel = F), NA)
  
  testthat::expect_error(marginalPlot(out, scale = F, histogram = T, plotPrior = T, singlePanel = T), NA)
  testthat::expect_error(marginalPlot(out, scale = T, histogram = T, plotPrior = T, singlePanel = T), NA)
  
  testthat::expect_error(marginalPlot(out, scale = F, histogram = T, plotPrior = F, singlePanel = F), NA)
  testthat::expect_error(marginalPlot(out, scale = T, histogram = T, plotPrior = F, singlePanel = F), NA)
  
  testthat::expect_error(marginalPlot(out, scale = F, histogram = T, plotPrior = F, singlePanel = T), NA)
  testthat::expect_error(marginalPlot(out, scale = T, histogram = T, plotPrior = F, singlePanel = T), NA)
  
  
  testthat::expect_error(marginalPlot(out, scale = F, histogram = F, plotPrior = T, singlePanel = F), NA)
  testthat::expect_error(marginalPlot(out, scale = T, histogram = F, plotPrior = T, singlePanel = F), NA)
  
  testthat::expect_error(marginalPlot(out, scale = F, histogram = F, plotPrior = T, singlePanel = T), NA)
  testthat::expect_error(marginalPlot(out, scale = T, histogram = F, plotPrior = T, singlePanel = T), NA)
  
  testthat::expect_error(marginalPlot(out, scale = F, histogram = F, plotPrior = F, singlePanel = F), NA)
  testthat::expect_error(marginalPlot(out, scale = T, histogram = F, plotPrior = F, singlePanel = F), NA)
  
  testthat::expect_error(marginalPlot(out, scale = F, histogram = F, plotPrior = F, singlePanel = T), NA)
  testthat::expect_error(marginalPlot(out, scale = T, histogram = F, plotPrior = F, singlePanel = T), NA)
  
})

