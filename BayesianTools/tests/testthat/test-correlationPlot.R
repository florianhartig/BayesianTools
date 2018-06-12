testthat::context("Test correlationPlot")

ll <- generateTestDensityMultiNormal(sigma = "no correlation")

bayesianSetup <- createBayesianSetup(likelihood = ll, lower = c(-10, -5, -7.5), upper = c(10, 7.5, 3))

settings = list(iterations = 2000, nrChains=2)

out_1 <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = list(iterations=2000))
out_2 <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
out_3 <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)




testthat::test_that("correlationPlot works for bayesianOutput", {
  testthat::expect_error(correlationPlot(out_1), NA)
  testthat::expect_error(correlationPlot(out_2), NA)
  testthat::expect_error(correlationPlot(out_3), NA)
  
})

testthat::test_that("correlationPlot works for various parameter combinations", {
  testthat::expect_error(correlationPlot(out_1, method = "spearman"), NA)
  testthat::expect_error(correlationPlot(out_1, method = "kendall"), NA)
  
  testthat::expect_error(correlationPlot(out_1, method = "spearman", density = "ellipse"), NA)
  testthat::expect_error(correlationPlot(out_1, method = "kendall", density = "corellipseCor"), NA)
  
  testthat::expect_error(correlationPlot(out_1, method = "spearman", density = "smooth", scaleCorText = F), NA)
  
  testthat::expect_error(correlationPlot(out_1, method = "kendall", density = "corellipseCor", whichParameters = c(1,2)), NA)
  
  testthat::expect_error(correlationPlot(out_1, method = "spearman", start = 200), NA)
  testthat::expect_error(correlationPlot(out_1, method = "kendall", start = 200, thin = 100), NA)
})