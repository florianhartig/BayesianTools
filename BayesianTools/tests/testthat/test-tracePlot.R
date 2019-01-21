testthat::context("Test tracePlot")

skip_on_cran()

ll <- generateTestDensityMultiNormal(sigma = "no correlation")

bayesianSetup <- createBayesianSetup(likelihood = ll, lower = c(-10, -5, -7.5), upper = c(10, 7.5, 3))

settings = list(iterations = 2000, nrChains=2)

out_1 <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = list(iterations=2000))
out_2 <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
out_3 <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)

coda_1 <- getSample(out_1, coda = T)
coda_2 <- getSample(out_2, coda = T)
coda_3 <- getSample(out_3, coda = T)

mat_1 <- getSample(out_1, coda = F)
mat_2 <- getSample(out_2, coda = F)
mat_3 <- getSample(out_3, coda = F)


testthat::test_that("tracePlot works for bayesianOutput", {
  testthat::expect_error(tracePlot(out_1), NA)
  testthat::expect_error(tracePlot(out_2), NA)
  testthat::expect_error(tracePlot(out_3), NA)
})

testthat::test_that("tracePlot works for coda", {
  testthat::expect_error(tracePlot(coda_1), NA)
  testthat::expect_error(tracePlot(coda_2), NA)
  testthat::expect_error(tracePlot(coda_3), NA)
})

testthat::test_that("tracePlot works for matrix", {
  testthat::expect_error(tracePlot(coda_1), NA)
  testthat::expect_error(tracePlot(coda_2), NA)
  testthat::expect_error(tracePlot(coda_3), NA)
})