testthat::context("Test plotDiagnostic")
bayesianSetup <- createBayesianSetup(likelihood = testDensityNormal, prior = createUniformPrior(lower = rep(-10,2), upper = rep(10,2)))

# settings = list(nrChains = 1, iterations = 100)
# outOneChain = runMCMC(bayesianSetup = bayesianSetup, settings = list(nrChains = 1, iterations = 100))
# outTwoChain = runMCMC(bayesianSetup = bayesianSetup, settings = list(nrChains = 2, iterations = 100))
# 
# outTwoChainDiffSamp = runMCMC(bayesianSetup = bayesianSetup, settings = list(nrChains = 2, iterations = 100), sampler = "DREAM")


testthat::test_that("plotDiagnostic works for bayesianOutput", {
  testthat::expect_error(plotDiagnostic(matrix(0, ncol = 4, nrow = 2)))
})