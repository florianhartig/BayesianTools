context("Test marginalLikelihood")

ll <- function(x, sum = TRUE){
  if(sum) sum(dnorm(x, log = T))
  else dnorm(x, log = T)
} 

setup <- createBayesianSetup(ll, lower = c(-10,-10,-10), upper = c(10,10,10))
settings <- list(iterations=10000)

out <- runMCMC(setup)

testthat::test_that("marginalLikelihood with method 'Chib' runs without error", {
  testthat::expect_error(marginalLikelihood(sampler = out, method = "Chib"), NA)
})

testthat::test_that("marginalLikelihood with method 'HM' runs without error", {
  testthat::expect_error(marginalLikelihood(sampler = out, method = "HM"), NA)
})

testthat::test_that("marginalLikelihood with method 'Bridge' runs without error", {
  testthat::expect_error(marginalLikelihood(sampler = out, method = "Chib"), NA)
})

testthat::test_that("marginalLikelihood with method 'prior' runs without error", {
  testthat::expect_error(marginalLikelihood(sampler = out, method = "Chib"), NA)
})
