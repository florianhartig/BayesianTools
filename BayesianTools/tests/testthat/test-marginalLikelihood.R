context("Test marginalLikelihood")

skip_on_cran()

ll <- function(x, sum = TRUE){
  if(sum) sum(dnorm(x, log = T))
  else dnorm(x, log = T)
} 

setup <- createBayesianSetup(ll, lower = c(-10,-10,-10), upper = c(10,10,10))
settings <- list(iterations=10000)

out <- runMCMC(setup)

testthat::test_that("marginalLikelihood with method 'Chib' runs without error", {
  marginalLikelihood(sampler = out, method = "Chib")
})

testthat::test_that("marginalLikelihood with method 'HM' runs without error", {
  marginalLikelihood(sampler = out, method = "HM")
})

testthat::test_that("marginalLikelihood with method 'Bridge' runs without error", {
  marginalLikelihood(sampler = out, method = "Bridge")
})

testthat::test_that("marginalLikelihood with method 'Prior' runs without error", {
  marginalLikelihood(sampler = out, method = "Prior")
  marginalLikelihood(sampler = setup, method = "Prior")
})
