context("Test marginalLikelihood")

skip_on_cran()

ll <- function(x, sum = TRUE){
  if(sum) sum(dnorm(x, log = T))
  else dnorm(x, log = T)
} 

setup <- createBayesianSetup(ll, lower = c(-10,-10,-10), upper = c(10,10,10))
settings <- list(iterations=1000)

out <- runMCMC(setup)

testthat::test_that("marginalLikelihood with method 'Chib' runs without error", {
  expect_lt(marginalLikelihood(sampler = out, method = "Chib")$ln.ML, 100)
})

testthat::test_that("marginalLikelihood with method 'HM' runs without error", {
  expect_lt(suppressWarnings(marginalLikelihood(sampler = out, method = "HM")$ln.ML), 100)
})

testthat::test_that("marginalLikelihood with method 'Bridge' runs without error", {
  expect_lt(suppressWarnings(marginalLikelihood(sampler = out, method = "Bridge")$ln.ML), 100)
})

testthat::test_that("marginalLikelihood with method 'Prior' runs without error", {
  expect_lt(marginalLikelihood(sampler = out, method = "Prior")$ln.ML, 100)
  expect_lt(marginalLikelihood(sampler = setup, method = "Prior")$ln.ML, 100)
})

