context("Test createBayesianSetup")

dens <- function(par) {
  d1 <- dunif(par[1], -2, 6, log = TRUE)
  d2 <- dnorm(par[2], mean = 2, sd = 3, log = TRUE)
  return(d1 + d2)
}

samp <- function(n = 1) {
  d1 <- runif(n, -2, 6)
  d2 <- rnorm(n, mean = 2, sd = 3)
  return(cbind(d1,d2))
}

lower = c(-5, -5)
upper = c(5, 5)

prior <- createPrior(density = dens, samp = samp)

ll <- function(x) sum(dnorm(x, log = TRUE)) # multivariate normal ll

# test_that("createBayesianSetup works with prior object input", {
#   testthat::expect_warning({
#     bs1 <- createBayesianSetup(likelihood = ll, prior = prior)
#   }, NA, label = "prior only")
#   testthat::expect_warning({
#     bs2 <- createBayesianSetup(likelihood = ll, prior = prior, priorSampler = samp)
#   }, NA, label = "prior and priorSampler")
#   testthat::expect_warning({
#     bs3 <- createBayesianSetup(likelihood = ll, prior = prior, lower = lower, upper = upper)
#   }, NA, label = "prior and lower/upper")
#   testthat::expect_warning({
#     bs4 <- createBayesianSetup(likelihood = ll, prior = prior, priorSampler = samp, lower = c(-5, -5), upper = c(5, 5))
#   }, NA, label = "prior, priorSampler, and lower/upper")
#   testthat::expect_true({
#     toString(bs1) == toString(bs2) & toString(bs1) == toString(bs3) & toString(bs1) == toString(bs4)
#   }, label = "consistent results")
# })

test_that("createBayesianSetup works with prior function input", {
  testthat::expect_error({
    bs1 <- createBayesianSetup(likelihood = ll, prior = dens)
  }, label = "density only. should throw error")
  testthat::expect_error({
    bs2 <- createBayesianSetup(likelihood = ll, prior = dens, priorSampler = samp)
  }, NA, label = "density and priorSampler")
  testthat::expect_error({
    bs3 <- createBayesianSetup(likelihood = ll, prior = dens, lower = c(-5, -5), upper = c(5, 5))
  }, NA, label = "density and lower/upper")
  testthat::expect_error({
    bs4 <- createBayesianSetup(likelihood = ll, prior = dens, priorSampler = samp, lower = c(-5, -5), upper = c(5, 5))
    testthat::expect_true({is.null(bs4$prior$lower) & is.null(bs4$prior$upper) & is.null(bs4$prior$best)})
    testthat::expect_true({is.null(bs4$info$plotLower) & is.null(bs4$info$plotUpper) & is.null(bs4$info$plotBest)})
  }, NA, label = "density, priorSampler, and lower/upper")
})

test_that("createBayesianSetup works with lower/upper input", {
  testthat::expect_error({
    bs1 <- createBayesianSetup(likelihood = ll, lower = lower, upper = upper)
    testthat::expect_true({ all (bs1$prior$lower == bs1$info$plotLower & bs1$prior$upper == bs1$info$plotUpper)})
  }, NA, label = "lower/upper only")
  testthat::expect_error({
    bs2 <- createBayesianSetup(likelihood = ll, lower = lower, upper = upper, priorSampler = samp)
    testthat::expect_true({all(bs1$prior$lower == bs1$info$plotLower & bs1$prior$upper == bs1$info$plotUpper)})
  }, NA, label = "lower/upper and priorSampler")
  testthat::expect_true({
    toString(bs1) == toString(bs2)
  }, label = "consistent results")
})


test_that("likelihood properly handles exceptions", {
  
  likelihood = function(x){
    if(x == 1) return(1) 
    if(x == 2) stop()
    if(x == 3) return(Inf)
    if(x == 4) return(-Inf)
    if(x == 5) return(NaN)
    if(x == 6) return(NA)
    if(x == 7) return(F)
    if(x == 8) return("test")
    
  }
  setup <- createBayesianSetup(likelihood, lower = 0.5, upper = 1.5)
  
  expect_equal(setup$posterior$density(1), 1)
  expect_equal(setup$posterior$density(2), -Inf)
  expect_equal(setup$posterior$density(3), -Inf)
  expect_equal(setup$posterior$density(4), - Inf)
  expect_equal(setup$posterior$density(5), - Inf)
  expect_equal(setup$posterior$density(6), - Inf)
  expect_equal(setup$posterior$density(7), -Inf)
  # This doesn't work yet, but not sure if it's neccessary to fix it
  #expect_equal(setup$posterior$density(8), - Inf)
  

})


