context("Test utility functions")

# sampleEquallySpaced
test_that("sampleEquallySpaced returns right number of samples", {
  # number of tests
  n_tests = 250
  
  # maximumum number of samples
  max_n_samples = 1000
  
  # the target length should sometimes be lower or higher than the max_n_samples to simulate errors
  # the higher the margin the more of these errors will occur
  margin <- round(max_n_samples * 0.05)
  
  # vectors
  for (i in 1:n_tests) {
    # create vector of random length
    vec <- rep(1, sample(1:max_n_samples, 1))
    # create random target length for sampling
    test_len <- sample((1 - margin):(max_n_samples + margin), 1)
    
    # try to sample vector to target length test_len
    result <- suppressWarnings(sampleEquallySpaced(vec, test_len))
    
    if (test_len <= 1) {
      expect_true(length(result) == 1)
    } else if (test_len >= length(vec)) {
      expect_true(length(result) == length(vec))
    } else {
      expect_true(length(result) == test_len)
    }
  }
  
  # matrices
  
  # max number of columns in test matrix
  max_n_col = 10;
  
  for (i in 1:n_tests) {
    # random number of columns
    n_col = sample(2:max_n_col, 1)
    # create matrix of random length
    vec <- rep(1, sample(1:max_n_samples, 1) * n_col)
    mat <- matrix(vec, ncol = n_col)
    
    # create random target length for sampling
    test_len <- sample((1 - margin):(max_n_samples + margin), 1)
    
    #try to sample matrix to target length test_len
    result <- suppressWarnings(sampleEquallySpaced(mat, test_len))
    
    if (test_len <= 1) {
      expect_true(length(result) == n_col)
    } else if (test_len >= nrow(mat)) {
      expect_true(nrow(result) == nrow(mat))  
    } else {
      expect_true(nrow(result) == test_len)
    }
  }
})


# test various prior/bayesianSetup print outputs
test_that("test bayesianSetup and prior print functions",{
  
  prior = createPrior(sampler =  function(n=1) return(cbind(rnorm(n),rnorm(n), rnorm(n))))
  ll = testDensityMultiNormal
  
  expect_output(print(createBayesianSetup(likelihood = ll,  lower = c(-10,3), upper =c(10,3))))
  expect_output(print(createBayesianSetup(likelihood = ll, prior = prior )))
  expect_output(print(createBayesianSetup(likelihood = ll, prior = prior, names = c("A","B","C"))))
  expect_error(print(createBayesianSetup(likelihood = ll)))
  expect_output(print(createPrior(sampler =  function(n=1) return(cbind(rnorm(n),rnorm(n), rnorm(n))))))
  expect_output(print(createUniformPrior(lower = c(0,0), upper = c(0,5))))
  expect_output(print(createTruncatedNormalPrior(c(0,0),c(0.4,5), lower = c(-2,-2), upper = c(1,1))))
  
})


