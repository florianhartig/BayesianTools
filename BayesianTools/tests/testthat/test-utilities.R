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


