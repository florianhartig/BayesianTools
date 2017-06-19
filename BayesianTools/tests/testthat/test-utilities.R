context("Test utility functions")

# sampleEquallySpaced
test_that("sampleEquallySpaced returns right number of samples", {
  max_n_samples = 1000
  margin <- round(max_n_samples * 0.05)
  n_tests = 250
  
  # vectors
  for (i in 1:n_tests) {
    vec <- rep(1, sample(1:max_n_samples, 1))
    test_len <- sample((1 - margin):(max_n_samples + margin), 1)
    
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
  max_n_col = 10;
  max_n_samples = 1000
  margin <- round(max_n_samples * 0.05)
  n_tests = 250
  
  for (i in 1:n_tests) {
    n_col = sample(2:max_n_col, 1)
    vec <- rep(1, sample(1:max_n_samples, 1) * n_col)
    mat <- matrix(vec, ncol = n_col)
    
    test_len <- sample((1 - margin):(max_n_samples + margin), 1)
    
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


