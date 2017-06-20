context("tests getSample and getSample.list")

test_that("getSample works for mcmc", {
  # number of test runs
  n_tests = 100
  
  # VECTOR
  
  # mcmc coda = FALSE
  for (i in 1:n_tests) {
    # create random vector
    dat <- coda::mcmc(rnorm(max(rnorm(n = 1, mean = 1000, sd = 1000), 10), 10, 5))
    len = length(dat)
    
    # "realistic" sample/thin
    thin <- sample(1:(len/5), 1)
    
    res_gs <- getSample(dat, coda = FALSE, thin = thin)
    res_coda <- window(x = dat, start = 1, end = len, thin = thin)
    expect_equal(res_gs, res_coda)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning == number of samples
    r <- getSample(dat, coda = FALSE, thin = len)
    expect_true(r[1] == dat[1] && length(r) == 1)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning too large
    r <- suppressWarnings(getSample(dat, coda = FALSE, thin = sample((len+1):(len + len/5), 1)))
    expect_true(r[1] == dat[1] && length(r) == 1)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning to low
    r <- suppressWarnings(getSample(dat, coda = FALSE, thin = sample((-len/5):0, 1)))
    expect_equal(r, dat)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning == 1
    r <- getSample(dat, coda = FALSE, thin = 1)
    expect_equal(r, dat)
    expect_true(class(res_gs) == "mcmc")
  }
  
  # TODO: Add tests for coda=FALSE
  
  for (i in 1:n_tests) {
    # create random vector
    dat <- mcmc(rnorm(max(rnorm(n = 1, mean = 1000, sd = 1000), 10), 10, 5))
    len = length(dat)
    
    # "realistic" sample
    thin <- sample(1:(len/5), 1)
    
    res_gs <- getSample(dat, coda = FALSE, thin = thin)
    res_coda <- window(x = dat, start = 1, end = len, thin = thin)
    expect_equal(length(res_gs), length(res_coda))
    expect_true(class(res_gs) == "numeric")
    
    # thinning == number of samples
    r <- getSample(dat, coda = FALSE, thin = len)
    expect_true(r[1] == dat[1])
    expect_true(class(r) == "numeric")
    
    # thinning too large
    r <- suppressWarnings(getSample(dat, coda = FALSE, thin = sample((len+1):(len + len/5), 1)))
    expect_true(r[1] == dat[1])
    expect_true(class(r) == "numeric")
    
    # thinning to low
    r <- suppressWarnings(getSample(dat, coda = FALSE, thin = sample((-len/5):0, 1)))
    expect_equal(length(r), length(dat))
    expect_true(class(r) == "numeric")
    
    # thinning == 1
    r <- getSample(dat, coda = FALSE, thin = 1)
    expect_equal(length(r), length(dat))
    expect_true(class(r) == "numeric")
  }
  
  
  # MATRIX
  
  # max number of columns in test matrix
  max_n_col = 10;
  
  # coda = FALSE
  for (i in 1:n_tests) {
    # random number of columns
    n_col = sample(2:max_n_col, 1)
    
    # create matrix of random length
    n_row = round(max(rnorm(n = 1, mean = 1000, sd = 1000), 10))
    vec <- rnorm(n_row * n_col, 10, 5)
    dat <- mcmc(matrix(vec, ncol = n_col))
    
    len = n_row
    
    # "realistic" sample/thin
    thin <- sample(1:(len/5), 1)
    
    res_gs <- getSample(dat, coda = FALSE, thin = thin)
    res_coda <- window(x = dat, start = 1, end = len, thin = thin)
    expect_equal(res_gs, res_coda)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning == number of samples
    r <- getSample(dat, coda = FALSE, thin = len)
    expect_true(r[1,] == dat[1,] && nrow(r) == 1)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning too large
    r <- suppressWarnings(getSample(dat, coda = FALSE, thin = sample((len+1):(len + len/5), 1)))
    expect_true(r[1,] == dat[1,] && nrow(r) == 1)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning to low
    r <- suppressWarnings(getSample(dat, coda = FALSE, thin = sample((-len/5):0, 1)))
    expect_equal(r, dat)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning == 1
    r <- getSample(dat, coda = FALSE, thin = 1)
    expect_equal(r, dat)
    expect_true(class(res_gs) == "mcmc")
  }
  
  # coda = FALSE
  for (i in 1:n_tests) {
    # random number of columns
    n_col = sample(1:max_n_col, 1)
    
    # create matrix of random length
    n_row = round(max(rnorm(n = 1, mean = 1000, sd = 1000), 10))
    vec <- rnorm(n_row * n_col, 10, 5)
    dat <- mcmc(matrix(vec, ncol = n_col))
    
    len = n_row
    
    # "realistic" sample/thin
    thin <- sample(1:(len/5), 1)
    
    res_gs <- getSample(dat, coda = FALSE, thin = thin)
    res_coda <- window(x = dat, start = 1, end = len, thin = thin)
    expect_equal(length(res_gs), length(res_coda))
    expect_true(class(res_gs) == "matrix")
    
    # thinning == number of samples
    r <- getSample(dat, coda = FALSE, thin = len)
    expect_true(r[1,] == dat[1,] && nrow(r) == 1)
    expect_true(class(res_gs) == "matrix")
    
    # thinning too large
    r <- suppressWarnings(getSample(dat, coda = FALSE, thin = sample((len+1):(len + len/5), 1)))
    expect_true(r[1,] == dat[1,] && nrow(r) == 1)
    expect_true(class(res_gs) == "matrix")
    
    # thinning to low
    r <- suppressWarnings(getSample(dat, coda = FALSE, thin = sample((-len/5):0, 1)))
    expect_equal(length(r), length(dat))
    expect_true(class(res_gs) == "matrix")
    
    # thinning == 1
    r <- getSample(dat, coda = FALSE, thin = 1)
    expect_equal(length(r), length(dat))
    expect_true(class(res_gs) == "matrix")
  }  
})
  

## mcmc.list

# TODO: Add tests for vector

# TODO: Add tests for coda=FALSE
# TODO: Add tests for coda=FALSE

# TODO: Add tests for matrix
# TODO: Add tests for coda=FALSE
# TODO: Add tests for coda=FALSE