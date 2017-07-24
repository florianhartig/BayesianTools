context("tests getSample and getSample.list")

test_that("getSample works for mcmc", {
  # number of test runs
  n_tests = 100
  
  # VECTOR
  
  # mcmc coda = TRUE
  for (i in 1:n_tests) {
    # create random vector
    dat <- coda::mcmc(rnorm(max(rnorm(n = 1, mean = 1000, sd = 1000), 10), 10, 5))
    len = length(dat)
    
    # "realistic" sample/thin
    thin <- sample(1:(len/5), 1)
    
    res_gs <- getSample(dat, coda = TRUE, thin = thin)
    res_coda <- window(x = dat, start = 1, end = len, thin = thin)
    expect_equal(res_gs, res_coda)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning == number of samples
    r <- getSample(dat, coda = TRUE, thin = len)
    expect_true(r[1] == dat[1] && length(r) == 1)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning too large
    r <- suppressWarnings(getSample(dat, coda = TRUE, thin = sample((len+1):(len + len/5), 1)))
    expect_true(r[1] == dat[1] && length(r) == 1)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning to low
    r <- suppressWarnings(getSample(dat, coda = TRUE, thin = sample((-len/5):0, 1)))
    expect_equal(r, dat)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning == 1
    r <- getSample(dat, coda = TRUE, thin = 1)
    expect_equal(r, dat)
    expect_true(class(res_gs) == "mcmc")
  }
  
  # TODO: Add tests for coda=FALSE
  
  for (i in 1:n_tests) {
    # create random vector
    dat <- coda::mcmc(rnorm(max(rnorm(n = 1, mean = 1000, sd = 1000), 10), 10, 5))
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
  
  # coda = TRUE
  for (i in 1:n_tests) {
    # random number of columns
    n_col = sample(2:max_n_col, 1)
    
    # create matrix of random length
    n_row = round(max(rnorm(n = 1, mean = 1000, sd = 1000), 10))
    vec <- rnorm(n_row * n_col, 10, 5)
    dat <- coda::mcmc(matrix(vec, ncol = n_col))
    
    len = n_row
    
    # "realistic" sample/thin
    thin <- sample(1:(len/5), 1)
    
    res_gs <- getSample(dat, coda = TRUE, thin = thin)
    res_coda <- window(x = dat, start = 1, end = len, thin = thin)
    expect_equal(res_gs, res_coda)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning == number of samples
    r <- getSample(dat, coda = TRUE, thin = len)
    expect_true(r[1,] == dat[1,] && nrow(r) == 1)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning too large
    r <- suppressWarnings(getSample(dat, coda = TRUE, thin = sample((len+1):(len + len/5), 1)))
    expect_true(r[1,] == dat[1,] && nrow(r) == 1)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning to low
    r <- suppressWarnings(getSample(dat, coda = TRUE, thin = sample((-len/5):0, 1)))
    expect_equal(r, dat)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning == 1
    r <- getSample(dat, coda = TRUE, thin = 1)
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
    dat <- coda::mcmc(matrix(vec, ncol = n_col))
    
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

test_that("getSample works for mcmc.list", {
  
  n_tests = 100
  
  # mcmc.list of mcmc objects containing vectors
  for (i in 1:n_tests) {
    
    # create random vectors
    len = sample(1:1000, 1)
    a <- coda::mcmc(rnorm(len, 0, 5))
    b <- coda::mcmc(rnorm(len, 0, 5))
    c <- coda::mcmc(rnorm(len, 0, 5))
    
    # create mcmc.list
    dat <- coda::mcmc.list(list(a,b,c))
    
    # coda = FALSE
    
    # "realistic" thinning
    thin <- round(runif(1, 1, len/3))
    r <- getSample(dat, thin = thin, coda=FALSE)
    expect_true(is.vector(r))
    
    # thinning too large
    thin <- length(dat) * len + 1
    r <- suppressWarnings(getSample(dat, thin = thin, coda=FALSE))
    expect_true(is.vector(r))
    expect_true(length(r) == 1)
    
    #thinning too small
    thin <- sample(-10:0, 1)
    r <- suppressWarnings(getSample(dat, thin = thin, coda=FALSE))
    expect_true(is.vector(r))
    expect_true(length(r) == len *length(dat))
    
    
    # coda = TRUE
    
    # "realistic" thinning
    thin <- round(runif(1, 1, len/3))
    r <- getSample(dat, thin = thin, coda=TRUE)
    expect_true(class(r) == "mcmc.list")
    
    # thinning too large
    thin <- length(dat) * len + 1
    r <- suppressWarnings(getSample(dat, thin = thin, coda=TRUE))
    expect_true(class(r) == "mcmc.list")
    expect_true(length(dat) == Reduce(function(acc, val) acc + length(val), r, 0))
    
    #thinning too small
    thin <- sample(-10:0, 1)
    r <- suppressWarnings(getSample(dat, thin = thin, coda=TRUE))
    expect_true(class(r) == "mcmc.list")
    expect_true(len * length(dat) == Reduce(function(acc, val) acc + length(val), r, 0))
  }
  
  
  # mcmc.list of mcmc objects containing matrices
  for (i in 1:n_tests) {
    
    # create random matrices
    len = sample(1:1000, 1)
    n_col = sample(2:10, 1)
    
    a <- coda::mcmc(matrix(rnorm(len * n_col, 0, 5), ncol = n_col))
    b <- coda::mcmc(matrix(rnorm(len * n_col, 0, 5), ncol = n_col))
    c <- coda::mcmc(matrix(rnorm(len * n_col, 0, 5), ncol = n_col))
    
    # create mcmc.list
    dat <- coda::mcmc.list(list(a,b,c))
    
    # coda = FALSE
    
    # "realistic" thinning
    thin <- round(runif(1, 1, len/3))
    r <- getSample(dat, thin = thin, coda=FALSE)
    expect_true(is.matrix(r))
    
    # thinning too large
    thin <- length(dat) * len + 1
    r <- suppressWarnings(getSample(dat, thin = thin, coda=FALSE))
    expect_true(is.matrix(r))
    expect_true(nrow(r) == 1)
    
    #thinning too small
    thin <- sample(-10:0, 1)
    r <- suppressWarnings(getSample(dat, thin = thin, coda=FALSE))
    expect_true(is.matrix(r))
    expect_true(nrow(r) == len * length(dat))
    
    
    # coda = TRUE
    
    # "realistic" thinning
    thin <- round(runif(1, 1, len/3))
    r <- getSample(dat, thin = thin, coda=TRUE)
    expect_true(class(r) == "mcmc.list")
    
    # thinning too large
    thin <- length(dat) * len + 1
    r <- suppressWarnings(getSample(dat, thin = thin, coda=TRUE))
    expect_true(class(r) == "mcmc.list")
    expect_true(length(dat) == Reduce(function(acc, val) acc + nrow(val), r, 0))
    
    #thinning too small
    thin <- sample(-10:0, 1)
    r <- suppressWarnings(getSample(dat, thin = thin, coda=TRUE))
    expect_true(class(r) == "mcmc.list")
    expect_true(len * length(dat) == Reduce(function(acc, val) acc + nrow(val), r, 0))
  }
  
})









