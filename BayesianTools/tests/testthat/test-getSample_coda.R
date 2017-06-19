context("tests getSample.mcmc and getSample.mcmc.list")

test_that("getSample works for mcmc and mcmc.list", {
  
  skip()
  
  n_tests = 100
  
  ## mcmc
  
  # vector
  
  # mcmc coda=TRUE
  for (i in 1:n_tests) {
    dat <- mcmc(rnorm(max(rnorm(n = 1, mean = 1000, sd = 1000), 10), 10, 5))
    len = length(dat)
    
    # "realistic" sample
    thin <- sample(1:(len/5), 1)
    
    res_gs <- getSample.mcmc(mcmc(dat), coda = TRUE, thin = thin)
    res_coda <- window(x = dat, start = 1, end = len, thin = thin)
    expect_equal(res_gs, res_coda)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning too large
    r <- suppressWarnings(getSample.mcmc(mcmc(dat), coda = TRUE, thin = sample((len+1):(len + len/5), 1)))
    expect_true(r[1] == dat[1])
    expect_true(class(res_gs) == "mcmc")
    
    # thinning == number of samples
    getSample.mcmc(mcmc(dat), coda = TRUE, thin = len)
    expect_true(r[1] == dat[1])
    expect_true(class(res_gs) == "mcmc")
    
    # thinning to low
    r <- suppressWarnings(getSample.mcmc(mcmc(dat), coda = TRUE, thin = sample((-len/5):0, 1)))
    expect_equal(r, dat)
    expect_true(class(res_gs) == "mcmc")
    
    # thinning == 1
    r <- getSample.mcmc(mcmc(dat), coda = TRUE, thin = 1)
    expect_equal(r, dat)
    expect_true(class(res_gs) == "mcmc")
  }
  
  # TODO: Add tests for coda=FALSE
  
  # TODO: Add tests for matrix

  
  ## mcmc.list
  
  # TODO: Add tests for vector
  
  # TODO: Add tests for coda=TRUE
  # TODO: Add tests for coda=FALSE
  
  # TODO: Add tests for matrix
  # TODO: Add tests for coda=TRUE
  # TODO: Add tests for coda=FALSE
})
  