prior <- createTruncatedNormalPrior(c(0,0),c(0.4,5), lower = c(-2,-2), upper = c(1,1))
prior$density(c(2,3))
prior$density(c(0.2,0.9))
prior$sampler()
