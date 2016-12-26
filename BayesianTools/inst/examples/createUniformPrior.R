prior <- createUniformPrior(c(0,0),c(0.4,5))
prior$density(c(2,3))
prior$density(c(0.2,2))
