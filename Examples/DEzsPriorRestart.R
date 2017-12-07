library(BayesianTools)

# This example show how DEzs acceptance can run very low when starting the sampler with a prior that is extremely wide compared to the posterior. I think what's happening is that the z Matrix is started up from the prior, leading to a lot of bad values in the Z matrix. As soon as the sampler gets some good values, there is very little acceptance, with the result that the z-matrix is getting better values, but always the same, which leads again to bad mixing. 

# I wonder if a cure could be to add some noice on the values in the z-matrix, potentially proportional to the current posterior width or something like that?

ll <- generateTestDensityMultiNormal()

# I tried to create the situation where we have very wide priors
bayesianSetup <- createBayesianSetup(likelihood = ll, lower = rep(-10000, 3), upper = rep(10000, 3))
out <- runMCMC(bayesianSetup = bayesianSetup)
plot(out)
# This actually succeeds in creating low acceptance in the final chain
plot(out, start = 1500)

# Just to show that this is a prior problem
bayesianSetup2 <- createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))
out2 <- runMCMC(bayesianSetup = bayesianSetup2)
plot(out2)
plot(out2, start = 1500)

# Re-starting DE sampler

# OK, the idea is to re-start the first sampler with a better guess of where the final posterior area is. 

x = getSample(out, start = 1500)
# because of the low sample size, I don't trust the correlations, will thus only look at means and the range, and use this as new values for the sampler

meansPost = apply(x, 2, mean)
sdPost = apply(x, 2, sd)
rangePost = apply(x, 2, range)

newZ = matrix(runif(1500, rangePost[1,], rangePost[2,]), ncol = 3, byrow = T)

settings = list( Z = newZ, startValue =  x[(nrow(x)-2):nrow(x), ])
out <- runMCMC(bayesianSetup = bayesianSetup,  sampler = "DEzs", settings = settings )
plot(out)
plot(out, start = 1500)

# Re-starting the Metropolis - not sure if I'm doing something wrong, but didn't work very well

settings = list( proposalGenerator = createProposalGenerator(rangePost[2,] - rangePost[1,]) , startValue =  meansPost, optimize = F)
out <- runMCMC(bayesianSetup = bayesianSetup,  sampler = "Metropolis", settings = settings )
plot(out)
plot(out, start = 1500)








