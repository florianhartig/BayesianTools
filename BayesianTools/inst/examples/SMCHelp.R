# creating a simple test BayesianSetup
bayesianSetup <- createBayesianSetup(likelihood = generateTestDensityMultiNormal(), 
                                     lower = rep(-10, 3),
                                     upper = rep(10, 3), parallel = "external")

# to directly call the smcSampler, use 
SMCout <- smcSampler(bayesianSetup, initialParticles = 100, iterations = 3)

# plotting results
marginalPlot(SMCout, type = "d", singlePanel = T)

# for options, see 
?smcSampler

# as for all MCMCs, you can also run the SMC via the general runMCMC function
SMCout<-runMCMC(bayesianSetup, sampler = "SMC", 
                settings = list(nrChains =2, initialParticles = 100, iterations = 3))
plot(SMCout)

