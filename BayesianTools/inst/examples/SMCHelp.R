## Example for the use of SMC 
# First we need a bayesianSetup - SMC makes most sense if we can  for demonstration,
# we'll write a function that puts out the number of model calls


bayesianSetup <- createBayesianSetup(likelihood = generateTestDensityMultiNormal(), 
                                     lower = rep(-10, 3),
                                     upper = rep(10, 3), parallel = "external")

SMCout <- smcSampler(bayesianSetup)
marginalPlot(SMCout, type = "d", singlePanel = T)



