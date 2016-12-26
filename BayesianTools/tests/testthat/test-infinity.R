# context("Test samplers on likelihoods with infinities")
# 
# 
# bayesianSetup = createBayesianSetup(likelihood = testDensityInfinity, lower = c(0, 0), upper = c(5, 5))
# 
# iter = 10000
# start = 500
# iterSMC = 400
# 
# test_that("sampler work correct for likelihoods with infinities", {
#   
#   skip_on_cran()
# 
# samp = getPossibleSamplerTypes()
# 
# for(i in 1:length(samp$BTname)){
#    # print(samp$BTname[i]) # Printing to console makes tests extremely slow
#   if(samp$univariatePossible[i] == T){
#     settings = list(iterations = iter, consoleUpdates = 1e+8)
#     if(samp$BTname[i] == "SMC") settings = list(iterations = iterSMC, consoleUpdates = 1e+8)
#     invisible(capture.output(suppressMessages(out <- runMCMC(bayesianSetup = setup, sampler = samp$BTname[i], settings = settings))))
#   
#   }
# #   plot(out)
# #   summary(out)
# #   marginalPlot(out)
# #   correlationPlot(out)
# #   DIC(out)
# #   marginalLikelihood(out)
#     
# #     x = getSample(out, numSamples  = 10000)
# #     y <- rnorm(10000)  ## TODO change
# #     for(z in 1:ncol(x)){
# #       
# #       #  ks <- ks.test(x[,z], pnorm)$p.value
# #       
# #       ks <- ks.boot(x[,z], y)$ks.boot.pvalue 
# #       
# #       # Test that distribution is not significally different from gaussian
# #       expect_true(ks>0.05)
# #       
# #     }
#     
# }
# }
# )
