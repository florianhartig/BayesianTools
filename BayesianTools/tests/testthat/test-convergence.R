# context("Test convergence of all samplers")
# 
# skip_on_cran()
# 
# 
# set.seed(1)
# library(BayesianTools)
# 
# 
# testfun <- list(testDensityNormal, testDensityMultiNormal, testDensityBanana)
# nrPar <- c(1,3,2)
# 
# samp <- getPossibleSamplerTypes()
# 
# 
# reference <- list(rnorm(1000), 
#                   mvtnorm::rmvnorm(1000, mean = rep(0, 3), sigma =  emulator::corr.matrix(cbind(c(0.2, 0.3, 0.3503)),scales=1) ),
#                   mvtnorm::rmvnorm(1000, mean = rep(0, 2), sigma = matrix(nrow = 2, data = c(1, 0.9, 0.9, 1))) )
# 
# 
# iter = 1000000
# settings <- list(iterations = iter, message = FALSE)
# 
# for(i in 1:length(testfun)){
#   
#   fun <- testfun[[i]]
#   setup <- createBayesianSetup(likelihood = fun, lower = rep(-10, nrPar[i]), 
#                                upper = rep(10, nrPar[i]))
#   
#   for(k in 1:(length(samp$BTname)-1)){
#     if(nrPar[i] != 1 | samp$univariatePossible[k] == T){ 
#     out <- runMCMC(setup, sampler = samp$BTname[k], settings = settings)
#     
#       
#     sample <- getSample(out, numSamples = 1000, start = 1000)
#     
#     if(nrPar[i] == 1){
#       ks <- ks.test(as.vector(sample), reference[[i]])$p.value
#       expect_true(ks>0.05)
#       
#     }else{
#     for(z in 1:ncol(sample)){
#      ks <- ks.test(sample[,z], reference[[i]][,z])$p.value
#      expect_true(ks>0.01) 
#     }
#     }
#     
#     #ksVal <-  ks.boot(getSample(out, numSamples = 10000), rnorm(10000))$ks.boot.pvalue 
#     } # if
#     
#     }
# 
#   }
#   
#  
#   
# 
# 
