context("Test convergence of all samplers")



set.seed(1)
library(BayesianTools)


testfun <- testDensityNormal # TODO replace 1d normal 
nrPar <- c(1,2)

samp <- c(rep("Metropolis", 4), "DE", "DEzs", "DREAM", "DREAMzs", "Twalk")



iter = 100000

settings <- list(list(iterations = iter, adapt = F, message = FALSE),
                 list(iterations = iter, adapt = T, message = FALSE),
                 list(iterations = iter, adapt = F, DRlevels = 2, message = FALSE),
                 list(iterations = iter, adapt = T, DRlevels = 2, message = FALSE),
                 list(iterations = iter, message = FALSE),
                 list(iterations = iter, message = FALSE),
                 list(iterations = iter, message = FALSE),
                 list(iterations = iter, message = FALSE),
                 list(iterations = iter, message = FALSE)
                 
                 )


test_that("All samplers sample correctly",{
      skip_on_cran()
      skip_on_travis()
          
for(i in 1:2){

  fun <- testfun
  setup <- createBayesianSetup(likelihood = fun, lower = rep(-10, nrPar[i]),
                               upper = rep(10, nrPar[i]))

  for(k in 1:(length(samp))){
    out <- runMCMC(setup, sampler = samp[k], settings = settings[[k]])

    sample <- suppressWarnings(getSample(out, numSamples = 1000, start = 1000))

    if(nrPar[i] == 1){
      ks <- suppressWarnings(ks.test(as.vector(sample), rnorm(1000))$p.value)
      expect_true(ks>0.01)

    }else{
    for(z in 1:ncol(sample)){
     ks <- suppressWarnings(ks.test(sample[,z], rnorm(1000))$p.value)
     expect_true(ks>0.01)
    }
    }

    #ksVal <-  ks.boot(getSample(out, numSamples = 10000), rnorm(10000))$ks.boot.pvalue

    }

  }}
)




