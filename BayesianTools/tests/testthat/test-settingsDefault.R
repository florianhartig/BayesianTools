context("Test Settings Default")

skip_on_cran()

set.seed(1)
library(BayesianTools)

test_that("Default works in principle",{
  skip_on_cran()
  
  settings <- list(iterations = 20000, adapt = T, DRlevels = 2, optimize = T, burnin=1000, adaptationInterval=10)
  applySettingsDefault(settings = settings, sampler = "Metropolis", check = FALSE)
  applySettingsDefault(settings = settings, sampler = "Metropolis", check = TRUE)

  }
)



test_that("Wrong inputs are caught",{
  skip_on_cran()

  #####################
  # Adaptation before burnin 
  
  settings <- list(iterations = 20000, adapt = T, DRlevels = 2, optimize = T, burnin=7000, adaptationInterval=10)  
  
  expect_error({
    # default adaptation is 3000, so the burnin is larger than that, should throw error
    applySettingsDefault(settings = settings, sampler = "Metropolis")
  })
  
  # should not throw an error if the sampler doesn't adapt 
  applySettingsDefault(settings = settings, sampler = "DE")    

  
  #####################
  # Burnin larger than iterations
  
  settings <- list(iterations = 2000, adapt = T, DRlevels = 2, optimize = T, burnin=2500, adaptationInterval=10)  
  
  expect_error( applySettingsDefault(settings = settings, sampler = "Metropolis") )
  expect_error( applySettingsDefault(settings = settings, sampler = "DE") )  
  


  
  
  
  
  
}
)



