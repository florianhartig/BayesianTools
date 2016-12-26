Sys.setenv("R_TESTS" = "")

library(BayesianTools)
library(testthat)

test_check("BayesianTools")
