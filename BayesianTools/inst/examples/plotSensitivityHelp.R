data = rnorm(20)

# create standardized likelihood density for Gaussian likelihood function
# fitting parameters mean and standard deviation
likelihood <- function(x) {
  predicted <- rep(x[1], length(data))
  LL = likelihoodIidNormal(predicted, data, x[2])
  return(LL)
}

setup = createBayesianSetup(likelihood, lower = c(-10, 0.01), upper = c(10,5))

# showing posterior response surface for parameter 1
# note that this is plotted at 
setup$prior$best
plotSensitivity(setup, selection = 2)

plotSensitivity(setup)
# when parameters are not on the same scale, we can use
plotSensitivity(setup, equalScale = FALSE)


