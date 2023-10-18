

data = rnorm(20)

# create standardized likelihood density for Gaussian likelihood function
# fitting parameters mean and standard deviation

likelihood <- function(x) {
  predicted <- rep(x[1], length(data))
  LL = likelihoodIidNormal(predicted, data, x[2])
  return(LL)
}


# showing posterior response surface for parameter 1
# note that this is plotted at 
bayesianSetup$prior$best
plotSensitivity(bayesianSetup, selection = 1)

plotSensitivity(bayesianSetup)
plotSensitivity(bayesianSetup, equalScale = F)


