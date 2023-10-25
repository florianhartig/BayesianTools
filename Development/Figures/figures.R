##### plotSensitivity #####

data = rnorm(20)

# create standardized likelihood density for Gaussian likelihood function
# fitting parameters mean and standard deviation
likelihood <- function(x) {
  predicted <- rep(x[1], length(data))
  LL = likelihoodIidNormal(predicted, data, x[2])
  return(LL)
}

setup = createBayesianSetup(likelihood, lower = c(-10, 0.01), upper = c(10,5))

png(filename = "BayesianTools/man/figures/plotSensitivity-ScaleDontMatch.png", width = 400, height = 300)
plotSensitivity(setup)
dev.off()

png(filename = "BayesianTools/man/figures/plotSensitivity-ScaleFalse.png", width = 400, height = 300)
plotSensitivity(setup, equalScale = F)
dev.off()