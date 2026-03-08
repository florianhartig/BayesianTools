
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



png(filename = "man/figures/betaDensity.png", width = 290, height = 200)
plot(density(rbeta(10000000,10,3)), main = "Beta Density with a = 10, b = 3")
dev.off()
png(filename = "man/figures/normalDensity.png", width = 290, height = 200)
plot(density(rnorm(10000000,0,1)), main = "Normal Density with mean = 0, sd = 1")
dev.off()
png(filename = "man/figures/uniformDensity.png", width = 290, height = 200)
plot(density(runif(10000000)), main = "Uniform Density")
dev.off()



png(filename = "vignettes/betaDensity.png", width = 200, height = 180)
plot(density(rbeta(10000000,10,3)), main = "Beta Density with \n a = 10, b = 3", xlab = "n = 10000000")
dev.off()
png(filename = "vignettes/normalDensity.png", width = 200, height = 180)
plot(density(rnorm(10000000,0,1)), main = "Normal Density with \n mean = 0, sd = 1", xlab = "n = 10000000")
dev.off()
png(filename = "vignettes/uniformDensity.png", width = 200, height = 180)
plot(density(runif(10000000)), main = "Uniform Density", xlab = "n = 10000000")
dev.off()

