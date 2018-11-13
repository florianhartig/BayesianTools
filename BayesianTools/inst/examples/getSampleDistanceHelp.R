library(mvtnorm)
library(BayesianTools)

sigma <- matrix(c(4,2,2,3), ncol=2)
X <- rmvnorm(n=2500, mean=c(1,2), sigma=sigma)
X = rbind(X,X)
Y <- rmvnorm(n=2500, mean=c(5,2), sigma=sigma)
Y = rbind(Y,X)

FNN::KL.dist(X, Y, k=10)

getSampleDistance(X,Y, type = "KL")
getSampleDistance(X,Y, type = "BH")
getSampleDistance(X,Y, type = "D")

# the distance functions are independent of the scale of the parameters, see e.g.

scale = 2
getSampleDistance(data.frame(rnorm(1000), rnorm(1000, sd = scale)), 
                  data.frame(rnorm(1000), rnorm(1000, mean = scale, sd = scale)), type = "BH")
getSampleDistance(data.frame(rnorm(1000), rnorm(1000, sd = scale)), 
                  data.frame(rnorm(1000), rnorm(1000, mean = scale, sd = scale)), type = "D")
