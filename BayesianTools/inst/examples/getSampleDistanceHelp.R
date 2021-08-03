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


# To get an idea about how the distance functions approach zero for a perfect sampler
# use the following code where we sample 2 distributions directly from their 
# rnorm function 

getD <- function(n){
  x1 = matrix(rnorm(n), ncol = 10)
  x2 = matrix(rnorm(n), ncol = 10)
  getSampleDistance(x1, x2, type = "D" )  
}

test = seq(1000, 30000, 1000)
plot(test, sapply(test, getD), ylim = c(0,0.4))
abline(h=0, col = "red")
