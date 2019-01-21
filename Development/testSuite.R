rosenbrock <- function (x) {
  if(!is.matrix(x)) x <- matrix(x, ncol = length(x))
  s <- rep(0, nrow(x))
  for(i in 1:(ncol(x) - 1)) {
    s <- s + (100 * (x[,i+1] - x[,i]^2) + (x[,i] - 1)^2)
  }
  return(s)
}

sphere <- function(x) {
  if(!is.matrix(x)) x <- matrix(x, ncol = length(x))
  return(apply(x^2, 1, sum))
}

m <- matrix(c(1,0.99,1.01,1,1,1.01), ncol = 2)

m
rosenbrock(m)

m2 <- matrix(c(1,2,3,1,2,3), ncol = 2)
apply(m2^2, 1, sum)

sphere(c(0,0,0))

# x1 <- seq(-2,2, length.out = 100)
# x2 <- x1
x1 <- runif(100, -3, 3)
x2 <- runif(100, -3, 3)
# x <- as.matrix(cbind(x1,x2))
# y <- matrix(0, nrow = length(x1), ncol = length(x2))
# for(i in 1:length(x1)) {
#   for(j in 1:length(x2)) {
#     y[i,j] <- rosenbrock(c(x1[i], x2[j]))
#   }
# }
y <- rosenbrock(as.matrix(cbind(x1, x2)))

contour(x1, x2, y)

ll <- function(x) {
  # pred <- dnorm(rosenbrock(x), 0, 0.05, log = T)
  pred <- 0 - rosenbrock(x)^2
  return(sum(pred))
}

ll2 <- function(x) {
  # pred <- dnorm(rosenbrock(x), 0, 0.05, log = T)
  pred <- 0 - sphere(x)^2
  return(sum(pred))
}

setup <- createBayesianSetup(ll2, lower = c(-10,-10), upper = c(10, 10))
settings <- list(iterations = 100000)
out <- runMCMC(setup, sampler = "DEzs", settings = settings)
plot(out)


