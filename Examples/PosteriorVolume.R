library(BayesianTools)

 

#x1 is technically unidentifiable, only change sensitivity to x2
#x2 has maximum at zero, regardless of x1
ll <- function(x){
  dens = (- x[2]^2 / exp(x[1])) 
}

x <- seq(-3,3,length.out=100)
y <- seq(-3,3,length.out=100)
z <- outer(x, y, Vectorize(function(x,y) ll(c(x,y))))

image(x, y, log(-z), xlab = "par1", ylab = "par2")


# fitting this likelihood with BT and flat priors
bayesianSetup <- createBayesianSetup(likelihood = ll, lower = rep(-10, 2), upper = rep(10, 2))
out <- runMCMC(bayesianSetup = bayesianSetup)

marginalPlot(out, histogram = F, singlePanel = T)
correlationPlot(out)
MAP(out)
