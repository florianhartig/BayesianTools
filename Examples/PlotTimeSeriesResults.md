# Demonstration of plotTimeSeriesResults
Tankred Ott  
3 Jul 2017  

# ATTENTION

This doesn't really make sense as the model used is not a time series


```r
set.seed(123)
library(BayesianTools)
```


# Creation of test case



```r
a <- 5
b <- 10
sigma <- 10
rsigma = 30
group = rep(1:11, each = 5)

x <- -27:27
y <- a * x + b + rnorm(55,0,sd = sigma)
plot(x,y, col = group, pch = 3)
```

![](PlotTimeSeriesResults_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

# Fitting the model with BayesianTools 


```r
likelihood <- function(par){
  
  ax = par[1]
  bx = par[2]
  sigmax <- par[3]
  
  llObservation = sum(dnorm(ax * x + bx - y , sd = sigmax, log = T))
  return(llObservation)
  
}
```



```r
library(BayesianTools)

setup <- createBayesianSetup(likelihood = likelihood, lower = c(-20,-20,0), upper = c(20,20,50))
settings <- list(iterations = 10000, message=F)

res <- runMCMC(bayesianSetup = setup, settings = settings)
summary(res)
```

```
## # # # # # # # # # # # # # # # # # # # # # # # # # 
## ## MCMC chain summary ## 
## # # # # # # # # # # # # # # # # # # # # # # # # # 
##  
## # MCMC sampler:  DEzs 
## # Nr. Chains:  3 
## # Iterations per chain:  3334 
## # Rejection rate:  0.761 
## # Effective sample size:  1012 
## # Runtime:  1  sec. 
##  
## # Parameter       MAP      2.5%    median   97.5% 
## #  par 1 :        5.011    4.838    5.011    5.168 
## #  par 2 :       10.602    6.448   10.474   12.946 
## #  par 3 :        8.865    7.585    9.252   12.027 
## 
## ## DIC:  414.021 
## ## Convergence 
##  Gelman Rubin multivariate psrf:   
##  
## ## Correlations 
##        par 1  par 2  par 3
## par 1  1.000 -0.087 -0.045
## par 2 -0.087  1.000 -0.721
## par 3 -0.045 -0.721  1.000
```

```r
marginalPlot(res)
```

![](PlotTimeSeriesResults_files/figure-html/unnamed-chunk-4-1.png)<!-- -->



```r
createPredictions <- function (par) {
  ax = par[1]
  bx = par[2]
  sigmax <- par[3]
  predicted <- ax * x  + bx
  return(predicted)
}

createError <- function(mean, par){
  return(rnorm(length(mean), mean = mean, sd = par[3]))
}
```

# Plot with DHARMa's simulated residuals


```r
plotTimeSeriesResults(sampler = res,
                      model = createPredictions,
                      observed = y,
                      error = createError,
                      start = 2500,
                      plotResiduals = TRUE)
```

![](PlotTimeSeriesResults_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

# Plot without DHARMa's simulated residuals


```r
plotTimeSeriesResults(sampler = res,
                      model = createPredictions,
                      observed = y,
                      error = createError,
                      start = 2500,
                      plotResiduals = FALSE)
```

![](PlotTimeSeriesResults_files/figure-html/unnamed-chunk-7-1.png)<!-- -->


