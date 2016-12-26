# Calibrating a forest model - a practical tutorial
Florian Hartig  
23 Oct 2015  

**Acknowledgements**: Many ideas in this script originate from tutorials / discussions of the COST FP1304 spring school on model calibration. Particularly noteworthy contributions are Francesco Minunno for providing a first tutorial on calibrating the Preles model, and Björn Reineking for the idea of displaying confidence / prediction bands for the fitted models. 

**Used packages:** In this example, we will use the following packages

* BayesianTools, containing MCMC sampler and plotting functions (obtain this via the PROFOUND github repo / contact Florian Hartig)
* Sensitivity, an R package for sensitivity analysis
* DEoptim, and R package for optimization


```r
library(BayesianTools)
library(sensitivity)
library(DEoptim)
```

```
## 
## DEoptim package
## Differential Evolution algorithm in R
## Authors: D. Ardia, K. Mullen, B. Peterson and J. Ulrich
```

```r
set.seed(123)
```

## A forest model and data 

We will calibrate the Preles mode, available in the package Rpreles maintained by Mikko Peltoniemi (available via the PROFOUND github repo). For a more detailed description, see the package.

If you want to run this script with your own forest model, you should place it here. 

* A tutorial on how to make your model callable from R is available [here](https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/09-BayesAndProcessBasedModels/LinkingAModelToR.md)
* TG13 TODO - extend this tutorial


```r
library(Rpreles)
```

Some flux data provided by Francesco Minunno


```r
load("/Users/Florian/Home/Teaching/Vorlesungen/USM/USM/Code/PRELES/Sensitivity/Boreal_sites.rdata")
load('/Users/Florian/Home/Teaching/Vorlesungen/USM/USM/Code/PRELES/Sensitivity/par.rdata')
```

I choose the parameters that we want to put under calibration


```r
parind <- c(5:8) #Indeces for PRELES parameters
nparModel <- length(parind)
print(par$name[parind])
```

```
## [1] "beta"   "tau"    "X[0]"   "S[max]"
```

Define a run Model function that will run the model with these parameters and return an output that we want to use for calibration, in this case GPP.



```r
runModel <- function(x){
  paras<-par$def
  paras[parind]<-x[1:nparModel]
  predicted<-PRELES(PAR=s1$PAR, TAir=s1$TAir, VPD=s1$VPD,Precip=s1$Precip, CO2=s1$CO2, fAPAR=s1$fAPAR, p=paras)$GPP
  return(predicted)
}
```


Define observations that we want to calibrate to


```r
observed = s1$GPPobs
```


## Build likelihood function / calibration target

Need a calibration target = distance between model and data.

For statistical calibration, this is always the likelihood = p(data| model, parameters). The likelihood is calculated based on probability distribution that the user has to specify and that will typically be also optimized during calibration. For an explanation of this approach, see Hartig, F.; Dyke, J.; Hickler, T.; Higgins, S. I.; O'Hara, R. B.; Scheiter, S. & Huth, A. (2012) Connecting dynamic vegetation models to data - an inverse perspective. J. Biogeogr., 39, 2240-2252.

For this example, I use a normal likelihood, but we will see later that this is not a perfect choice. More flexible / sophisticated likelihood functions are provided in the BayesianTools package.

* TODO TG13 - tutorial on the choice of the likelihood

#### Likelihood definition


```r
likelihood <- function(x){
  predicted = runModel(x[1:nparModel])
  ll <- likelihoodIidNormal(predicted = predicted, observed = observed, x[nparModel+1])
  return(ll)
}

x <-c(par$def[parind],1) 
likelihood(x)
```

```
## [1] -834.0919
```

## Sensitivity analysis

Sensitivity analysis (SA) allows us to see which parameters affect the ouptut most strongly. If you don't know your model well, it makes sense to run a sensitivity analysis in advance of calibration / MCMC to see which parameters have a strong influence on the output. We could apply it on any output, but as we are interested in the likelihood, it makes sense to apply it on the likelihood directly.

There are a large number of local and global SA methods. Here, I use the Morris screening, which gives a global estimate of importance (sensitivity) and nonlinearities (sigma) of each parameter. For details on (global) SA, see Saltelli, A. (2004) Sensitivity analysis in practice: a guide to assessing scientific models. John Wiley & Sons Inc, 

* TODO TG13 - tutorial on SA?

#### Running the SA


```r
# This function of BayesianTools generates the form of the likelihood that is needed for the sensitivity package
parLL <- generateParallelExecuter(likelihood)

mins <- c(par$min[parind] , 0.001)
maxs <- c(par$max[parind], 20)
parNames <- c(par$name[parind] , "sd")

system.time(morrisOut <- morris(model = parLL, factors = parNames, r = 500, design = list(type = "oat", levels = 5, grid.jump = 3), binf = mins, bsup = maxs, scale = TRUE)) 
```

```
## Warning in morris(model = parLL, factors = parNames, r = 500, design =
## list(type = "oat", : keeping 497 repetitions out of 500
```

```
##    user  system elapsed 
##   1.739   0.043   1.807
```

```r
plot(morrisOut)
```

![](ForestModel_files/figure-html/unnamed-chunk-8-1.png) 

It seems the sd of the likelihood as well as beta and X0 are the most important parameters.

## Optimization

Running a global optimization, trying to find the point of highest likelihood, and plotting the result.


```r
out<- DEoptim(function(x)-likelihood(x), lower = c(par$min[parind],0), upper =  c(par$max[parind],20) )
paras<-par$def
paras[parind]<-out$optim$bestmem[1:nparModel]

pred = PRELES(PAR=s1$PAR, TAir=s1$TAir, VPD=s1$VPD, Precip=s1$Precip, CO2=s1$CO2, fAPAR=s1$fAPAR, p=paras)$GPP

plotTimeSeries(observed = s1$GPPobs, predicted = pred)
```

![](ForestModel_files/figure-html/unnamed-chunk-9-1.png) 

TODO TG13 - tutorial on optimization?

## Bayesian calibration / MCMC

Running an MCMC for Bayesian calibration. The current examples assumes flat priors, but you can provide a different prior distribution to the MCMC. 

* If you want to know more about Bayesian calibration in general the [Learning Bayes Website](http://florianhartig.github.io/LearningBayes/) TODO TG13 - tuturial on Bayesian Analysis.
  
* If you want to know more about prior choice see [here](https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/01-Principles/Priors.md) TODO TG13 - tutorial on prior choice.

* If you want to know more about MCMC algorithms see the [tuturial MCMC](https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/02-Samplers/MCMC/Metropolis.md). The [tutorial rejection sampling](https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/02-Samplers/Rejection/ExampleRejectionSampler.md) and [tutorial SMC](https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/02-Samplers/SMC/SMC.md) demonstrate some alternatives to MCMC that are easier parallelizable. TODO TG13 - tutorial on MCMC / sampling in general.


#### Running the MCMC




```r
proposalGenerator = createProposalGenerator((maxs - mins)/100)
sampler <-mcmcSampler(likelihood = likelihood, startvalue = out$optim$bestmem, proposalGenerator = proposalGenerator, optimize = F)
```

```
## Loading required package: compiler
```

```r
sampler<- getSamples(sampler, 2000)
```

```
## Done.
```

```r
sampler<- getSamplesAdaptive(sampler, 5, 2000)
```

```
## Done.
## Done.
## Done.
## Done.
## Done.
```


#### MCMC results

* see [tutorial posterior interpretation](https://github.com/florianhartig/LearningBayes/blob/master/CommentedCode/01-Principles/Posterior.md)
* TODO TG13 - tutorial on posterior interpretation



```r
summary(sampler)
```

```
## 
## Iterations = 1:12001
## Thinning interval = 1 
## Number of chains = 1 
## Sample size per chain = 12001 
## 
## 1. Empirical mean and standard deviation for each variable,
##    plus standard error of the mean:
## 
##         Mean       SD  Naive SE Time-series SE
## 1     0.6997 0.004627 4.223e-05      0.0002509
## 2     8.9561 0.533872 4.873e-03      0.0331458
## 3    -2.8261 0.300364 2.742e-03      0.0240074
## 4    17.0606 0.422124 3.853e-03      0.0318971
## 5     0.6159 0.015010 1.370e-04      0.0007258
## LL -679.8794 1.710210 1.561e-02      0.1445382
## LP -679.8794 1.710210 1.561e-02      0.1445382
## 
## 2. Quantiles for each variable:
## 
##         2.5%       25%       50%       75%     97.5%
## 1     0.6914    0.6965    0.6993    0.7026    0.7094
## 2     8.0034    8.5575    8.9961    9.2318   10.1416
## 3    -3.4526   -2.9766   -2.8764   -2.6442   -2.2523
## 4    16.2817   16.8167   17.1145   17.2847   17.9537
## 5     0.5876    0.6085    0.6123    0.6254    0.6476
## LL -684.0195 -680.9041 -679.7048 -678.6784 -677.7282
## LP -684.0195 -680.9041 -679.7048 -678.6784 -677.7282
```

```r
plot(sampler)
```

![](ForestModel_files/figure-html/unnamed-chunk-11-1.png) ![](ForestModel_files/figure-html/unnamed-chunk-11-2.png) 

```r
correlationPlot(sampler)
```

![](ForestModel_files/figure-html/unnamed-chunk-11-3.png) 

```r
marginalPlot(sampler, bounds = cbind(mins, maxs))
```

```
## Loading required package: vioplot
## Loading required package: sm
## Package 'sm', version 2.2-5.4: type help(sm) for summary information
## 
## Attaching package: 'sm'
## 
## The following object is masked from 'package:MASS':
## 
##     muscle
```

![](ForestModel_files/figure-html/unnamed-chunk-11-4.png) 

```r
errorFunction <- function(mean, par) rnorm(length(mean), mean = mean, sd = par[nparModel+1])

plotTimeSeriesAuto(mcmcSampler = sampler, model = runModel, observed = observed, error = errorFunction)
```

![](ForestModel_files/figure-html/unnamed-chunk-11-5.png) 





