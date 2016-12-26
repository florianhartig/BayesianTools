# BayesianTools

R package for performing Bayesian inference, including various MCMC and SMC sampling algorithms!

## Getting BayesianTools

### From CRAN 

BayesianTools is on [CRAN](https://cran.r-project.org/web/packages/BayesianTools/index.html). So, to install the latest major release, just run 

```{r}
install.packages("DHARMa")
```

To get an overview about its functionality once the package is installed, run

```{r}
library(BayesianTools)
?BayesianTools
vignette("DHARMa", package="BayesianTools")
```

### Development release 

If you want to install the current (development) version from this repository, run

```{r}
devtools::install_github(repo = "BayesianTools", username = "florianhartig", subdir = "BayesianTools", dependencies = T)
```
Below the status of the automatic Travis CI tests on the master branch (if this doesn load see [here](https://travis-ci.org/florianhartig/BayesianTools))

[![Build Status](https://travis-ci.org/florianhartig/BayesianTools.svg?branch=master)](https://travis-ci.org/florianhartig/BayesianTools)

### Older releases

To install a specific (older) release, decide for the version number that you want to install in [https://github.com/florianhartig/BayesianTools/releases](https://github.com/florianhartig/BayesianTools/releases) (version numbering corresponds to CRAN, but there may be smaller releases that were not pushed to CRAN) and run 

```{r}
devtools::install_github(repo = "BayesianTools", username = "florianhartig", subdir = "BayesianTools", ref = "v0.0.10")
```
with the appropriate version number. 








