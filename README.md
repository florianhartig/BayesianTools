# BayesianTools

R package for performing Bayesian inference, including various MCMC and SMC sampling algorithms!

## Getting BayesianTools

BayesianTools is on CRAN. To install the latest CRAN release, type

```{r}
install.packages("BayesianTools")
```

To get an overview about its functionality once the package is installed, run

```{r}
library(BayesianTools)
?BayesianTools
vignette("BayesianTools", package="BayesianTools")
```

### Development release 

If you want to install the current (development) version from this repository, run

```{r}
devtools::install_github(repo = "florianhartig/BayesianTools", subdir = "BayesianTools", dependencies = T, build_vignettes = T)
```
Below the status of the automatic Travis CI tests on the development version in the master branch (if this doesn load see [here](https://travis-ci.org/florianhartig/BayesianTools))

[![Build Status](https://travis-ci.org/florianhartig/BayesianTools.svg?branch=master)](https://travis-ci.org/florianhartig/BayesianTools)

### Older releases

To install a specific (older) release, decide for the version number that you want to install in [https://github.com/florianhartig/BayesianTools/releases](https://github.com/florianhartig/BayesianTools/releases) (version numbering corresponds to CRAN, but there may be smaller releases that were not pushed to CRAN) and run, e.g.  

```{r}
devtools::install_github(repo = "florianhartig/BayesianTools", subdir = "BayesianTools", ref = "v0.0.10", dependencies = T, build_vignettes = T)
```
with v0.0.10 replaced by the appropriate version number. 

## Acknowledgements

Work on this package was facilicated through meetings of [Cost Action FP 1304 Profound](http://www.cost.eu/COST_Actions/fps/FP1304). 






