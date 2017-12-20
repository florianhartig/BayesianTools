[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/BayesianTools)](https://cran.r-project.org/package=BayesianTools)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.1.2-6666ff.svg)](https://cran.r-project.org/)

Status master [![Build Status](https://travis-ci.org/florianhartig/BayesianTools.svg?branch=master)](https://travis-ci.org/florianhartig/BayesianTools)

Status development [![Build Status](https://travis-ci.org/florianhartig/BayesianTools.svg?branch=development)](https://travis-ci.org/florianhartig/BayesianTools)

# BayesianTools

R package for performing Bayesian inference, including various MCMC and SMC sampling algorithms!

## Getting BayesianTools

BayesianTools is on CRAN (see [here](https://cran.r-project.org/web/packages/BayesianTools/index.html)). To install the latest CRAN release, type

```{r}
install.packages("BayesianTools")
```

To get an overview about its functionality once the package is installed, run

```{r}
library(BayesianTools)
?BayesianTools
vignette("BayesianTools", package="BayesianTools")
```

As for every R package, you can get the suggested citation via

```{r}
citation("BayesianTools")
```

### Development release 

If you want to install the current (development) version from this repository, run

```{r}
devtools::install_github(repo = "florianhartig/BayesianTools", subdir = "BayesianTools", dependencies = T, build_vignettes = T)
```
Below the status of the automatic Travis CI tests on the development version in the master branch (if this doesn load see [here](https://travis-ci.org/florianhartig/BayesianTools))

[![Build Status](https://travis-ci.org/florianhartig/BayesianTools.svg?branch=master)](https://travis-ci.org/florianhartig/BayesianTools)

Windows users: the package contains c++ code, so if you compile yourself, you need [RTools](https://cran.r-project.org/bin/windows/Rtools/) installed. 

### Older releases

To install a specific (older) release, decide for the version number that you want to install in [https://github.com/florianhartig/BayesianTools/releases](https://github.com/florianhartig/BayesianTools/releases) (version numbering corresponds to CRAN, but there may be smaller releases that were not pushed to CRAN) and run, e.g.  

```{r}
devtools::install_github(repo = "florianhartig/BayesianTools", subdir = "BayesianTools", ref = "v0.0.10", dependencies = T, build_vignettes = T)
```
with v0.0.10 replaced by the appropriate version number. 


## Getting help

We highly welcome questions by users, so don't be shy - any questions, even if it feels "stupid", helps us to understand how we can improve the interface, documentation, or code of the package. 

If you want to ask a question or report a bug, the most convenient way for us would be to provide a [reproducible example](http://stackoverflow.com/questions/5963269/how-to-make-a-great-r-reproducible-example) via the GitHub [issues](https://github.com/florianhartig/BayesianTools/issues)

## Acknowledgements

Work on this package was facilicated through meetings of [Cost Action FP 1304 Profound](http://www.cost.eu/COST_Actions/fps/FP1304). 







