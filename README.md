[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/BayesianTools)](https://cran.r-project.org/package=BayesianTools)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.1.2-6666ff.svg)](https://cran.r-project.org/)
[![](https://cranlogs.r-pkg.org/badges/BayesianTools)](https://cran.r-project.org/package=BayesianTools)

# BayesianTools

R package for performing Bayesian inference, including various MCMC and SMC sampling algorithms! A google scholar list of studies citing BT is [here](https://scholar.google.com/scholar?oi=bibs&hl=en&cites=15191862660990668802,4910024394983334104&as_sdt=5).

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

### Installing the development version from GitHub

If you want to install our development version from GitHub, use 

```{r}
devtools::install_github(repo = "florianhartig/BayesianTools", subdir = "BayesianTools", dependencies = T, build_vignettes = T)
```

New developments will be done in extra branches and will be tested before merging in the developtment branch, so the developmet version should usually be usable (consider it in a beta stage), while feature branches should be considered alpha. 


Status master development branch 

* GH Action [![R-CMD-check](https://github.com/florianhartig/BayesianTools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/florianhartig/BayesianTools/actions/workflows/R-CMD-check.yaml)

Windows users: the package contains c++ code, so if you compile yourself, you need [RTools](https://cran.r-project.org/bin/windows/Rtools/) installed. 

### Older releases

To install a specific (older) release or a branch, decide for the version number / branch that you want to install in [https://github.com/florianhartig/BayesianTools/releases](https://github.com/florianhartig/BayesianTools/releases) (version numbering corresponds to CRAN, but there may be smaller releases that were not pushed to CRAN) and run, e.g.  

```{r}
devtools::install_github(repo = "florianhartig/BayesianTools", subdir = "BayesianTools", ref = "v0.0.10", dependencies = T, build_vignettes = T)
```
with v0.0.10 replaced by the appropriate version number or branch. 

## Getting help

We highly welcome questions by users, so don't be shy - any questions, even if it feels "stupid", helps us to understand how we can improve the interface, documentation, or code of the package. 

If you want to ask a question or report a bug, the most convenient way for us would be to provide a [reproducible example](http://stackoverflow.com/questions/5963269/how-to-make-a-great-r-reproducible-example) via the GitHub [issues](https://github.com/florianhartig/BayesianTools/issues)









