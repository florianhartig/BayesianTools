% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmcRun.R
\name{runMCMC}
\alias{runMCMC}
\title{Main wrapper function to start MCMCs, particle MCMCs and SMCs}
\usage{
runMCMC(bayesianSetup, sampler = "DEzs", settings = NULL)
}
\arguments{
\item{bayesianSetup}{either a BayesianSetup (see \code{\link{createBayesianSetup}}) or a BayesianOutput created by runMCMC. The latter allows to continue a previous MCMC run. See details for how to restart a sampler.}

\item{sampler}{sampling algorithm to be run. Default is DEzs. Options are "Metropolis", "AM", "DR", "DRAM", "DE", "DEzs", "DREAM", "DREAMzs", "SMC". For details see the help of the individual functions.}

\item{settings}{list with settings for each sampler (see help of sampler for details). If a setting is not provided, defaults (see \code{\link{applySettingsDefault}}) will be used. You can see the default values by running \code{\link{applySettingsDefault}} with the respective sampler name, or in the help of the samplers.}
}
\value{
The function returns an object of class mcmcSampler (if one chain is run) or mcmcSamplerList. Both have the superclass bayesianOutput. It is possible to extract the samples as a coda object or matrix with \code{\link{getSample}}.
It is also possible to summarize the posterior as a new prior via \code{\link{createPriorDensity}}.
}
\description{
Main wrapper function to start MCMCs, particle MCMCs and SMCs
}
\details{
The runMCMC function can be started with either one of a) an object of class BayesianSetup with prior and likelihood function (recommended, see \code{\link{createBayesianSetup}}), b) a log posterior or other target function, or c) an object of class BayesianOutput created by runMCMC. The latter allows to continue a previous MCMC run. If a bayesianSetup is provided, check if appropriate parallization options are used - many samplers can make use of parallelization if this option is activated when the class is created.

For details about the different MCMC samplers, make sure you have read the Vignette (run vignette("BayesianTools", package="BayesianTools"). Also, see \code{\link{Metropolis}} for Metropolis based samplers, \code{\link{DE}} and \code{\link{DEzs}} for standard differential evolution samplers, \code{\link{DREAM}} and \code{\link{DREAMzs}} for DREAM sampler, \code{\link{Twalk}} for the Twalk sampler, and \code{\link{smcSampler}} for rejection and Sequential Monte Carlo sampling.\cr

The samplers "AM", "DR", and "DRAM" are special cases of the "Metropolis" sampler and are shortcuts for predefined settings ("AM": adapt=TRUE; "DR": DRlevels=2; "DRAM": adapt=True, DRlevels=2).

The settings list allows to change the settings for the MCMC samplers and some other options. For the MCMC sampler settings, see their help files. Global options that apply for all MCMC samplers are: iterations (number of MCMC iterations), and nrChains (number of chains to run). Note that running several chains is not done in parallel, so if time is an issue it will be better to run the MCMCs individually and then combine them via \code{\link{createMcmcSamplerList}} into one joint object.

Startvalues: all samplers allow to provide explicit startvalues. If startvalues are not provided, they are sampled from the prior. Usually, this is a good choice, so don't feel compelled to provide startvalues.

Note that DE and DREAM variants as well as SMC and T-walk require a population to start, which should be provided as a matrix. Default (NULL) sets the population size for DE to 3 x dimensions of parameters, for DREAM to 2 x dimensions of parameters and for DEzs and DREAMzs to three, sampled from the prior. Note also that the zs variants of DE and DREAM require two populations, the current population and the z matrix (a kind of memory) - if you want to set both, provide a list with startvalue$X and startvalue$Z.

setting startValue for sampling with nrChains > 1 : if you want to provide different start values for the different chains, provide them as a list
}
\examples{
## Generate a test likelihood function. 
ll <- generateTestDensityMultiNormal(sigma = "no correlation")

## Create a BayesianSetup object from the likelihood 
## is the recommended way of using the runMCMC() function.
bayesianSetup <- createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))

## Finally we can run the sampler and have a look
settings = list(iterations = 1000, adapt = FALSE)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)

## out is of class bayesianOutput. There are various standard functions 
# implemented for this output

plot(out)
correlationPlot(out)
marginalPlot(out)
summary(out)

## additionally, you can return the sample as a coda object, and make use of the coda functions
# for plotting and analysis

codaObject = getSample(out, start = 500, coda = TRUE)
}
\seealso{
\code{\link{createBayesianSetup}}
}
\author{
Florian Hartig
}
