% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmcDE.R
\name{DE}
\alias{DE}
\title{Differential-Evolution MCMC}
\usage{
DE(
  bayesianSetup,
  settings = list(startValue = NULL, iterations = 10000, f = -2.38, burnin = 0, thin = 1,
    eps = 0, consoleUpdates = 100, blockUpdate = list("none", k = NULL, h = NULL, pSel =
    NULL, pGroup = NULL, groupStart = 1000, groupIntervall = 1000), currentChain = 1,
    message = TRUE)
)
}
\arguments{
\item{bayesianSetup}{a BayesianSetup with the posterior density function to be sampled from}

\item{settings}{list with parameter settings}

\item{startValue}{(optional) eiter a matrix with start population, a number to define the number of chains that are run or a function that samples a starting population.}

\item{iterations}{number of function evaluations.}

\item{burnin}{number of iterations treated as burn-in. These iterations are not recorded in the chain.}

\item{thin}{thinning parameter. Determines the interval in which values are recorded.}

\item{f}{scaling factor gamma}

\item{eps}{small number to avoid singularity}

\item{blockUpdate}{list determining whether parameters should be updated in blocks. For possible settings see Details.}

\item{message}{logical determines whether the sampler's progress should be printed}
}
\description{
Differential-Evolution MCMC
}
\details{
For blockUpdate the first element in the list determines the type of blocking.
Possible choices are
\itemize{
\item{"none"}{ (default), no blocking of parameters}
\item{"correlation"} { blocking based on correlation of parameters. Using h or k (see below)}
\item{"random"} { random blocking. Using k (see below)}
\item{"user"} { user defined groups. Using groups (see below)}
}
Further seven parameters can be specified. "k" determnined the number of groups, "h" the strength
of the correlation used to group parameter and "groups" is used for user defined groups.
"groups" is a vector containing the group number for each parameter. E.g. for three parameters
with the first two in one group, "groups" would be c(1,1,2).
Further pSel and pGroup can be used to influence the choice of groups. In the sampling process
a number of groups is randomly drawn and updated. pSel is a vector containing relative probabilities
for an update of the respective number of groups. E.g. for always updating only one group pSel = 1.
For updating one or two groups with the same probability pSel = c(1,1). By default all numbers
have the same probability.
The same principle is used in pGroup. Here the user can influence the probability of each group
to be updated. By default all groups have the same probability.
Finally "groupStart" defines the starting point of the groupUpdate and "groupIntervall" the intervall
in which the groups are evaluated.
}
\examples{
library(BayesianTools)

ll <- generateTestDensityMultiNormal(sigma = "no correlation")
bayesianSetup <- createBayesianSetup(likelihood = ll, 
                                     lower = rep(-10, 3), 
                                     upper = rep(10, 3))

settings = list(iterations = 200)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
summary(out)

# DE family samplers are population MCMCs that run a number of internal chains
# in parallel. Here examples how to change the internal chains
# note that internal chains can be executedi n parallel
settings = list(startValue = 4, iterations = 200)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
summary(out)

# Modify the start values of the internal chains (note that this is a matrix
# of dim nChain * nPar)
settings = list(startValue = matrix(rnorm(12), nrow = 4, ncol = 3), 
                iterations = 200)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
summary(out)

# In the DE sampler family with Z matrix, the previous chains are written in 
# a common matrix, from which proposals are generated. Per default this matrix
# is started with samples from the prior, but we can change this. Often useful
# to improve sampler convergence, 
# see  https://github.com/florianhartig/BayesianTools/issues/79
settings = list(startValue = matrix(rnorm(12), nrow = 4, ncol = 3),
                Z = matrix(rnorm(300), nrow = 100, ncol = 3),
                iterations = 200)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
summary(out)


}
\references{
Braak, Cajo JF Ter. "A Markov Chain Monte Carlo version of the genetic algorithm Differential Evolution: easy Bayesian computing for real parameter spaces." Statistics and Computing 16.3 (2006): 239-249.
}
\seealso{
\code{\link{DEzs}}
}
\author{
Francesco Minunno and Stefan Paul
}
