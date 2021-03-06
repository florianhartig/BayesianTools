% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codaFunctions.R
\name{combineChains}
\alias{combineChains}
\title{Function to combine chains}
\usage{
combineChains(x, merge = T)
}
\arguments{
\item{x}{a list of MCMC chains}

\item{merge}{logical determines whether chains should be merged}
}
\value{
combined chains
}
\description{
Function to combine chains
}
\note{
to combine several chains to a single McmcSamplerList, see \code{\link{createMcmcSamplerList}}
}
\keyword{internal}
