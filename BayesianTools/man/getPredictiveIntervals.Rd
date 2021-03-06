% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmcConfidence.R
\name{getPredictiveIntervals}
\alias{getPredictiveIntervals}
\title{Calculates Bayesian credible (confidence) and predictive intervals based on parameter sample}
\usage{
getPredictiveIntervals(
  parMatrix,
  model,
  numSamples = 1000,
  quantiles = c(0.025, 0.975),
  error = NULL
)
}
\arguments{
\item{parMatrix}{matrix of parameter values}

\item{model}{model / function to calculate predictions. Outcome should be a vector}

\item{numSamples}{number of samples to be drawn}

\item{quantiles}{quantiles to calculate}

\item{error}{function with signature f(mean, par) that generates error expectations from mean model predictions. Par is a vector from the matrix with the parameter samples (full length). f needs to know which of these parameters are parameters of the error function. If supplied, will calculate also predictive intervals additional to credible intervals}
}
\description{
Calculates Bayesian credible (confidence) and predictive intervals based on parameter sample
}
\details{
If numSamples is greater than the number of rows in parMatrix, or NULL, or FALSE, or less than 1 all samples in parMatrix will be used.
}
\seealso{
\code{\link{getPredictiveDistribution}} \cr
\code{\link{getCredibleIntervals}} \cr
}
\author{
Florian Hartig
}
