% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-Numerics.R
\name{logSumExp}
\alias{logSumExp}
\title{Funktion to compute log(sum(exp(x))}
\usage{
logSumExp(x, mean = F)
}
\arguments{
\item{x}{values at log scale}

\item{mean}{logical, determines whether the mean should be used instead of the sum}
}
\description{
Funktion to compute log(sum(exp(x))
}
\details{
This function computes log(sum(exp(x)), using the offset trick to avoid numeric overflow, see, e.g. http://jblevins.org/notes/log-sum-exp. The mean option allows calculating logMeanExp
}
\author{
Florian Hartig
}
\keyword{internal}
