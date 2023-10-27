
#' @title BayesianTools
#' @name BayesianTools
#' @docType package
#' @useDynLib BayesianTools, .registration = TRUE
#' @description This package contains general-purpose Markov Chain Monte Carlo 
#' (MCMC) and Sequential Monte Carlo (SMC) samplers, along with diagnostic 
#' functions and plots used for Bayesian statistics.
#' @details This package works particularly for process-based models. 
#' 
#' The package includes two primary functions, \code{\link{createBayesianSetup}}, 
#' which creates a standardized Bayesian setup with likelihood and priors, and 
#' \code{\link{runMCMC}}, which allows various MCMC and SMC samplers to be run.
#' 
#' Additionally, the package can be used for general (non-Bayesian) target 
#' functions.
#' 
#' To use this package, start by creating a BayesianSetup using 
#' \code{\link{createBayesianSetup}}. Generally, a BayesianSetup contains priors
#'  and likelihood densities, or in general, a target function.
#'
#' You can sample with \code{\link{runMCMC}} function. This function can call 
#' several general
#' purpose Metropolis samplers including the \code{\link{Metropolis}}. This 
#' function allows you to specify various popular Metropolis variants such as 
#' adaptive and/or delayed rejection Metropolis; two variants of differential 
#' evolution MCMC - \code{\link{DE}} and \code{\link{DEzs}}, two variants of DREAM - 
#' \code{\link{DREAM}} and \code{\link{DREAMzs}}, the \code{\link{Twalk}} MCMC, 
#' and a Sequential Monte Carlo sampler \code{\link{smcSampler}}. 
#'
#' If a single run is performed, the output of runMCMC is of class mcmcSampler / 
#' smcSampler, otherwise, it would be of class mcmcSamplerList / smcSamplerList.
#' To plot, compare models (DIC, marginal likelihood) or use the output as a 
#' new prior, several functions are available.  
#'
#' To learn how to use the package, run the following command: 
#' vignette("BayesianTools", package="BayesianTools").
#' 
#' Run citation("BayesianTools") to cite the package.
#'
#' To report bugs or ask for help, post at
#' \href{https://stackoverflow.com/questions/5963269/how-to-make-a-great-r-reproducible-example}{reproducible example} 
#' via the BayesianTools \href{https://github.com/florianhartig/BayesianTools/issues}{issue tracker} 
#' on GitHub. 
#'
#'Acknowledgements: The creation and maintenance of this package profited from 
#'funding and collaboration through Cost Action FP 1304 PROFOUND, DFG	DO 786/12-1 CONECT, EU FP7 ERA-NET Sumforest REFORCE and Bayklif Project BLIZ. 
NULL