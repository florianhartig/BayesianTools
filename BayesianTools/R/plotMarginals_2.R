library(BayesianTools)

ll <- function(x) sum(mvtnorm::dmvnorm(x, mean = c(5, -2), sigma = diag(2), log = T))

bs <- createBayesianSetup(ll, lower = c(-10, -10), upper = c(20, 20))

out <- runMCMC(bs, sampler = 'Metropolis', settings = list(nrChains = 1))

plot(out)

s <- getSample(out, numSamples = 100)


#' Plot MCMC marginals
#' @param x bayesianOutput, or matrix or data.frame containing with samples as rows and parameters as columns
#' @param prior logical determining whether the prior should be plotted, or if x is matrix oder data.frame, a matrix of prior draws with draws as rows and parameters as columns
#' @param xrange vector or matrix of plotting ranges for the x axis. If matrix, the rows must be parameters and the columns min and max values.
#' @param type character determining the plot type. Either 'd' for density plot, 'h' for histogram, or 'v' for violin plot
#' @param settings optional list of additional settings for density plot, histogram, and violin plot, respectively
#' @param ... additional arguments passed to getSample
marginalPlot <- function(x, prior = TRUE, xrange = NULL, type = 'd', settings = NULL, ...) {
  nPriorDraws <- 10000
  
  posteriorMat <- getSample(x, parametersOnly = TRUE, ...)
  
  priorMat <- if (!is.null(prior)) {
    if ('bayesianOutput' %in% class(x)) bs <- getSetup(x)$prior$sampler(nPriorDraws) # draw prior from bayesianSetup
    else if (any(c('data.frame', 'matrix') %in% class(prior))) {
      if (ncol(posteriorMat) == ncol(prior)) prior
      else stop('posterior and prior must have the same number of parameters/columns')
    } else stop('prior must be matrix or data.frame, or NULL, if x is matrix/data.frame')
  } else NULL
  
  colnames(priorMat) <- colnames(posteriorMat)
  nPar <- ncol(posteriorMat)
  
  
  if (!is.null(xrange)) {
    if (!any(c('numeric', 'matrix') %in% class(xrange))) stop('xrange must be numeric or matrix, or NULL')
    if ('numeric' %in% class(xrange)) xrange <- matrix(rep(xrange), nPar, nrow = 2)
    else if ('matrix' %in% class(xrange)) {
      if (ncol(xrange) != ncol(posteriorMat)) stop('xrange must have as many colums as there are parameterss')
      else if (nrow(xrange) != 2) stop('xrange must have two rows (min, max)')
    }
  } else {
    posteriorRanges <- apply(posteriorMat, 2, range)
    priorRanges <- if(!is.null(priorMat)) apply(priorMat, 2, range) else NULL
    
    xrange <- if (is.null(priorRanges)) posteriorRanges else apply(rbind(priorRanges, posteriorRanges), 2, range)
  }
  
  
  if (any(c('d', 'dens', 'density') == type)) type <- 'd'
  else if (any(c('h', 'hist', 'histogram') == type)) type <- 'h'
  else if (any(c('v', 'violin') == type)) type <- 'v'
  else stop('type must be one of "d", "h", "v"')
  
  
  .args <- c(list(posteriorMat, priorMat, xrange), settings)
  if (type == 'd') do.call(marginalPlotDensity, .args)
  else if (type == 'h') do.call(marginalPlotHistogram, .args)
  else if (type == 'v') do.call(marginalPlotViolin, .args)
}


#' Plot marginals as densities
#' @param posteriorMat matrix with samples as rows and parameters as columns
#' @param priorMat matrix (optional) with samples as rows and parameters as columns
#' @param xrange vector or matrix (optional), determining the plotting range, with parameters as columns and min, max as rows
#' @author Tankred Ott
#' @keywords internal
marginalPlotDensity <- function(posteriorMat, priorMat = NULL, xrange = NULL, col=c("#FF5000A0","#4682B4A0"), ...) {
  print('d')
  
  nPar <- ncol(posteriorMat)
  parNames <- colnames(posteriorMat)
  
  if (is.null(xrange)) {
    posteriorRanges <- apply(posteriorMat, 2, range)
    priorRanges <- if(!is.null(priorMat)) apply(priorMat, 2, range) else NULL

    xrange <- if (is.null(priorRanges)) posteriorRanges else apply(rbind(priorRanges, posteriorRanges), 2, range)
  }

  # print(xrange)
  
  posteriorDensities <- lapply(1:ncol(posteriorMat),
                               function(i) density(posteriorMat[,i], from = xrange[1,i], to = xrange[2,i], ...))
  priorDensities <- if (!is.null(priorMat)) lapply(1:ncol(priorMat),
                                                   function(i) density(priorMat[,i], from = xrange[1,i], to = xrange[2,i], ...))
                    else NULL
  
  mfrow <- BayesianTools:::getPanels(nPar)
  op <- par(mfrow = mfrow, mar=c(4.5, 4.5, 5, 3))
  on.exit(par(op))
  for (i in 1:length(posteriorDensities)) {
    postX <- posteriorDensities[[i]]$x
    postX <- c(postX[1], postX, postX[length(postX)])
    postY <- posteriorDensities[[i]]$y
    postY <- c(0, postY, 0)
    
    priorX <- NULL
    priorY <- NULL
    if (!is.null(priorDensities[[i]])) {
      priorX <- priorDensities[[i]]$x
      priorX <- c(priorX[1], priorX, priorX[length(priorX)])
      priorY <- priorDensities[[i]]$y
      priorY <- c(0, priorY, 0)
    }
    
    yrange <- if (is.null(priorX)) range(postY) else range(c(postY, priorY))
    
    plot(NULL, NULL, xlim = xrange[,i], ylim = yrange, main = parNames[i],
         xlab = parNames[i], ylab = 'density')
  
    polygon(postX, postY, col = col[1], border = NA)
    if (!is.null(priorX)) polygon(priorX, priorY, col = col[2], border = NA)
    
    legend('topright',
           legend = if(is.null(priorX)) 'posterior' else c('posterior', 'prior'),
           fill = col, bty = 'n', border = NA, inset = c(0,0))
  }
  
  # TODO: LABELS, AXIS, MAIN, LEGEND
}

#' @keywords internal
marginalPlotHistogram <- function(posteriorMat, priorMat, xrange, ...) {
  stop('Histogram not implemented')
}

#' @keywords internal
marginalPlotViolin <- function(posteriorMat, priorMat, xrange, ...) {
  stop('Violin plot not implemented')
}




marginalPlot(out, xrange = matrix(c(-1, 10, -8, 5), nrow = 2), type='d', settings = list(n = 512, cut = T))
marginalPlot(out, type='d', settings = list(n = 256, cut = F))
marginalPlot(matrix(c(1,2,3,4,5, 2,3,4,5,6), ncol = 2),
             prior = NULL,
             type='d')
marginalPlot(matrix(c(1,2,3,4,5, 2,3,4,5,6), ncol = 2),
             prior = matrix(c(-5,2,3,4,15,6,7,8, -3,3,4,5,9,0,1,7), ncol = 2),
             type='d')
# marginalPlot(out, )

plot(ds[[1]])
