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
#' @param singlePanel logical, determining whether the parameter should be plotted in a single panel or each in its own panel
#' @param settings optional list of additional settings for density plot, histogram, and violin plot, respectively
#' @param ... additional arguments passed to getSample
marginalPlot <- function(x, prior = TRUE, xrange = NULL, type = 'd', singlePanel = TRUE, settings = NULL, ...) {
  nPriorDraws <- 10000
  
  posteriorMat <- getSample(x, parametersOnly = TRUE, ...)
  
  priorMat <- if (!is.null(prior) & prior != FALSE) {
    if ('bayesianOutput' %in% class(x)) bs <- getSetup(x)$prior$sampler(nPriorDraws) # draw prior from bayesianSetup
    else if (any(c('data.frame', 'matrix') %in% class(prior))) {
      if (ncol(posteriorMat) == ncol(prior)) prior
      else stop('posterior and prior must have the same number of parameters/columns')
    } else stop('prior must be matrix or data.frame, or NULL, if x is matrix/data.frame')
  } else NULL
  
  if (!is.null(priorMat)) colnames(priorMat) <- colnames(posteriorMat)
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
  
  
  .args <- c(list(posteriorMat = posteriorMat, priorMat = priorMat, xrange = xrange, singlePanel = singlePanel),
             settings)
  if (type == 'd') do.call(marginalPlotDensity, .args)
  else if (type == 'h') do.call(marginalPlotHistogram, .args)
  else if (type == 'v') do.call(marginalPlotViolin, .args)
}


#' Plot marginals as densities
#' @param posteriorMat matrix with samples as rows and parameters as columns
#' @param priorMat matrix (optional) with samples as rows and parameters as columns
#' @param xrange vector or matrix (optional), determining the plotting range, with parameters as columns and min, max as rows
#' @param col vector of colors for posterior and
#' @param singlePanel logical, determining whether the parameter should be plotted in a single panel or each in its own panel
#' @param ... further options
#' @author Tankred Ott
#' @keywords internal
marginalPlotDensity <- function(posteriorMat, priorMat = NULL, xrange = NULL, col=c('#FF5000A0','#4682B4A0'),
                                singlePanel = TRUE, ...) {
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
  
  postXY <- lapply(posteriorDensities, function(d) {
    xy <- cbind(c(d$x[1], d$x, d$x[length(d$x)]),
                c(0, d$y, 0))
    colnames(xy) <- c('x', 'y')
    xy
  })
  
  priorXY <- if (!is.null(priorDensities)) lapply(priorDensities, function(d) {
    xy <- cbind(c(d$x[1], d$x, d$x[length(d$x)]),
                c(0, d$y, 0))
    colnames(xy) <- c('x', 'y')
    xy
  }) else NULL
  
  
  # TODO: put singlePanel FALSE and TRUE in same loop, only change plot() parameters (do.call() ?) 
  
  if (singlePanel) {
    op <- par(mfrow = c(nPar,1), mar = c(2, 5, 2, 2), oma = c(4, 4, 4, 0))
    on.exit(par(op))
    
    for (i in 1:length(posteriorDensities)) {
      postX <- postXY[[i]][,1]
      postY <- postXY[[i]][,2]
      
      priorX <- if (!is.null(priorXY[[i]])) priorXY[[i]][,1] else NULL
      priorY <- if (!is.null(priorXY[[i]])) priorXY[[i]][,2] else NULL
      
      yrange <- if (is.null(priorX)) range(postY) else range(c(postY, priorY))
      
      plot(NULL, NULL, xlim = xrange[,i], ylim = yrange, main = NA,
           xlab = NA, ylab = NA, bty = 'n', yaxt = 'n')
      mtext(sprintf('%20s', parNames[i]), 2, las = 1, adj = 1.25)
      
      
      polygon(postX, postY, col = col[1], border = NA)
      if (!is.null(priorX)) polygon(priorX, priorY, col = col[2], border = NA)
      
    }
    
    mtext('Marginal parameter uncertainity', outer = TRUE, cex = 1.5)
    
    # overlay plot with empty plot to be able to place the legends freely
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
    
    legend('bottom', if (!is.null(priorX)) c('posterior', 'prior') else 'posterior', xpd = TRUE, horiz = TRUE,
           inset = c(0, 0), bty = 'n', pch = 15, col = col, cex = 1.5)
    
    
  } else {
    mfrow <- if (nPar <= 16) BayesianTools:::getPanels(nPar) else c(4,4)
    print(mfrow)
    
    op <- par(mfrow = mfrow, mar = c(4.5, 4.5, 5, 3), oma=c(3, 0, 2, 0), xpd=TRUE)
    
    on.exit(par(op))
    for (i in 1:length(posteriorDensities)) {
      postX <- postXY[[i]][,1]
      postY <- postXY[[i]][,2]
      
      priorX <- if (!is.null(priorXY[[i]])) priorXY[[i]][,1] else NULL
      priorY <- if (!is.null(priorXY[[i]])) priorXY[[i]][,2] else NULL
      
      yrange <- if (is.null(priorX)) range(postY) else range(c(postY, priorY))
      
      plot(NULL, NULL, xlim = xrange[,i], ylim = yrange, main = parNames[i],
           xlab = parNames[i], ylab = 'density')
      
      polygon(postX, postY, col = col[1], border = NA)
      if (!is.null(priorX)) polygon(priorX, priorY, col = col[2], border = NA)
      
    }
    mtext('Marginal parameter uncertainity', outer = TRUE, cex = 1.5)
    
    # overlay plot with empty plot to be able to place the legends freely
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
    
    legend('bottom', if (!is.null(priorX)) c('posterior', 'prior') else 'posterior', xpd = TRUE, horiz = TRUE, inset = c(0, 0),
           bty = 'n', pch = 15, col = col, cex = 1.5)
  }
  
}

#' @keywords internal
marginalPlotHistogram <- function(posteriorMat, priorMat, xrange, ...) {
  stop('Histogram not implemented')
}

#' @keywords internal
marginalPlotViolin <- function(posteriorMat, priorMat, xrange, ...) {
  stop('Violin plot not implemented')
}



singlePanel <- T
type <- 'd'
marginalPlot(out, xrange = matrix(c(-1, 10, -8, 5), nrow = 2), type=type, settings = list(n = 512, cut = T), singlePanel = singlePanel)
marginalPlot(out, prior = F, type=type, settings = list(n = 256, cut = F), singlePanel = singlePanel)
marginalPlot(matrix(c(1,2,3,4,5, 2,3,4,5,6), ncol = 2),
             prior = NULL,
             type=type, singlePanel = singlePanel)
marginalPlot(matrix(c(1,2,3,4,5, 2,3,4,5,6), ncol = 2),
             prior = matrix(c(-5,2,3,4,15,6,7,8, -3,3,4,5,9,0,1,7), ncol = 2),
             type=type, singlePanel = singlePanel)
# marginalPlot(out, )

plot(ds[[1]])
