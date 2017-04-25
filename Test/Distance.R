#' calculate distance between samples 
#' 
#' 
#' @param sample1 a matrix. If the distance function is not symmetric, this is assumed to be a sample from the target (reference) distribution 
#' @param sample2 a matrix. If the distance function is not symmetric, this is assumed to be a sample from the "other" distribution
#' @param type the distance function
#' 
#' @details Currently, the following distance functions are implemented
#' 
#' KL = Kullback-Leibler Divergence
#' BH = Bhattacharyya distance
#' D =  normalized Euclidean distance between mean and standard deviation of sample and target. Not symmetric. Target is sample 1. This was described in eq. 10 in Laloy, E., and J. A. Vrugt. 2012. High-dimensional posterior exploration of hydrologic models using multiple-try DREAM(ZS) and high-performance computing. Water Resour. Res. 48(1)
#' 
getSampleDistance <- function(sample1, sample2, type = "KL"){
  
  if(type == "KL"){
    library(FNN) # I'm keeping this at the moment to throw an error if FNN is not installed
    x = FNN::KL.dist(sample1, sample2, k=10)
    out = mean(x) # FH: no idea if the mean is a good idea. KL.dist returns one value per cluster size, I don't know which value is best chosen
    
  } else if(type == "BH" ){
    Sigma1 = cov(sample1)
    Sigma2 = cov(sample2)
    mu1 = colMeans(sample1)
    mu2 = colMeans(sample2)
    
    # The following code is copied from package fpc
    
    aggregatesigma <- (Sigma1+Sigma2)/2
    d1 <- mahalanobis(mu1,mu2,aggregatesigma)/8
    d2 <- log(det(as.matrix(aggregatesigma))/sqrt(det(as.matrix(Sigma1))*
                                                    det(as.matrix(Sigma2))))/2
    out <- d1+d2
    
    # end fpc
    
    out
    
  } else if(type == "D" ){
    
    sd1 = apply(sample1, 2, sd)
    sd2 = apply(sample2, 2, sd)
    mu1 = colMeans(sample1)
    mu2 = colMeans(sample2)
    dev = sum( ((mu1 - mu2)/sd1)^2 + ((mu1 - mu2)/sd1)^2  )
    
    out = sqrt( 1/(2*length(mu1)) * dev )
    
    
  } else if(type == "xxxx" ){
    
  } else if(type == "xxxx" ){
    
  } else if(type == "xxxx" ){
    
  } else stop("unrecognized argument for type in BayesianTools::getDistanceDistributions")
  
  return(out)
  
}

#' 
#' KL Distance 
#' There are at least three options in R
#' 1.https://artax.karlin.mff.cuni.cz/r-help/library/LaplacesDemon/html/KLD.html
#' 2. http://svitsrv25.epfl.ch/R-doc/library/flexmix/html/KLdiv.html
#' 3. FNN::KL.dist, which is used here
#' 
#' Mahalanobis distance
#' 
#' is implemented in stats. 
#' 


library(mvtnorm)

sigma <- matrix(c(4,2,2,3), ncol=2)
X <- rmvnorm(n=5000, mean=c(1,2), sigma=sigma)
Y <- rmvnorm(n=5000, mean=c(5,2), sigma=sigma)


getSampleDistance(X,Y, type = "KL")
getSampleDistance(X,Y, type = "BH")
getSampleDistance(X,Y, type = "D")



