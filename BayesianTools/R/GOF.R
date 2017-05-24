#' Standard GOF metrics
#' Startvalues for sampling with nrChains > 1 : if you want to provide different start values for the different chains, provide a list
#' @author Florian Hartig
#' @param predicted predicted values
#' @param observed observed values
#' @param plot should a plot be created
#' 
#' @details The function considers predicted ~ observed and calculates
#' 
#' 1) rmse = root mean squared error
#' 2) mae = mean absolute errorr
#' 3) a linear regression with slope, intercept and coefficient of determination R2
#' 
#' For the linear regression, the x axis is centered, meaning that the intercept is the difference between observed / predicted for the MEAN predicted value. This setting avoids a correlation between slope and intercept (that the intercept is != 0 as soon as the slope is !=0)   
#' 
#' @return A list with the following entries: rmse = root mean squared error, mae = mean absolute error, slope = slope of regression, offset = intercept of regression, R2 = R2 of regression  
#' @example /inst/examples/GOF.R
#' @export
GOF<- function(observed, predicted, plot = F){
  
  # root mean squared error
  rmse <- sqrt( mean( (predicted - observed)^2, na.rm = T) )
  
  # mean absolute error
  mae <- mean( abs(predicted - observed), na.rm = TRUE) 
  
  # linear regression 
  linReg = lm( predicted ~ scale(observed, scale = F) ) 
  

  
  #   scale(predicted) ~ scale(observed)     scale(predicted) ~ scale(observed)
  
  
  #ssq <- sum( (predicted - observed)^2, na.rm= T) 
  
  if(plot == T){
    plot(predicted ~ scale(observed, scale = F))
    abline(linReg, col = "red")
    abline(v = 0, lty = 2)
    abline(h = 0 , lty = 2)
  }
  
  out = list(rmse = rmse, mae = mae, slope = as.numeric(coefficients(linReg)[2]), offset = as.numeric(coefficients(linReg)[1]), R2 = summary(linReg)$r.squared)
  return(out)
}
