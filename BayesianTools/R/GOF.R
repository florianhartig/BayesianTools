#' Standard GOF metrics
#' Startvalues for sampling with nrChains > 1 : if you want to provide different start values for the different chains, provide a list
#' @author Florian Hartig
#' @param predicted predicted values
#' @param observed observed values
#' @param plot should a plot be created
#' @param centered if T, variables are centered to the mean of the observations, i.e. the intercept is for the mean value of the observation
#' 
#' @details The function considers observed ~ predicted and calculates
#' 
#' 1) rmse = root mean squared error
#' 2) mae = mean absolute errorr
#' 3) a linear regression with slope, intercept and coefficient of determination R2
#' 
#' For the linear regression, centered = T means that variables will be centered around the mean value of the observation. This setting avoids a correlation between slope and intercept (that the intercept is != 0 as soon as the slope is !=0)   
#' 
#' @note In principle, it is possible to plot observed ~ predicted and predicted ~ observed. However, if we assume that the error is mainly on the y axis (observations), i.e. that observations scatter around the true (ideal) value, we should plot observed ~ predicted. See Pineiro et al. (2008). How to evaluate models: observed vs. predicted or predicted vs. observed?. Ecological Modelling, 216(3-4), 316-322.
#' @return A list with the following entries: rmse = root mean squared error, mae = mean absolute error, slope = slope of regression, offset = intercept of regression, R2 = R2 of regression  
#' @example /inst/examples/GOF.R
#' @export
GOF<- function(observed, predicted, plot = F, centered = T){
  
  # root mean squared error
  rmse <- sqrt( mean( (predicted - observed)^2, na.rm = T) )
  
  # mean absolute error
  mae <- mean( abs(predicted - observed), na.rm = TRUE) 
  
  #ssq <- sum( (predicted - observed)^2, na.rm= T) 

  # linear regression 
  
  if(centered == T){
    meanObs = mean(observed, na.rm = T)
    observed = observed - meanObs
    predicted = predicted - meanObs
  }
  
  linReg = lm( observed ~ predicted) 
  
  if(plot == T){
    plot(observed ~ predicted)
    abline(linReg, col = "red")
    abline(v = 0, lty = 2)
    abline(h = 0 , lty = 2)
    #loessMod <- mgcv::gam(predicted ~ s(observed))
    #ord = order(predicted)
    #lines(predicted[ord], predict(loessMod)[ord], col = "green")
  }
  
  out = list(rmse = rmse, mae = mae, slope = as.numeric(coefficients(linReg)[2]), offset = as.numeric(coefficients(linReg)[1]), R2 = summary(linReg)$r.squared)
  return(out)
}
