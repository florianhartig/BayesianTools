
#' Very simple ecosystem model
#' @description A very simple ecosystem model, based on three carbon pools and a basic LUE model 
#' @param pars a parameter vector with parameters and initial states
#' @param PAR Forcing, photosynthetically active radiation (PAR) MJ /m2 /day
#' @param C switch to choose whether to use the C or R version of the model. C is much faster. 
#' @return a matrix with columns NEE, CV, CR and CS units and explanations see details
#' @import Rcpp
#' @useDynLib BayesianTools, .registration = TRUE
#' @details This Very Simple Ecosystem Model (VSEM) is a 'toy' model designed to be very simple but yet bear some resemblance to deterministic processed based ecosystem models (PBMs) that are commonly used in forest modelling.
#' 
#' The model determines the accumulation of carbon in the plant and soil from the growth of the plant via photosynthesis and senescence to the soil which respires carbon back to the atmosphere.
#' 
#' The model calculates Gross Primary Productivity (GPP) using a very simple light-use efficiency (LUE) formulation multiplied by light interception. Light interception is calculated via Beer's law with a constant light extinction coefficient operating on Leaf Area Index (LAI).
#' 
#' A parameter (GAMMA) determines the fraction of GPP that is autotrophic respiration. The Net Primary Productivity (NPP) is then allocated to above and below-ground vegetation via a fixed allocation fraction. Carbon is lost from the plant pools to a single soil pool via fixed turnover rates. Heterotrophic respiration in the soil is determined via a soil turnover rate.
#' 
#' The model equations are 
#' 
#' -- Photosynthesis
#' 
#' \deqn{LAI = LAR*Cv}
#' \deqn{GPP = PAR * LUE * (1 - \exp^{(-KEXT * LAI)})}
#' \deqn{NPP = (1-GAMMA) * GPP}
#' 
#' -- State equations
#' \deqn{dCv/dt  = Av * NPP - Cv/tauV}
#' \deqn{dCr/dt  = (1.0-Av) * NPP - Cr/tauR}
#' \deqn{dCs/dt  = Cr/tauR + Cv/tauV - Cs/tauS}
#' 
#' The model time-step is daily.
#' 
#' -- VSEM inputs:
#' 
#' PAR    Photosynthetically active radiation (PAR) MJ /m2 /day
#' 
#' -- VSEM parameters:
#' 
#' KEXT   Light extinction coefficient m2 ground area / m2 leaf area
#' 
#' LAR    Leaf area ratio m2 leaf area / kg aboveground vegetation
#' 
#' LUE    Light-Use Efficiency (kg C MJ-1 PAR)
#' 
#' GAMMA  Autotrophic respiration as a fraction of GPP
#' 
#' tauV   Longevity of aboveground vegetation days
#' 
#' tauR   Longevity of belowground vegetation days
#' 
#' tauS   Residence time of soil organic matter d
#' 
#' -- VSEM states:
#' 
#' Cv     Above-ground vegetation pool kg C / m2
#' 
#' Cr     Below-ground vegetation pool kg C / m2
#' 
#' Cs     Carbon in organic matter kg C / m2
#' 
#' -- VSEM fluxes:
#' 
#' G     Gross Primary Productivity kg C /m2 /day
#' 
#' NPP   Net Primary Productivity kg C /m2 /day
#' 
#' NEE   Net Ecosystem Exchange kg C /m2 /day
#' @seealso \code{\link{VSEMgetDefaults}}, \code{\link{VSEMcreatePAR}}, , \code{\link{VSEMcreateLikelihood}} 
#' @example /inst/examples/VSEMHelp.R
#' @export
#' @author David Cameron, R and C implementation by Florian Hartig
VSEM <- function(pars =  c(KEXT = 0.5,
                          LAR  = 1.5,
                          LUE = 0.002,
                          GAMMA = 0.4,
                          tauV = 1440,
                          tauS = 27370,
                          tauR = 1440,
                          Av = 0.5,
                          Cv = 3,
                          Cs = 15,
                          Cr = 3), 
                 PAR, C = TRUE){ 
  
  if (C == T){
    out <- vsemC(pars, PAR)
    colnames(out) = c("NEE", "Cv", "Cs", "CR")
    
    return(out)    

  } else {
    
    numObs = length(PAR)
    KEXT  = pars[1]
    LAR   = pars[2]
    LUE   = pars[3]
    GAMMA = pars[4]
    tauV  = pars[5]
    tauS  = pars[6]
    tauR  = pars[7] 
    Av    = pars[8] 
    Cv    = pars[9] 
    Cs    = pars[10] 
    Cr    = pars[11] 
    out   = matrix(nrow = numObs, ncol = 4 )
    colnames(out) = c("NEE", "Cv", "Cs", "CR")
    for (i in 1:numObs){
      G   = PAR[i] * LUE * (1 - exp(-KEXT*LAR*Cv)) 
      NPP = (1-GAMMA)*G
      Cv  = Cv + Av*NPP - Cv/tauV
      Cr  = Cr + (1.0-Av)*NPP - Cr/tauR
      Cs  = Cs + Cr/tauR + Cv/tauV - Cs/tauS
      NEE = (Cs/tauS + GAMMA*G) - G    
      out[i, ] = c(NEE, Cv, Cs, Cr)
    }
    return(out)
  }
}

#' returns the default values for the VSEM
#' @export
#' @return a data.frame
VSEMgetDefaults <- function(){
  
  best = list(    KEXT  = 0.5,
                          LAR   = 1.5,
                          LUE   = 0.002,
                          GAMMA = 0.4,
                          tauV  = 1440,
                          tauS  = 27370,
                          tauR  = 1440, 
                          Av    = 0.5,
                          Cv    = 3.0,
                          Cs    = 15,
                          Cr    = 3.0
  )
  def = data.frame(best = unlist(best))
  def$lower = c(0.2,0.2,0.0005, 0.2, 500,4000,500, 0.2, 0,0,0)
  def$upper= c(1,3,0.004, 0.6, 3000,50000,3000, 1, 400,1000,200)  
  return(def)
}

#' Allows to mix a given parameter vector with a default parameter vector
#' @param pars vector with new parameter values
#' @param defaults vector with defaukt parameter values
#' @param locations indices of the new parameter values
#' @rdname package-deprecated
#' @description This function is deprecated and will be removed by v0.2. 
#' @export
createMixWithDefaults <- function(pars, defaults, locations){
  .Deprecated(package = "BayesianTools")
  out = defaults
  out[locations] = pars
  return(out)
}


#' Create a random radiation (PAR) time series
#' @author David Cameron, R implementation by Florian Hartig
#' @param days days to calculate the PAR for 
#' @export
VSEMcreatePAR <- function(days = 1:(3*365)){
  PAR = (abs (sin(days/365 * pi)+ rnorm(length(days)) *0.25)) *10
  return(PAR)
}


#' Create an example dataset, and from that a likelihood or posterior for the VSEM model
#' @author Florian Hartig
#' @param likelihoodOnly logical, decides whether to create only a likelihood, or a full bayesianSetup with uniform priors.
#' @param plot logical, decides whether data should be plotted
#' @param selection vector containing the indices of the selected parameters
#' @details The purpose of this function is to be able to conveniently create a likelihood for the VSEM model for demonstration purposes. The function creates example data --> likelihood --> BayesianSetup, where the latter is the 
#' @export
VSEMcreateLikelihood <- function(likelihoodOnly = F, plot = F, selection =  c(1:6, 12)){
  
  # create radiation input 
  PAR <- VSEMcreatePAR(1:1000)
  plotTimeSeries(observed = PAR)
  
  # create reference parameters and add one row for the SD of the observed data
  
  refPars <- VSEMgetDefaults()
  refPars[12,] <- c(0.1, 0.001, 0.5)
  rownames(refPars)[12] <- "error-sd"
  
  # create reference data 
  
  referenceData <- VSEM(refPars$best[1:11], PAR) 
  obs <- referenceData + rnorm(length(referenceData), sd = (abs(referenceData) + 1E-7) * refPars$best[12])
  
  # plot if that switch is on
  
  if(plot == T){
    oldpar <- par(mfrow = c(2,2))
    for (i in 1:4) plotTimeSeries(observed = obs[,i], predicted = referenceData[,i], main = colnames(referenceData)[i])
    par(oldpar)
  }
  
  # Create likelihood for reference data
  
  likelihood <- function(x){
    mix = refPars$best
    mix[selection] = x
    
    predicted <- VSEM(mix[1:11], PAR)
    predicted[,1] <- predicted[,1]
    diff <- c(predicted[,1:3] - obs[,1:3])
    llValues <- dnorm(diff, sd = (abs(c(predicted[,1:3])) + 0.0000001) * mix[12], log = T) 
    if (sum == FALSE) return(llValues)
    else return(sum(llValues))
  }

  if(likelihoodOnly == T) return(likelihood)
  else{

  bayesianSetup <- createBayesianSetup(likelihood, lower = refPars$lower[selection], upper = refPars$upper[selection] , best = refPars$best[selection], names = rownames(refPars)[selection])

  return(bayesianSetup)
  
  }
}


