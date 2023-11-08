## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(fig.width=5, fig.height=5, warning=FALSE, cache = F)

## ----echo = F, message = F----------------------------------------------------
set.seed(123)

## ----eval = F-----------------------------------------------------------------
#  runMyModel(par)

## ----eval = F-----------------------------------------------------------------
#  dyn.load(model)
#  
#  runMyModel(par){
#    out = # model call here
#    # process out
#    return(out)
#  }

## ----eval = F-----------------------------------------------------------------
#  runMyModel(par){
#  
#    # Create here a string with what you would write to call the model from the command line
#    systemCall <- paste("model.exe", par[1], par[2])
#  
#    out = system(systemCall,  intern = TRUE) # intern indicates whether to capture the output of the command as an R character vector
#  
#    # write here to convert out in the apprpriate R classes
#  
#  }

## ----eval = F-----------------------------------------------------------------
#  runMyModel(par, returnData = NULL){
#  
#    writeParameters(par)
#  
#    system("Model.exe")
#  
#    if(! is.null(returnData)) return(readData(returnData)) # The readData function will be defined later
#  
#  }
#  
#  writeParameters(par){
#  
#    # e.g.
#    # read template parameter fil
#    # replace strings in template file
#    # write parameter file
#  }

## ----eval = F-----------------------------------------------------------------
#  setUpModel <- function(parameterTemplate, site, localConditions){
#  
#    # create the runModel, readData functions (see later) here
#  
#    return(list(runModel, readData))
#  
#  }

## ----eval = F-----------------------------------------------------------------
#  getData(type = X){
#  
#    read.csv(xxx)
#  
#    # do some transformation
#  
#    # return data in desidered format
#  }

## ----eval = F-----------------------------------------------------------------
#  par = c(1,2,3,4 ..)
#  
#  runMyModel(par)
#  
#  output <- getData(type = DesiredType)
#  
#  plot(output)

## -----------------------------------------------------------------------------
mymodel<-function(x){
  output<-0.2*x+0.1^x
  return(output)
}

## ----eval = F-----------------------------------------------------------------
#  
#  library(parallel)
#  cl <- makeCluster(2)
#  
#  runParallel<- function(parList){
#    parSapply(cl, parList, mymodel)
#  }
#  
#  runParallel(c(1,2))

## ----eval = F-----------------------------------------------------------------
#  library(BayesianTools)
#  parModel <- generateParallelExecuter(mymodel)

## -----------------------------------------------------------------------------
library(BayesianTools)

ll <- generateTestDensityMultiNormal(sigma = "no correlation")
bayesianSetup <- createBayesianSetup(likelihood = ll, lower = rep(-10, 3), upper = rep(10, 3))

settings = list(iterations = 200)

# run the several MCMCs chains either in seperate R sessions, or via R parallel packages
out1 <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
out2 <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)

res <- createMcmcSamplerList(list(out1, out2))
plot(res)

