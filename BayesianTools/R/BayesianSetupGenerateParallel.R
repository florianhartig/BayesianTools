#' Factory to generate a parallel executer of an existing function
#' 
#' @author Florian Hartig
#' @param fun function to be changed to parallel exectution
#' @param parallel should a parallel R cluster be used or not. If set to T, cores will be detected automatically and n-1 of the available n cores of the machine will be used. Alternatively, you can set the number of cores used by hand
#' @param parallelOptions list containing three lists. First "packages" determines the R packages necessary to run the likelihood function. Second "variables" the objects in the global envirnment needed to run the likelihood function and third "dlls" the DLLs needed to run the likelihood function (see Details).
#' @note Can also be used to make functions compatible with library sensitivity
#' @details For parallelization, option T means that an automatic parallelization via R is attempted, or "external", in which case it is assumed that the likelihood is already parallelized. In this case it needs to accept a matrix with parameters as columns.
#' Further you can specify the packages, objects and DLLs that are exported to the cluster. 
#' By default a copy of your workspace is exported. However, depending on your workspace this can be very inefficient. 
#'
#' Alternatively you can specify the environments and packages in the likelihood function (e.g. BayesianTools::VSEM() instead of VSEM()).
#' @export
#' @example /inst/examples/generateParallelExecuter.R

generateParallelExecuter <- function(fun, parallel = F, parallelOptions = list(variables = "all", packages = "all", dlls = NULL)){
  
  if (parallel == F){
    parallelFun <- function(parMat, ...){
      res <- apply(parMat, 1, fun, ...)
      if(! is.null(dim(res))) res = t(res) # to have results row-wise if multiple results are returned
      return(res)
    } 
    cl <- "Cluster not defined for bayesianSetup if parallel = FALSE"
  }else{
    
    #library(foreach)
    #library(iterators)
   # library(parallel)
    
    if (parallel == T | parallel == "auto"){
      cores <- parallel::detectCores() - 1
    } else if (is.numeric(parallel)){
      cores <- parallel
      if (cores > parallel::detectCores()) stop("BayesianTools: error - more cores specified than available on this machine")
    } else stop("BayesianTools: error wrong argument to parallel") 
    
    # get variables, packages, dlls in current workspace here if defaults are set in parameters
    
    cl <- parallel::makeCluster(cores)
    
    
    # update the parallelOptions based on user settings.
    defaultParallelOptions <- list(variables = "all", packages = "all", dlls = NULL)
    parallelOptions <- modifyList(defaultParallelOptions, parallelOptions)
    
    
    # get loaded packages
    if(is.null(parallelOptions$packages[1])) packages <- parallelOptions$packages
    else if(parallelOptions$packages[1] == "all") packages <- (.packages())
    else packages <- parallelOptions$packages
    
    # get loaded DLLs
    if(is.null(parallelOptions$dlls[1])) dlls <- parallelOptions$dlls
    else if(parallelOptions$dlls[1] == "all"){
      tmpdlls <- getLoadedDLLs()
      dlls <- vector(mode = "character", length = length(tmpdlls))
      counter <- 0
      for(i in tmpdlls){
        counter <- counter+1
        dlls[counter] <- i[[2]]
      }
    }else dlls <- unlist(parallelOptions$dlls)
 
    
    # get objects in global environment
    if(is.null(parallelOptions$variables[1])) objects = NULL
    else if(parallelOptions$variables[1] == "all") objects <- ls(envir = .GlobalEnv)
    else objects <- unlist(parallelOptions$variables)

    # function to export packages and dlls
    packageFun <- function(packages = NULL, dlls = NULL) {
      if(!is.null(packages)){
      for(i in packages) library(i, character.only = TRUE)
      }
      if(!is.null(dlls)){
       for(i in dlls) try(dyn.load(i), silent = T)
      }
    }
  
    # export packages, dlls and objects to cluster
    parallel::clusterCall(cl, packageFun, packages, dlls)
    parallel::clusterExport(cl, varlist = objects)
    
    #doParallel::registerDoParallel(cl)
    

    
    
    parallelFun <- function(parMat, ...){
      res = parallel::parApply(cl = cl, parMat, 1, fun, ...)
      if(! is.null(dim(res))) res = t(res) # to have results row-wise if multiple results are returned
      return(res)
    }
#     parallelFun <- function(parMat){
#       res <- foreach::foreach(parMat=iter(parMat, by='row'), .combine = "rbind", .packages = "BayesianTools")%dopar%{
#         fun
#       }
#       if(! is.null(dim(res))) res = t(res) # to have results row-wise if multiple results are returned
#       return(res)
#    }
    message("parallel function execution created with", cores, "cores.")
    
  } 
  
  return(list(parallelFun = parallelFun, cl = cl))
}















