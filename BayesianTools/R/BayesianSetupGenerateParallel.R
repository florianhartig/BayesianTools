#' Factory to generate a parallel executor of an existing function
#' 
#' @author Florian Hartig
#' @param fun function to be changed to parallel execution
#' @param parallel should a parallel R cluster be used? If set to T, the operating system will automatically detect the available cores and n-1 of the available n cores will be used. Alternatively, you can manually set the number of cores to be used
#' @param parallelOptions a list containing three lists. \n First, "packages": determines the R packages required to run the likelihood function. \n Second, "variables": the objects in the global environment needed to run the likelihood function \n Third, "dlls": the DLLs needed to run the likelihood function (see Details).
#' @note can be used to make functions compatible with library sensitivity
#' @details For parallelization, if option T is selected, an automatic parallelization is tried via R. Alternatively, "external" can be selected on the assumption that the likelihood has already been parallelized. In the latter case, a matrix with parameters as columns must be accepted. You can also specify which packages, objects and DLLs are exported to the cluster. By default, a copy of your workspace is exported, but depending on your workspace, this can be inefficient. As an alternative, you can specify the environments and packages in the likelihood function (e.g. BayesianTools::VSEM() instead of VSEM()).
#' @export
#' @example /inst/examples/generateParallelExecuter.R
#' 
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















