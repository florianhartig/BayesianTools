# 
# ############### MCMC RESULT #############################
# 
# 
# #' runs an automatic MCMC
# #' 
# #' @export
# automaticMCMC <- function(likelihood, prior = NULL, startvalues, maxiter=20000, steplength = 1000, optimize = T){
#   
#   mcmc1 <- mcmcSampler(likelihood = likelihood, prior = prior, startvalue = startvalues[[1]], optimize = optimize)
#   mcmc2 <- mcmcSampler(likelihood = likelihood, prior = prior, startvalue = startvalues[[2]], optimize = optimize)
#   mcmc3 <- mcmcSampler(likelihood = likelihood, prior = prior, startvalue = startvalues[[3]], optimize = optimize)
#   
#   steps = 0
#   conv = F
#   convTemp = F
#   while(nrow(mcmc1$chain) < maxiter){
#     mcmc1<- getSamples(mcmc1, steplength)
#     mcmc2<- getSamples(mcmc2, steplength)    
#     mcmc3<- getSamples(mcmc3, steplength)
#     
#     sel <- round(steps/3):steps
#     parcol <- 1:mcmc1$numPars
#     parLL <- mcmc1$numPars + 1
#     parLP <- mcmc1$numPars + 1
#     
#     res <- mcmc.list(mcmc(mcmc1$chain[sel, parcol]), mcmc(mcmc2$chain[sel, parcol]), mcmc(mcmc3$chain[sel, parcol]))
#     
#     likelihoods <- c(mcmc1$chain[sel, parLL], mcmc2$chain[sel, parLL], mcmc3$chain[sel, parLL])
#     posteriors <- c(mcmc1$chain[sel, parLP], mcmc2$chain[sel, parLP], mcmc3$chain[sel, parLP])    
#     
#     currentConv <- tryCatch(
# {
#   x <- all(gelman.diag(res)$psrf[1:mcmc1$numPars,] < 1.05)
#   ifelse(is.na(x), F, x)
# },
# error=function(cond) {
#   F
# })
# 
# print(currentConv)
# 
# if (convTemp & currentConv ) {
#   conv = T
#   break
# }  
# convTemp = currentConv
# currentConv = F
# steps = steps + steplength
#   }
# if (! conv) print("Algorithm not converged")
# 
# 
# calc.dic <- function(x,lik,lik.fun,...)
#   dic = calc.dic(x = combineChains(res),lik = likelihoods,lik.fun = mcmc1$catchingLikelihood)
# 
# # Not working yet
# #marginalLik = marginal.likelihood(x = res, lik = likelihoods,lik.fun = mcmc1$catchingLikelihood, prior.fun = mcmc1$catchingPrior)
# 
# 
# classFields = list(
#   mcmclist = res,
#   dic = dic
# )
# 
# class(classFields) <- append(class(classFields),"mcmcResult")
# return(classFields)
# }
