#' Determine the groups of correlated parameters
#' @author Stefan Paul
#' @param chain MCMC chain including only the parameters (not logP,ll, logP)
#' @param blockSettings a list with settings
#' @return groups
#' @keywords internal
updateGroups <- function(chain,blockSettings){
  
  settings <- getBlockSettings(blockSettings)
  blockUpdateType  <- settings$blockUpdateType
  
  switch(blockUpdateType,
         "correlation" = {
         ## (Pair wise) Correlation in the parameters
         cormat <- abs(cor(chain[,1:(ncol(chain)-3),sample(1:dim(chain)[3],1)]))
         diag(cormat) <- 0
         # Correct for NA and Inf values as this could cause error in as.dist()
         cormat[c(which(is.na(cormat)),which(cormat == Inf),which(cormat == -Inf)) ] <- 0
         tree <- hclust(as.dist(1-cormat))  # get tree based on distance(dissimilarity = 1-cor).
         cT <- cutree(tree, k = settings$k, h = settings$h) # get groups. With h we can manipulate the strength of the interaction.
         },
         "user" = { 
          cT <-  settings$groups
           },
         "random" = {
           pool <- c(1:settings$k, sample(1:settings$k, (ncol(chain)-3-settings$k)))
           cT <- sample(pool)
         }
  )
  
  pSel <- settings$pSel
  if(is.null(pSel) && is.null(settings$pGroup)) pSel = rep(1,ncol(chain)-3)
  return(list(cT = cT, pGroup = settings$pGroup, pSel = pSel))
}


#' Determine the parameters in the block update
#' @param blockSettings settings for block update
#' @return vector containing the parameter to be updated
#' @keywords internal
getBlock <- function(blockSettings){
  groups <- blockSettings$cT
  pGroup <- blockSettings$pGroup
  pSel <- blockSettings$pSel
  
  
  nGroups = max(groups)
  if(nGroups == 1) return(1:length(groups))
  if (is.null(pGroup)) pGroup = rep(1,nGroups)
  if(length(pSel) > nGroups) pSel <- pSel[1:nGroups]
  pSel = c(pSel, rep(0,nGroups - length(pSel)))
  groupsToSample = sample.int(nGroups, 1, prob = pSel)
  
  selectedGroups = sample.int(nGroups,groupsToSample, prob = pGroup[1:nGroups])
  GroupMember  <- which(is.element(groups,selectedGroups))
  return(GroupMember)
  
}


#' getblockSettings
#' @description Transforms the original settings to settings used in the model runs
#' @param blockUpdate input settings
#' @return list with block settings
#' @keywords internal
getBlockSettings <- function(blockUpdate){
  
    h <- k <- pSel <- pGroup <- groups <- NULL  
    blockUpdateType <- blockUpdate[[1]]
    
    switch(blockUpdateType,
           "correlation" = {
             h <- blockUpdate$h
             k <- blockUpdate$k
             pSel <- blockUpdate$pSel
             pGroup <- blockUpdate$pGroup 
           },
           "random"={
             k <- blockUpdate$k 
           },
           "user"= {
             groups <- blockUpdate$groups
             pSel <- blockUpdate$pSel
             pGroup <- blockUpdate$pGroup 
           })
    
  return(list(blockUpdateType = blockUpdateType, h = h, k = k, pSel = pSel,
              pGroup = pGroup, groups = groups))
  }


