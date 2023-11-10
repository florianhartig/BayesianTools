######
# Twalk helper functions
######


#' Wrapper for step function
#' @param Npar number of parameters
#' @param FUN log posterior density
#' @param x parameter vector of chain 1
#' @param Eval last evaluation of x
#' @param x2 parameter vector of chain 2
#' @param Eval2 last evaluation of x
#' @param at "traverse" move proposal parameter.
#' @param aw "walk" move proposal parameter.
#' @param pn1 Probability determining the number of parameters that are changed.
#' @param Ptrav Move probability of "traverse" moves, default to 0.4918
#' @param Pwalk Move probability of "walk" moves, default to 0.4918
#' @param Pblow Move probability of "blow" moves, default to 0.0082
#' @param Phop Move probability of "hop" moves
#' @references Christen, J. Andres, and Colin Fox. "A general purpose sampling algorithm for continuous distributions (the t-walk)." Bayesian Analysis 5.2 (2010): 263-281.
#' @keywords internal
TwalkMove <- function (Npar, FUN, x, Eval, x2, Eval2, at = 6, aw = 1.5, pn1 = min(Npar, 4)/Npar,
                  Ptrav = 0.4918, Pwalk = 0.4918, Pblow = 0.0082, Phop = 0.0082) 
{

  p <- sample(4,1, prob = c(Ptrav,Pwalk,Pblow,Phop))
  
  if(p == 1)case <- "traverse"
  else if(p ==2) case <- "walk"
  else if(p ==3) case <- "blow"
  else case <- "hop"
  
  
  out <- Twalksteps(case = case, Npar = Npar, FUN = FUN, x = x,
               Eval = Eval, x2 = x2, Eval2 = Eval2, at = at, aw = aw, pn1 = pn1)
  
  
  return(list(y = out$y, val = out$val, y2 = out$y2, val2 = out$val2, alpha = out$alpha))
}


 
#' Main function that is executing and evaluating the moves
#' @param case Type of Twalk move. Either "walk", "traverse", "hop" or "blow"
#' @param Npar number of parameters
#' @param FUN Log posterior density
#' @param x parameter vector of chain 1
#' @param Eval last evaluation of x
#' @param x2 parameter vector of chain 2
#' @param Eval2 last evaluation of x
#' @param at "traverse" move proposal parameter.
#' @param aw "walk" move proposal parameter.
#' @param pn1 Probability determining the number of parameters that are changed.
#' @references Christen, J. Andres, and Colin Fox. "A general purpose sampling algorithm for continuous distributions (the t-walk)." Bayesian Analysis 5.2 (2010): 263-281.
#' @keywords internal
Twalksteps <- function(case, Npar, FUN, x,
                  Eval, x2, Eval2, at, aw, pn1){
  
  val <- NULL
  val2 <- NULL
  p <- runif(1)
  
  switch(case,
         "traverse" = { #Traverse
           if (p < 0.5) {
             beta <- betaFun(at)
             tmp <- propFun(case, Npar = Npar, pn1 = pn1, x2 = x, x = x2, beta = beta)
             y2 <- tmp$prop
             npSel <- tmp$npSel
             y <- x
             val <- Eval
             val2 <- FUN(y2, returnAll = T)
             
             if (npSel == 0) alpha <- 1
             else alpha <- exp((- Eval2[1] + val2[1]) + (npSel - 2) * log(beta))
             
           }else{
             beta <- betaFun(at)
             tmp <- propFun(case, Npar = Npar, pn1 = pn1, x = x, x2 = x2, beta = beta)
             y <- tmp$prop
             npSel <- tmp$npSel
             y2 <- x2
             val2 <- Eval2
             
             val <- FUN(y, returnAll = T)
    
             if (npSel == 0) alpha <- 1
             else alpha <- exp((-Eval[1] + val[1]) + (npSel - 2) * log(beta))
             
           }}, # End traverse
         "walk" = { # walk
           if (p < 0.5) {
             tmp <- propFun(case, Npar = Npar, pn1 = pn1, aw = aw, x2 = x, x = x2)
             y2 <- tmp$prop
             npSel <- tmp$npSel
             y <- x
             val <- Eval
             if ( (all(abs(y2 - y) > 0))) {
               val2 <- FUN(y2, returnAll = T)

               alpha <- exp(-Eval2[1] + val2[1])
             }
             else {
               alpha <- 0
             }
           }else{
             tmp <- propFun(case, Npar = Npar, pn1 = pn1, aw = aw, x = x, x2 = x2)
             y <- tmp$prop
             npSel <- tmp$npSel
             y2 <- x2
             val2 <- Eval2
             if ( (all(abs(y2 - y) > 0))) {
               val <- FUN(y, returnAll = T)

               alpha <- exp(-Eval[1] + val[1])
             }
             else {
               alpha <- 0
             }
         }}, # End walk
         "blow" = { #blow
           if (p < 0.5) {
             tmp <- propFun(case, Npar = Npar, pn1 = pn1, x = x2, x2 = x)
             y2 <- tmp$prop
             npSel <- tmp$npSel
             pSel <- tmp$pSel
             y <- x
             val <- Eval
             if ( all(y2 != x)) {
               val2 <- FUN(y2, returnAll = T)

               G1 <- Gfun(case, npSel, pSel, y2, x2, x)
               G2 <- Gfun(case, npSel, pSel, x2, y2, x)
               alpha <- exp((-Eval2[1] + val2[1]) + (G1 - G2))
             }
             else {
               alpha <- 0
             }
           }else{
             tmp <- propFun(case, Npar = Npar, pn1 = pn1, x = x, x2 = x2)
             y <- tmp$prop
             npSel <- tmp$npSel
             pSel <- tmp$pSel
             y2 <- x2
             val2 <- Eval2
             if (all(y != x2)) {
               val <- FUN(y, returnAll = T)

               G1 <- Gfun(case, npSel, pSel, y, x, x2)
               G2 <- Gfun(case, npSel, pSel, x, y, x2)
               alpha <- exp((-Eval[1] + val[1]) + (G1 - G2))
             }
             else {
               alpha <- 0
             }
           }
         }, # End blow
         "hop" = { #hop
           if (p < 0.5) {
             tmp <- propFun(case, Npar = Npar, pn1 = pn1, x2 = x, x = x2)
             y2 <- tmp$prop
             npSel <- tmp$npSel
             pSel <- tmp$pSel
             y <- x
             val <- Eval
             if ( all(y2 != x)) {
               val2 <- FUN(y2, returnAll = T)

               G1 <- Gfun(case, npSel, pSel, y2, x2, x)
               G2 <- Gfun(case, npSel, pSel, x2, y2, x)
               alpha <- exp((-Eval2[1] + val2[1]) + (G1 - G2))
             }
             else {
               alpha <- 0
             }
           }else{
             tmp <- propFun(case, Npar = Npar, pn1 = pn1, x = x, x2 = x2)
             y <- tmp$prop
             npSel <- tmp$npSel
             pSel <- tmp$pSel
             y2 <- x2
             val2 <- Eval2
             if ( all(y != x2)) {
               val <- FUN(y, returnAll = T)

               G1 <- Gfun(case, npSel, pSel, y, x, x2)
               G2 <- Gfun(case, npSel, pSel, x, y, x2)
               alpha <- exp((-Eval[1] + val[1]) + (G1 - G2))
             }
             else {
               alpha <- 0
             }
           
         }}) # End hop and end switch
  return(list(y = y, val = val, y2 = y2, val2 = val2, alpha = alpha, 
              npSel = npSel))
}







################## Helper functions 
###############################################################

#' Helper function for sum of x*x
#' @param x vector of values
#' @keywords internal
sumSquare <- function(x){return(sum(x*x))}

 
#' Helper function to create proposal
#' @param case  Type of Twalk move. Either "walk", "traverse", "hop" or "blow"
#' @param Npar number of parameters
#' @param pn1 Probability determining the number of parameters that are changed.
#' @param aw "walk" move proposal parameter.
#' @param beta parameter for "traverse" move proposals.
#' @param x parameter vector of chain 1
#' @param x2 parameter vector of chain 2
#' @references Christen, J. Andres, and Colin Fox. "A general purpose sampling algorithm for continuous distributions (the t-walk)." Bayesian Analysis 5.2 (2010): 263-281.
#' @keywords internal
propFun <- function(case, Npar, pn1, x, x2, beta = NULL, aw = NULL){
  
  switch(case,
         "traverse"={
           pSel <- (runif(Npar) < pn1) 
           prop <- NULL 
           for (i in 1:Npar){
             if (pSel[i]) prop <- c( prop, x2[i] + beta*(x2[i] - x[i]))
             else prop <- c( prop, x[i]) 
           }
           return(list(prop=prop, npSel=sum(pSel))) 
         },
         "walk"={
           u <- runif(Npar) 
           pSel <- (runif(Npar) < pn1) 
           z <- (aw/(1+aw))*(aw*u^2 + 2*u -1) 
           z <- z*pSel 
           return(list( prop=x + (x - x2)*z, npSel=sum(pSel))) 
         },
         "blow"={
           pSel <- (runif(Npar) < pn1) 
           sigma <- max(pSel*abs(x2 - x)) 
           return(list( prop=x2*pSel + sigma*rnorm(Npar)*pSel + x*(1-pSel), npSel=sum(pSel), pSel=pSel)) 
           
         },
         "hop"={
           pSel <- (runif(Npar) < pn1) 
           sigma <- max(pSel*abs(x2 - x))/3 
           prop <- NULL 
           for (i in 1:Npar){
             if (pSel[i]) prop <- c( prop, x[i] + sigma*rnorm(1))
             else prop <- c( prop, x[i]) 
           }
           return(list( prop=prop, npSel=sum(pSel), pSel=pSel)) 
           
         }
         
  )
}

 
#' Helper function for calculating beta
#' @param at "traverse" move proposal parameter.
#' @keywords internal
betaFun <- function(at)
{
  if (runif(1) < (at-1)/(2*at)) return(exp(1/(at + 1)*log(runif(1))))
  else return(exp(1/(1 - at)*log(runif(1)))) 
}

 
#' Helper function for blow and hop moves
#' @param case  Type of Twalk move. Either "hop" or "blow"
#' @param npSel number of parameters that are changed.
#' @param pSel vector containing information about which parameters are changed.
#' @param h Parameter for "blow" and hop moves
#' @param x parameter vector of chain 1
#' @param x2 parameter vector of chain 2
#' @references Christen, J. Andres, and Colin Fox. "A general purpose sampling algorithm for continuous distributions (the t-walk)." Bayesian Analysis 5.2 (2010): 263-281.
#' @keywords internal
Gfun <- function(case, npSel, pSel, h, x, x2){
  switch(case, 
         "blow"= {
           sigma <- max(pSel*abs(x2 - x))
           if(npSel > 0) return((npSel/2)*log(2*pi) + npSel*log(sigma) + 0.5*sumSquare(h - x2)/(sigma^2))
           else return(0)
         },
         "hop" = {
           sigma <- max(pSel*abs(x2 - x))/3 
           if (npSel > 0) return((npSel/2)*log(2*pi) - npSel*log(3) + npSel*log(sigma) + 0.5*9*sumSquare((h - x))/(sigma^2))
           else return(0)
         })
  
}



