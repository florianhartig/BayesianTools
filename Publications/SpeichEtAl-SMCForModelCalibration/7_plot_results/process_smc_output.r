process.out <- function(dir, bin.width=20, bin.max=8000){
  
  # This function takes in the output of SMC experiments (including the MCMC reference runs) and returns
  # a list with the following objects:
  #   settings.conv: a data.frame containing the settings (SMC parameters), time and performance of the SMC runs
  #                  (only for those runs which reached convergence)
  #   bins: vector of equal-interval bins used for the display of the MCMC reference curves
  #   min.ref: vector containing the lower limit of the MCMC "polygon" (lowest Ds value at each time step / bin)
  #   max.ref: vector containing the upper limit of the MCMC "polygon" (highest Ds value at each time step / bin)
  
  # Arguments:
  #  dir: path to directory where 
  
  library(BayesianTools)
  library(coda)
  
  # Get a list of all settings combinations
  
  runs <- list.files(paste0(dir,"5/"), pattern = ".RData", full.names = FALSE)
  
  # Read in the outputs (samples) for all settings combinations (there are 5 runs for each combination).
  # Obtain Gelman-Rubin convergence diagnostic between the 5 runs of each combination.
  
  conv <- list()
  
  for(i in 1:length(runs)){
    print(runs[i])
    outs <- list()
    
    for(j in 1:5){
      lpath <- paste0(dir, j, "/")
      loaded <- load(file=paste0(lpath,runs[i]))
      outs[[j]] <- mcmc(get(loaded)$particles)
    }
    
    conv[[i]] <- gelman.diag(outs)  # From coda package
  }
  
  # For each settings combination, get the maximum Gelman-Rubin score for point estimate and upper CI, as well as multivariate Gelman-Rubin.
  
  max.gr.point <- max.gr.uci <- gr.multi <- vector("numeric", length(runs))
  
  for(i in 1:length(runs)){
    max.gr.point[i] <- max(conv[[i]]$psrf[,1])
    max.gr.uci[i] <- max(conv[[i]]$psrf[,2])
    gr.multi[i] <- conv[[i]]$mpsrf
  }
  
  settings <- expand.grid(mutateSteps = c(10, 2, 20, 30, 5), proposalScale = c(0.01, 0.1, 0.333, 0.5), ess.limit = c(0.5, 0.75, 0.9),  particles = c(100000, 20000, 50000, 5000))
  
  settings$max.gr.point <- max.gr.point
  settings$max.gr.uci <- max.gr.uci
  settings$gr.multi <- gr.multi
  settings$ind <- 1:nrow(settings)
  
  
  # Keep only those settings combinations where the runs have converged
  
  settings.conv <- settings[settings$max.gr.point < 1.05 & settings$max.gr.uci < 1.05 & settings$gr.multi < 1.2, ]
  
  # Read in MCMC reference results and SMC results
  mcmc.out.np1 <- read.table(paste0(dir,"1/mcmcRef_np.txt"), header=TRUE)
  mcmc.out.np2 <- read.table(paste0(dir,"2/mcmcRef_np.txt"), header=TRUE)
  mcmc.out.np3 <- read.table(paste0(dir,"3/mcmcRef_np.txt"), header=TRUE)
  mcmc.out.np4 <- read.table(paste0(dir,"4/mcmcRef_np.txt"), header=TRUE)
  mcmc.out.np5 <- read.table(paste0(dir,"5/mcmcRef_np.txt"), header=TRUE)
                            
  smc.out <- list()
  for(j in 1:5){
    smc.out[[j]] <- read.table(paste0(dir, j, "/smc_out.txt"), header=TRUE)
  }
  
  avg.time <- min.time <- max.time <- avg.d <- min.d <- max.d <- vector("numeric", nrow(settings.conv))
  
  
  
  for(i in 1:nrow(settings.conv)){
    time.vec <- d.vec <- vector("numeric", 5)
    ind <- which(smc.out[[1]]$particles == settings.conv[i,]$particles & smc.out[[1]]$ess.limit == settings.conv[i,]$ess.limit & smc.out[[1]]$proposalScale == settings.conv[i,]$proposalScale & smc.out[[1]]$mcmcSteps == settings.conv[i,]$mutateSteps)
    for(j in 1:5){
      print(c(settings.conv[i,]$particles, settings.conv[i,]$ess.limit, settings.conv[i,]$proposalScale, settings.conv[i,]$mutateSteps))
      time.vec[j] <- smc.out[[j]]$time[ind]
      d.vec[j] <- smc.out[[j]]$d[ind]
    }
    avg.time[i] <- mean(time.vec)
    min.time[i] <- min(time.vec)
    max.time[i] <- max(time.vec)
    avg.d[i] <- mean(d.vec)
    min.d[i] <- min(d.vec)
    max.d[i] <- max(d.vec)
  }
  
  settings.conv$avg.time <- avg.time
  settings.conv$min.time <- min.time
  settings.conv$max.time <- max.time
  settings.conv$avg.d <- avg.d
  settings.conv$min.d <- min.d
  settings.conv$max.d <- max.d
  
  ###################
  # Aggregate the 5 MCMC reference runs
  bins <- seq(0,bin.max,by=bin.width)
  min.ref <- max.ref <- vector("numeric", (length(bins)-1))
  
  for(i in 1:length(min.ref)){
    
    all.inbin <- c(mcmc.out.np1$distance[mcmc.out.np1$time > bins[i] & mcmc.out.np1$time <= bins[i+1]],
                   mcmc.out.np2$distance[mcmc.out.np2$time > bins[i] & mcmc.out.np2$time <= bins[i+1]],
                   mcmc.out.np3$distance[mcmc.out.np3$time > bins[i] & mcmc.out.np3$time <= bins[i+1]],
                   mcmc.out.np4$distance[mcmc.out.np4$time > bins[i] & mcmc.out.np4$time <= bins[i+1]],
                   mcmc.out.np5$distance[mcmc.out.np5$time > bins[i] & mcmc.out.np5$time <= bins[i+1]])
    
    min.ref[i] <- min(all.inbin)
    max.ref[i] <- max(all.inbin)
  }
  
  out <- list(settings.conv, bins, min.ref, max.ref)
  return(out)
}

#################################################################
#################################################################

plot.smc.out <- function(out.list, xlim = c(0,7000), ylim = c(0,0.2), main="Model", bars=TRUE, symbol.scale = 1, xlog = "no", polygon.col = adjustcolor("blue", alpha.f = 0.2)){
 
# This function plots one panel of Fig. 3 in the manuscript of Speich et al. 
  
  settings.conv <- out.list[[1]]
  bins <- out.list[[2]]
  min.ref <- out.list[[3]]
  max.ref <- out.list[[4]]
  
  
  settings.conv$col <- NA
  settings.conv$col[settings.conv$proposalScale==0.01] <- "yellow"
  settings.conv$col[settings.conv$proposalScale==0.1] <- "orange"
  settings.conv$col[settings.conv$proposalScale==0.333] <- "red"
  settings.conv$col[settings.conv$proposalScale==0.5] <- "black"
  
  
  settings.conv$cex <- NA
  settings.conv$cex[settings.conv$particles==5000] <- 0.25
  settings.conv$cex[settings.conv$particles==20000] <- 0.5
  settings.conv$cex[settings.conv$particles==50000] <- 1
  settings.conv$cex[settings.conv$particles==100000] <- 2
  settings.conv$cex <- settings.conv$cex * symbol.scale
  
  settings.conv$fill <- NA
  settings.conv$fill[settings.conv$ess.limit==0.5] <- "#edf8fb"
  settings.conv$fill[settings.conv$ess.limit==0.75] <- "#b2e2e2"
  settings.conv$fill[settings.conv$ess.limit==0.9] <- "#66c2a4"
  #settings.conv$fill[settings.conv$ess.factor==0.99] <- "#238b45"
  
  settings.conv$pch <- NA
  settings.conv$pch[settings.conv$mutateSteps == 2] <- 21
  settings.conv$pch[settings.conv$mutateSteps == 5] <- 22
  settings.conv$pch[settings.conv$mutateSteps == 10] <- 23
  settings.conv$pch[settings.conv$mutateSteps == 20] <- 24
  settings.conv$pch[settings.conv$mutateSteps == 30] <- 25
  
  if(xlog=="ln"){
    avg.time <- log(settings.conv$avg.time)
    if(bins[1]==0) bins[1] <- 1
    bins <- log(bins)
    min.time <- log(settings.conv$min.time)
    max.time <- log(settings.conv$max.time)
    if(xlim[2] >= 10000){
      x.labels <- seq(round(xlim[1],-3), xlim[2], by=2000)
      x.ticks <- log(seq(round(xlim[1],-2), xlim[2], by=200))
    } else {
      x.labels <- seq(round(xlim[1],-3), xlim[2], by=1000)
      x.ticks <- log(seq(round(xlim[1],-2), xlim[2], by=100))
    }
    x.pos <- log(x.labels)
    xlim <- log(xlim)
  } else if(xlog=="log10"){
    avg.time <- log10(settings.conv$avg.time)
    if(bins[1]==0) bins[1] <- 1
    bins <- log10(bins)
    min.time <- log10(settings.conv$min.time)
    max.time <- log10(settings.conv$max.time)
    if(xlim[2] >= 10000){
      x.labels <- seq(round(xlim[1],-3), xlim[2], by=2000)
      x.ticks <- log10(seq(round(xlim[1],-2), xlim[2], by=200))
    } else {
      x.labels <- seq(round(xlim[1],-3), xlim[2], by=1000)
      x.ticks <- log10(seq(round(xlim[1],-2), xlim[2], by=100))
    }
    x.pos <- log10(x.labels)
    xlim <- log10(xlim)
  } else if(xlog=="no"){
    avg.time <- settings.conv$avg.time
    min.time <- settings.conv$min.time
    max.time <- settings.conv$max.time
    if(xlim[2] >= 20000){
      x.labels <- seq(round(xlim[1],-3), xlim[2], by=2000)
      x.ticks <- seq(round(xlim[1],-2), xlim[2], by=200)
    } else {
      x.labels <- seq(round(xlim[1],-3), xlim[2], by=1000)
      x.ticks <- seq(round(xlim[1],-2), xlim[2], by=200)
      }
    x.pos <- x.labels
  }
  
  plot(avg.time, settings.conv$avg.d, type="n", ylim=ylim, xlim=xlim, xlab="Runtime [s]", ylab="Distance from reference", main = main, las=1, axes=FALSE, cex.main=0.8)
  #polygon.col <- adjustcolor("blue", alpha.f = 0.2)
  polygon(x = c(bins[1:(length(bins)-1)],rev(bins[1:(length(bins)-1)])), y=c(min.ref, rev(max.ref)), col=polygon.col, border=NA)
  if(bars){
    segments(x0=min.time, y0=settings.conv$avg.d, x1=max.time, lwd=0.5)
    segments(x0 = avg.time, y0 = settings.conv$min.d, y1 = settings.conv$max.d, lwd=0.5)
  }
  points(avg.time, settings.conv$avg.d, pch=settings.conv$pch, col=settings.conv$col, cex=settings.conv$cex, bg=settings.conv$fill)
  
  axis(side = 1, at=x.pos, labels = x.labels, las=3, cex.axis = 0.8)
  axis(side=1, lwd = 0, at=x.ticks, labels = rep("", length(x.ticks)), tck = -0.02, lwd.ticks = 0.5)
  axis(side = 2, las = 1, cex.axis = 0.8)
}


processed.VSEMa <- process.out("/Users/mspeich/Documents/bayes/output/VSEMa")
processed.VSEMb <- process.out("/Users/mspeich/Documents/bayes/output/VSEMb")
processed.threePGN <- process.out("/Users/mspeich/Documents/bayes/output/threePGN")
processed.threePGN_sleep <- process.out("/Users/mspeich/Documents/bayes/output/threePGN_sleep", bin.max=30000, bin.width=300)

proc.out <- list(processed.VSEMa, processed.VSEMb, processed.threePGN, processed.threePGN_sleep)
save(proc.out, file="/Users/mspeich/Documents/bayes/output/processed.RData")

plot.smc.out(processed.VSEMa, main="VSEM - Strong correlation", xlim=c(10,2000), bars=TRUE, symbol.scale = 1.5, xlog="ln")
plot.smc.out(processed.VSEMa, main="VSEM - Strong correlation", xlim=c(10,2000), bars=TRUE, symbol.scale = 1.5, xlog="log10")
plot.smc.out(processed.VSEMb, main="VSEM - No strong correlation", xlim=c(10,5000), bars=FALSE, symbol.scale=1.5, xlog="log10")
plot.smc.out(processed.threePGN, main="3-PGN", xlim=c(200,10000), bars=TRUE, symbol.scale = 1.5, xlog = "log10")

processed.threePGN_sleep <- process.out("/Users/mspeich/Documents/bayes/output/threePGN_sleep")


###################################
### The following code generates the figure shown as Fig. 3 in Speich et al.

pdf(file ="/Users/mspeich/Documents/bayes/figs/smc_compare_results_linTime.pdf", height = 8, width = 9)
#png(file ="/Users/mspeich/Documents/bayes/figs/smc_compare_results_linTime.png", height = 8, width = 9,units="in")

#par(mfrow=c(2,2))
par(mar=c(3.1, 3.1, 1.1, 1.1))
par(fig=c(0.05, 0.525, 0.6, 1))
#plot.smc.out(processed.VSEMa, main="VSEM - Strong correlation", xlim=c(10,5000), bars=FALSE, symbol.scale = 1, xlog="log10", polygon.col = "lightblue")
plot.smc.out(processed.VSEMa, main="VSEM - Strong correlation", xlim=c(10,5000), bars=FALSE, symbol.scale = 1, xlog="no", polygon.col = "lightblue")

par(fig=c(0.525, 1, 0.6, 1), new=TRUE)
#plot.smc.out(processed.VSEMb, main="VSEM - No strong correlation", xlim=c(10,5000), bars=FALSE, symbol.scale=1, xlog="log10", polygon.col = "lightblue")
plot.smc.out(processed.VSEMb, main="VSEM - No strong correlation", xlim=c(10,5000), bars=FALSE, symbol.scale=1, xlog="no", polygon.col = "lightblue")

par(fig=c(0.05, 0.525, 0.2, 0.6), new=TRUE)
#plot.smc.out(processed.threePGN, main="3-PGN", xlim=c(20,10000), bars=FALSE, symbol.scale = 1, xlog = "log10", polygon.col = "lightblue")
plot.smc.out(processed.threePGN, main="3-PGN", xlim=c(20,8000), bars=FALSE, symbol.scale = 1, xlog = "no", polygon.col = "lightblue")


par(fig=c(0.525, 1, 0.2, 0.6), new=TRUE)
#plot.smc.out(processed.sleep, main="3-PGN sleep", xlim=c(10,20000), bars=FALSE, symbol.scale = 1, xlog = "log10", polygon.col = "lightblue") # TEMP
plot.smc.out(processed.threePGN_sleep, main="3-PGN sleep", xlim=c(10,24000), bars=FALSE, symbol.scale = 1, xlog = "no", polygon.col = "lightblue") # TEMP

## LEGEND
par(fig=c(0, 1, 0, 0.15),new=TRUE)
par(mar=c(0,0,0,0))
plot(1,1, type = "n", axes=FALSE, xlab="", ylab="", xlim = c(0,7000), ylim=c(-0.01, 0.05))

points(x=0, y=0.03, pch=21, cex=0.25, col="black", bg="white")
points(x=0, y=0.02, pch=21, cex=0.5, col="black", bg="white")
points(x=0, y=0.01, pch=21, cex=1, col="black", bg="white")
points(x=0, y=0, pch=21, cex=2, col="black", bg="white")
text(x=200, y=0.04, labels="N", font=2)
text(x=400, y=0.03, labels="5000", adj=0, cex=0.8)
text(x=400, y=0.02, labels="20000", adj=0, cex=0.8)
text(x=400, y=0.01, labels="50000", adj=0, cex=0.8)
text(x=400, y=0, labels="100000", adj=0, cex=0.8)

points(x=1500, y=0.03, pch=21, col="black", bg="#edf8fb")
points(x=1500, y=0.015, pch=21, col="black", bg="#b2e2e2")
points(x=1500, y=0, pch=21, col="black", bg="#66c2a4")
text(x=1700, y=0.04, labels="a", font=2)
text(x=1900, y=0.03, labels = "0.5", adj=0, cex = 0.8)
text(x=1900, y=0.015, labels = "0.75", adj=0, cex = 0.8)
text(x=1900, y=0, labels = "0.9", adj=0, cex = 0.8)

points(x=3000, y = 0.03, pch = 21, col = "black", bg = "white")
points(x=3000, y = 0.03 * 0.75, pch = 22, col = "black", bg = "white")
points(x=3000, y = 0.03 * 0.5, pch = 23, col = "black", bg = "white")
points(x=3000, y = 0.03 * 0.25, pch = 24, col = "black", bg = "white")
points(x=3000, y = 0, pch = 25, col = "black", bg = "white")
text(x=3200, y=0.04, labels = "S", font = 2)
text(x=3400, y = 0.03, labels = "2", adj = 0, cex = 0.8)
text(x=3400, y = 0.03 * 0.75, labels = "5", adj = 0, cex = 0.8)
text(x=3400, y = 0.03 * 0.5, labels = "10", adj = 0, cex = 0.8)
text(x=3400, y = 0.03 * 0.25, labels = "20", adj = 0, cex = 0.8)
text(x=3400, y = 0, labels = "30", adj = 0, cex = 0.8)

points(x = 4500, y = 0.03, pch = 21, col = "yellow", bg = "white")
points(x = 4500, y = 0.02, pch = 21, col = "orange", bg = "white")
points(x = 4500, y = 0.01, pch = 21, col = "red", bg = "white")
points(x = 4500, y = 0.00, pch = 21, col = "black", bg = "white")
text(x = 4700, y = 0.04, labels = expression(gamma), font = 2)
text(x = 4900, y = 0.03, labels = "0.01", adj = 0, cex = 0.8)
text(x = 4900, y = 0.02, labels = "0.1", adj = 0, cex = 0.8)
text(x = 4900, y = 0.01, labels = "0.333", adj = 0, cex = 0.8)
text(x = 4900, y = 0, labels = "0.5", adj = 0, cex = 0.8)

rect(xleft=6000, xright = 7000, ybottom = 0.02, ytop=0.03, border = NA, col = "lightblue")
text(x = 6000, y = 0.015, labels = "MCMC reference", cex = 0.8, adj = 0)


## Annotations
par(fig=c(0,1,0,1), new=TRUE)
par(mar=c(0,0,0,0))
plot(1, 1, type="n", axes = FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
text(x = 0.52, y = 0.15, labels = "Runtime [s]")
text(x = 0.01, y = 0.6, labels = "Distance from reference [-]", srt = 90)

dev.off()