library(BayesFactor)

## Sleep data from t test example
data(sleep)
plot(extra ~ group, data = sleep)

## Calculating BF for difference values of the prior 

scale = exp(seq(-2,2, len = 50))
bf = rep(NA, 50)
for(i in 1:50) bf[i] = ttestBF(x = sleep$extra[sleep$group==1], y = sleep$extra[sleep$group==2], paired=TRUE, rscale = scale[i])@bayesFactor$bf
plot(scale, bf, xlab= "prior scale", ylab = "BF")


ModelB4 <- generalTestBF(SDI ~ C + A + RS + RT + Gender + Context +
                           
                           C:Gender:Context +
                           
                           A:Gender:Context+
                           
                           RT:Gender:Context+
                           
                           RS:Gender:Context, data= df.exlcude,
                         
                         whichRandom = ID)

