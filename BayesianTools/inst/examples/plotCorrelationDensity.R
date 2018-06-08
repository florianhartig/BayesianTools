dat = generateTestDensityMultiNormal(sigma = "no correlation", sample = TRUE)
correlationPlot(dat(100000))
correlationPlot(dat(100000), scaleCorText = FALSE)
