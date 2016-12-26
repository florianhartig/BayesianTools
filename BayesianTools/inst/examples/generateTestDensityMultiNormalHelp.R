# sampling from the test function
x = generateTestDensityMultiNormal(sample  = TRUE, n = 1000)(1000)
correlationPlot(x)
marginalPlot(x)

# generating the the density
density = generateTestDensityMultiNormal(sample  = FALSE)
density(x[1,])
