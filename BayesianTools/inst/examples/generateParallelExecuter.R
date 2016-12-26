
testDensityMultiNormal <- generateTestDensityMultiNormal()


parDen <- generateParallelExecuter(testDensityMultiNormal)$parallelFun
x = matrix(runif(9,0,1), nrow = 3)
parDen(x)

