
x = runif(500,-1,1)
y = 0.2 + 0.9  *x + rnorm(500, sd = 0.5)

summary(lm(y ~ x))

GOF(x,y)

GOF(x,y, plot = TRUE)
