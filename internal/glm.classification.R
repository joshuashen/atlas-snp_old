

y <- read.table("lists.10000", header=TRUE)


z <- glm(y$class ~ y$swap + y$homo  + y$qual + y$dist3 , family=binomial )

summary(z)

beta = z$coefficients

y.n <- read.table("lists.10000.2", header=TRUE)
ymatrix = cbind(c(rep(1,10000)), y.n$swap, y.n$homo, y.n$qual, y.n$dist3)
ypredict = ymatrix %*% beta

pp = exp(ypredict) / (1 + exp(ypredict))

pp.t <- pp[y.n$class > 0 ]
pp.f <- pp[y.n$class < 1 ]


