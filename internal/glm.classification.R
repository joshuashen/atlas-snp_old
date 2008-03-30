

y <- read.table("train.list.2", header=TRUE)


z <- glm(y$class ~ y$swap + y$homo + y$qual + y$dist3 + y$dist5, family=binomial )

summary(z)

beta = z$coefficients

y.n <- read.table("test.list.2", header=TRUE)
dd = length(y.n$class)
ymatrix = cbind(c(rep(1,dd)), y.n$swap, y.n$homo, y.n$qual, y.n$dist3, y.n$dist5)
ypredict = ymatrix %*% beta

pp = exp(ypredict) / (1 + exp(ypredict))

pp.t <- pp[y.n$class > 0 ]
pp.f <- pp[y.n$class < 1 ]


