
tms <- function(x, trim=0.25) {
  xs <- sort(x^2)
  n <- length(x)
  return( mean( xs[1:floor(n*(1-trim))] ) )
}

source('fasts.R')
data(Boston, package='MASS')
x <- model.matrix(medv ~ ., data=Boston)
y <- Boston$medv
# sest <- fast.s(x=x, y=y, int=0, N=500, k=2, best.r=5, seed=123)
n <- dim(x)[1]

# 10-fold CV
k <- 10
set.seed(123)
ii <- sample( (1:n) %% k )
pr.s <- pr.ls <- rep(NA, n)
for(j in 1:k) {
  trs <- (ii != j)
  tr.x <- x[ trs, ]
  tr.y <- y[ trs ]
  ses <- fast.s(x=tr.x, y=tr.y, int=0, N=500, k=2, best.r=5, seed=123)
  lses <- lm(medv ~ ., data=Boston, subset = trs)
  pr.ls[ !trs ] <- predict(lses, newdata = Boston[ !trs, ])
  pr.s[ !trs ] <- as.vector( x[ !trs, ] %*% ses[[1]] )
}
tms( (y - pr.ls) )
tms( (y - pr.s) )


