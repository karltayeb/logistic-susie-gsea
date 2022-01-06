library(mr.ash.alpha)

fit.mr.ash.lasso <- function(X, y, ...){
  lasso <- cv.glmnet(X, y, family = 'gaussian', alpha=1.0, lower.limits = 0)
  c <- coef(lasso, s="lambda.min")
  X <- as.matrix(X)
  y <- as.matrix(y)

  mr.ash.fit <- mr.ash(
    X, y,
    beta.init = tail(as.vector(c), -1), ...)

  mr.ash.post <- mr.ash.alpha::get.full.posterior(mr.ash.fit)
  pip <- 1 - mr.ash.post$phi[, 1]
  beta <- mr.ash.alpha::coef.mr.ash(mr.ash.fit)
  names(beta) <- c('Intercept', names(pip))
  return(tibble(model=list(mr.ash.fit), pip=list(pip), beta=list(beta)))
}

fit.mr.ash <- function(X, y, ...){
  X <- as.matrix(X)
  y <- as.matrix(y)
  mr.ash.fit <- mr.ash(
    X, y, ...)

  mr.ash.post <- mr.ash.alpha::get.full.posterior(mr.ash.fit)
  pip <- 1 - mr.ash.post$phi[, 1]
  beta <- mr.ash.alpha::coef.mr.ash(mr.ash.fit)
  names(beta) <- c('Intercept', names(pip))
  return(tibble(model=list(mr.ash.fit), pip=list(pip), beta=list(beta)))
}
