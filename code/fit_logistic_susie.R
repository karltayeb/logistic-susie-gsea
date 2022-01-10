library(VEB.Boost)
#source("~/Research/susieR_logistic_wflow/code/logistic_susie_VB_functions.R")

# Logistic SuSiE
get.logistic.susie.alpha <- function(veb.fit) {
  return(t(do.call(cbind, lapply(veb.fit$leaves, function(x) x$currentFit$alpha))))
}

# Logistic SuSiE
get.logistic.susie.expected.effect <- function(veb.fit) {
  return(rowSums(do.call(cbind, lapply(veb.fit$leaves, function(x) x$currentFit$alpha * x$currentFit$mu))))
}

get.logistic.susie.pip <- function(veb.fit){
  alpha <- get.logistic.susie.alpha(veb.fit)
  return(apply(alpha, 2, function (x) 1 - exp(sum(log(1 - x + 1e-10)))))
}

get.logistic.coef <- function(veb.fit){
  ALPHAS <- do.call(
    cbind, lapply(veb.fit$leaves, function(x) x$currentFit$alpha))
  MUS <- do.call(
    cbind, lapply(veb.fit$leaves, function(x) x$currentFit$mu))
  coef <- rowSums(ALPHAS * MUS)
  return(coef)
}

get.pip <- function(alpha){
  return(apply(alpha, 2, function (x) 1 - exp(sum(log(1 - x + 1e-10)))))
}

fit.logistic.susie.veb.boost <- function(X, y, ...) {
  veb.fit = veb_boost_stumps(
    X, y, family = 'binomial',
    include_stumps=FALSE,
    growTree=FALSE, changeToConstant=F, ...)
  alpha <- get.logistic.susie.alpha(veb.fit)
  colnames(alpha) <- colnames(X)
  pip <- get.pip(alpha)
  cs <- alpha2cs(alpha)
  coef <- get.logistic.coef(veb.fit)
  return(tibble(
    model=list(veb.fit),
    alpha=list(alpha),
    cs=list(cs),
    pip=list(pip),
    coef=list(coef)))
}
