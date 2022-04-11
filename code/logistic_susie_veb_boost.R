#TODO: add maxit argument

logistic.susie.veb.boost <- function(X, y, L=10){
  require(VEB.Boost)
  veb.fit <- VEB.Boost::veb_boost_stumps(
    X, y, k=L,
    family='binomial',
    include_stumps = F,
    growTree=F,
    changeToConstant=F,
    max_log_prior_var = 35,
    scale_X = 'NA'
  )
  alpha <- t(do.call(cbind, lapply(veb.fit$leaves, function(x) x$currentFit$alpha)))
  mu <- t(do.call(cbind, lapply(veb.fit$leaves, function(x) x$currentFit$mu)))
  mu2 <- t(do.call(cbind, lapply(veb.fit$leaves, function(x) x$currentFit$mu2)))
  elbo <- veb.fit$ELBO_progress[[2]]
  res <- list(alpha=alpha, mu=mu, mu2=mu2, elbo=elbo, veb.fit=veb.fit)
  
  class(res) <- 'susie'
  colnames(res$alpha) <- colnames(X)
  colnames(res$mu) <- colnames(X)
  res$pip <- susieR::susie_get_pip(res)
  names(res$pip) <- colnames(X)
  res$sets <- susieR::susie_get_cs(res, X=X)
  return(res)
}
