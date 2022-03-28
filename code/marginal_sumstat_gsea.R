sigmoid = function(x){
  1 / (1 + exp(-x))
}
#' make all the functions for performing summary stat gsea
#' beta_i | se_i, x, theta, alpha  ~ pi_i N_1 + (\pi_i - 1) N_0
summary.stat.gsea = function(x, beta, se){
  # add column of ones for intercept
  X <- cbind(rep(1, length(x)), x)

  #' likelihood function (component assignment marginalized out)
  #' params a list with alpha, sigma0 and theta
  #' theta, alpha, sigma0 are optional arguments to mask value in params
  likelihood = function(params=NULL, theta=NULL, alpha=NULL, sigma0=NULL){
    if(is.null(theta)){
      theta = params$theta
    }
    if(is.null(alpha)){
      alpha = params$alpha
    }
    if(is.null(sigma0)){
      sigma0 = params$sigma0
    }
    pi = sigmoid(theta[1] + theta[2] * x)
    sigma = se^alpha * sigma0
    null = dnorm(beta, mean=0, sd=se, log=TRUE) + log(1-pi)
    non_null = dnorm(beta, mean=0, sd=sqrt(se^2 + sigma^2), log=TRUE) + log(pi)
    ll = sum(purrr::map2_dbl(null, non_null, ~ matrixStats::logSumExp(c(.x, .y))))
    return(ll)
  }

  #' Q function which is maximized in M step (expectation over assignments)
  #' params a list with alpha, sigma0 and theta
  #' theta, alpha, sigma0 are optional arguments to mask value in params
  Q = function(params=NULL, theta=NULL, alpha=NULL, sigma0=NULL, responsibilities=NULL){
    if(is.null(theta)){
      theta = params$theta
    }
    if(is.null(alpha)){
      alpha = params$alpha
    }
    if(is.null(sigma0)){
      sigma0 = params$sigma0
    }
    if(is.null(responsibilities)){
      if(is.null(params$responsibilities)){
        params <- update.responsibilities(params)
      }
      responsiblities <- params$responsibilities
    }
    logit.pi = theta[1] + theta[2] * x
    pi = sigmoid(logit.pi)
    sigma = se^alpha * sigma0
    null = dnorm(beta, mean=0, sd=se, log=TRUE) + log(1-pi)
    non_null = dnorm(beta, mean=0, sd=sqrt(se^2 + sigma^2), log=TRUE) + log(pi)
    return(sum(responsiblities * non_null + (1-responsiblities) * null))
  }
  
  #' compute posterior probability of coming from the non-null component
  #' conditioned on current parameter settings
  update.responsibilities = function(params, log.odds=FALSE){
    logit.pi = params$theta[1] + params$theta[2] * x
    sigma = se^params$alpha * params$sigma0
    null = dnorm(beta, mean=0, sd=se, log=TRUE)
    non_null = dnorm(beta, mean=0, sd=sqrt(se^2 + sigma^2), log=TRUE)
    responsibilities = (non_null - null) + logit.pi
    if(!log.odds){
      responsibilities = sigmoid(responsibilities)
    }
    # update parameters
    new.params = params
    new.params$responsibilities = responsibilities
    return(new.params)
  }
  
  #' factory returns function
  #' optimization objective for theta with fixed responsibilities
  #' and other parameters fixex
  theta.obj.factory = function(params){
    y <- params$responsibilities
    theta.obj = function(theta) {
      logit.pi = params$theta[1] + params$theta[2] * x
      pi = sigmoid(logit.pi)
      #null = dnorm(beta, mean=0, sd=se, log=TRUE) + log(1 - pi)
      #non_null = dnorm(beta, mean=0, sd=sqrt(se^2 + sigma^2), log=TRUE) + log(pi)
      #ll = sum(purrr::map2_dbl(null, non_null, ~ matrixStats::logSumExp(c(.x, .y))))
      #return(ll)
      return(sum(y * log(pi) + (1-y) * log(1-pi)))
    }
    return(theta.obj)
  }

  #' factory returns function
  #' gradient of function returned by `theta.lik.factory`
  theta.grad.factory = function(params){
    y <- params$responsibilities
    theta.grad = function(theta){
      # compute grad
      logit.pi = params$theta[1] + params$theta[2] * x
      pi = sigmoid(logit.pi)
      grad <- c(sum(y - pi), sum(x * (y - pi)))
      return(grad)
    }
  }

  #' factory returns function
  #' hessian of function returned by `theta.lik.factory`
  theta.hess.factory = function(params){
    y <- params$responsibilities
    theta.hess = function(theta){
      # compute grad
      logit.pi = params$theta[1] + params$theta[2] * x
      pi = sigmoid(logit.pi)
      W = pi * (1-pi)
      hess = - t(X) %*% apply(X, 2, function(x) x*pi)
      return(hess)
    }
  }

  #' routine for updating theta
  #' implemented with optim and above objective/grad function
  update.theta.grad = function(params){
    theta.opt = optim(
      params$theta,
      theta.obj.factory(params),
      theta.grad.factory(params),
      control = list(fnscale=-1)
    )
    # update parameters
    new.params = params
    new.params$theta = theta.opt$par
    return(new.params)
  }
  
  #' routine for updating theta w.r.t Q
  update.theta.Q = function(params){
    # function in terms of log alpha, log sigma0
    lik.theta = function(par){
      theta = par
      return(Q(params, theta=theta))
    }
    # optimize
    res = optim(
      params$theta,
      lik.theta,
      control = list(fnscale=-1)
    )
    # update parameters
    new.params = params
    new.params$theta <- res$par
    browser()
    return(new.params)
  }
  update.theta = update.theta.Q


  #' routine for updating alpha and sigma0 w.r.t marginal likelihood
  update.alpha.sigma0 = function(params){
    # function in terms of log alpha, log sigma0
    lik.lnalpha.lntheta = function(par){
      lnalpha = par[1]
      lnsigma0 = par[2]
      return(likelihood(params, alpha=exp(lnalpha), sigma0=exp(lnsigma0)))
    }
    # optimize
    res = optim(
      c(-1, -1),
      lik.lnalpha.lntheta,
      control = list(fnscale=-1)
    )
    # update parameters
    new.params = params
    new.params$alpha <- exp(res$par[1])
    new.params$sigma0 <- exp(res$par[2])
    return(new.params)
  }

  #' routine for updating alpha w.r.t marginal likelihood
  update.alpha = function(params){
    # function in terms of log alpha, log sigma0
    lik.lnalpha= function(par){
      lnalpha = par[1]
      return(likelihood(params, alpha=exp(lnalpha)))
    }
    # optimize
    res = optim(
      c(-1),
      lik.lnalpha,
      control = list(fnscale=-1)
    )
    # update parameters
    new.params = params
    new.params$alpha <- exp(res$par[1])
    return(new.params)
  }

  #' routine for updating sigma0 w.r.t marginal likelihood
  update.sigma0 = function(params){
    # function in terms of log sigma0, log sigma0
    lik.lnsigma0= function(par){
      lnsigma0 = par[1]
      return(likelihood(params, sigma0=exp(lnsigma0)))
    }
    # optimize
    res = optim(
      c(-1),
      lik.lnsigma0,
      control = list(fnscale=-1)
    )
    # update parameters
    new.params = params
    new.params$sigma0 <- exp(res$par[1])
    return(new.params)
  }

  #' routine for updating alpha and sigma0 w.r.t marginal likelihood
  update.alpha.sigma0.Q = function(params){
    # function in terms of log alpha, log sigma0
    lik.lnalpha.lntheta = function(par){
      lnalpha = par[1]
      lnsigma0 = par[2]
      return(Q(params, alpha=exp(lnalpha), sigma0=exp(lnsigma0)))
    }
    # optimize
    res = optim(
      c(-1, -1),
      lik.lnalpha.lntheta,
      control = list(fnscale=-1)
    )
    # update parameters
    new.params = params
    new.params$alpha <- exp(res$par[1])
    new.params$sigma0 <- exp(res$par[2])
    return(new.params)
  }

  #' routine for updating alpha w.r.t marginal likelihood
  update.alpha.Q = function(params){
    # function in terms of log alpha, log sigma0
    lik.alpha= function(par){
      alpha = par[1]
      return(Q(params, alpha=alpha))
    }
    # optimize
    new.alpha = optimize(lik.alpha, c(0, 1), tol=1e-4, maximum=T)$maximum

    # update parameters
    new.params = params
    new.params$alpha <- new.alpha
    return(new.params)
  }
  
  #' routine for updating sigma0 w.r.t marginal likelihood
  update.sigma0.Q = function(params){
    # function in terms of log sigma0, log sigma0
    lik.sigma0= function(par){
      sigma0 = par[1]
      return(Q(params, sigma0=sigma0))
    }
    new.sigma0 = optimize(lik.sigma0, c(0, 100), tol=1e-4, maximum=T)$maximum
    # update parameters
    new.params = params
    new.params$sigma0 <- new.sigma0
    return(new.params)
  }
  
  #' full em routine
  #' optimize marginal: if TRUE, optimize alpha and sigma0 w.r.t marginal likelihood
  expectation.maximization = function(params.init, n.iter=50, tol=1e-4, update.alpha=TRUE, update.sigma0=TRUE){
    params.fit <- params.init
    params.fit$lik.history = -Inf
    pb <- progress::progress_bar$new(total = n.iter)
    
    for(i in 1:n.iter){
      pb$tick()
      params.fit <- update.responsibilities(params.fit)
      params.fit <- update.theta(params.fit)
      if(update.alpha){
        params.fit <- update.alpha.Q(params.fit)
      }
      if(update.sigma0){
        params.fit <- update.sigma0.Q(params.fit)
      }
      params.fit$lik.history = c(params.fit$lik.history, likelihood(params.fit))
      if(abs(diff(tail(params.fit$lik.history, 2))) < tol){
        params.fit$lik.history = tail(params.fit$lik.history, -1)
        break
      }
    }
    return(params.fit)
  }
  
  return(list(
    update.responsibilities=update.responsibilities,
    theta.obj=theta.obj.factory,
    theta.grad=theta.grad.factory,
    update.theta=update.theta.Q,
    update.alpha = update.alpha.Q,
    update.sigma0 = update.sigma0.Q,
    expectation.maximiztion = expectation.maximization,
    likelihood=likelihood,
    Q=Q
    )
  )
}
