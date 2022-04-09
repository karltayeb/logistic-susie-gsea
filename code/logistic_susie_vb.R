### VB LOWER BOUND SUSIE BELOW

### NOTE: This is basic code, and I did not attempt to mirror the level of numerical sophistication in the susie functions
### If this idea is worth pursuing further, then this code can be improved

calc_Q = function(X, Sigma2, Mu, Alpha, Z, delta) {
  # X is data matrix
  # Sigma2[j, l] is the posterior variance for b_l when entry j is selected, p x L
  # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
  # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
  # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)
  # delta is current estimate for effects of Z variables
  ASU2 = Alpha * (Sigma2 + Mu^2) # [j, l] = alpha[j, l] * (Sigma2[j, l] + Mu[j, l]^2)
  Q = rowSums(X^2 %*% ASU2) # start w/ sum_l E_(ql)[(x_i' b_l)^2]

  # now add 2 sum_l sum_{k>l} (x_i' b_l_post)(x_i' b_k_post)
  b_post_mat = Mu * Alpha # each column is posterior mean for b_l, p x L
  X_b_post = X %*% b_post_mat # [i, l] = x_i' b_l_post
  Q = Q + rowSums(X_b_post)^2 - rowSums(X_b_post^2)

  # now, add other terms with Z and delta
  Q = Q + as.numeric(2 * (X %*% rowSums(b_post_mat)) * (Z %*% delta)) + as.numeric((Z %*% delta)^2)
  return(Q)
}

g = function(x) { # ilogit function
  1 / (1 + exp(-x))
}


update_xi = function(X, Sigma2, Mu, Alpha, Z, delta) {
  # X is data matrix
  # Sigma2[j, l] is the posterior variance for b_l when entry j is selected, p x L
  # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
  # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
  # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)
  # delta is current estimate for effects of Z variables

  Q = calc_Q(X, Sigma2, Mu, Alpha, Z, delta)

  xi = sqrt(Q)

  return(xi)

}


update_b_l = function(X, y, xi, prior_weights, V, Sigma2, Mu, Alpha, l, Z, delta) {
  # X is data matrix
  # y is binary response
  # xi is lower-bound approximation parameters
  # prior_weights is prior probabilities for selecting j (p-vector)
  # V is prior variance of current effect (scalar)
  # Sigma2[j, l] is the posterior variance for b_l when entry j is selected, p x L
  # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
  # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
  # l is index to update
  # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)
  # delta is current estimate for effects of Z variables

  b_post_mat = Mu * Alpha # each column is posterior mean for b_l, p x L
  b_post_not_l = rowSums(as.matrix(b_post_mat[, -l], nrow = nrow(Mu))) # posterior, sum_(k != l) b_k_post
  g_xi = g(xi) # vector of g(xi_i), pre-compute once
  g_xi_5_xi = ((g_xi - .5) / xi) # vector of (g(xi_i) - .5) / xi_i, pre-compute once
  g_xi_5_xi[xi == 0] = .25 # case of 0/0 (e.g. x_i is all 0), use L'Hopital

  common_denoms = (1 / V) + (t(X^2) %*% g_xi_5_xi) # appears many times, compute once, posterior precisions

  # update Alpha[, l]
  nums = t(X) %*% (y - .5 -  (g_xi_5_xi * ((X %*% b_post_not_l) + (Z %*% delta)))) # numerator in exp
  Alpha[, l] = as.numeric(log(prior_weights) + (nums^2 / (2*common_denoms)) - (1/2)*log(common_denoms))
  Alpha[, l] = exp(Alpha[, l] - max(Alpha[, l])) # remove max for stability, everything still proportional
  Alpha[, l] = Alpha[, l] / sum(Alpha[, l]) # normalize, sum to 1

  # update Mu[, l]
  Mu[, l] = as.numeric(nums / common_denoms)

  # update Sigma[, l]
  Sigma2[, l] = as.numeric(1 / common_denoms)

  return(list(Sigma2 = Sigma2, Mu = Mu, Alpha = Alpha))

}


update_delta = function(X, y, xi, Mu, Alpha, Z) {
  # X is data matrix
  # y is binary response
  # xi is lower-bound approximation parameters
  # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
  # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
  # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)

  g_xi = g(xi) # vector of g(xi_i), pre-compute once
  g_xi_5_xi = ((g_xi - .5) / xi) # vector of (g(xi_i) - .5) / xi_i, pre-compute once
  g_xi_5_xi[xi == 0] = .25

  #D = Matrix::Diagonal(x = g_xi_5_xi / 2)
  D = Matrix::Diagonal(x = g_xi_5_xi)
  ZtDZ = t(Z) %*% D %*% Z
  DXB = D %*% X %*% rowSums(Mu * Alpha)
  RHS = t(Z) %*% (y - .5 - DXB) # RHS of system LL'delta = RHS

  # NOTE: the below system can/should be solved w/ Cholesky and forward/backward substitution (but if Z has small # of columns, shouldn't matter)
  #Lt = chol(ZtDZ)
  #u = forwardsolve(t(Lt), RHS)
  #delta = backsolve(Lt, u)
  ## COULD USE chol2inv!!!!
  #delta = solve(ZtDZ, RHS)
  delta = solve(as.matrix(ZtDZ), as.matrix(RHS))

  return(as.numeric(delta))

}

update_all = function(X,
                      y,
                      xi,
                      prior_weights,
                      V,
                      Sigma2,
                      Mu,
                      Alpha,
                      Z,
                      delta,
                      estimate_prior_variance,
                      share_prior_variance,
                      intercept) {
  # X is data matrix
  # y is binary response
  # xi is lower-bound approximation parameters
  # prior_weights is prior probabilities for selecting j (p-vector)
  # V is prior variance (scalar if the same for all L effects, or vector of length L if different)
  # Sigma2[j, l] is the posterior variance for b_l when entry j is selected, p x L
  # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
  # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
  # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)
  # delta is current estimate for effects of Z variables
  # estimate_prior_variance is logical for if prior variance, V, should be estimated

  L = ncol(Mu)

  # update delta
  if (intercept) { # if covariates and/or intercept
    delta = update_delta(X, y, xi, Mu, Alpha, Z)
  }

  if (length(V) == L) {
    prior_vars = V
  } else if (length(V) == 1) {
    prior_vars = rep(V, L)
    V = prior_vars
  } else {
    stop("Argument 'V' must be of length either 1 or L")
  }

  # now, iterate over l = 1:L
  for (l in 1:L) {
    res_l = update_b_l(X, y, xi, prior_weights, prior_vars[l], Sigma2, Mu, Alpha, l, Z, delta)
    Sigma2 = res_l$Sigma2
    Mu = res_l$Mu
    Alpha = res_l$Alpha
  }

  # now, update xi
  xi = update_xi(X, Sigma2, Mu, Alpha, Z, delta)

  if (estimate_prior_variance == TRUE) { # if estimating prior variance
    ASU2 = Alpha * (Sigma2 + Mu^2) # [j, l] = alpha[j, l] * (Sigma2[j, l] + Mu[j, l]^2)
    if (share_prior_variance) { # if common prior variance for all effects
      V = sum(ASU2) / ncol(Mu) # ncol(Mu) = L = sum_l sum_j alpha_jl
    } else { # if different for all effects (already know length is L, checked a few lines before)
      V = colSums(ASU2)
    }
  }


  return(list(Sigma2 = Sigma2, Mu = Mu, Alpha = Alpha, delta = delta, xi = xi, V = V))

}

# y is vector of binary response (n x 1, can be sparse)
# X is matrix of variables (n x p, can be sparse)
# L is number of non-zero effects, positive integer
# V is prior variance on non-zero effect size, scalar or L-vector
# prior_weights is a vector of prior inclusion probabilities (p x 1)
# intercept is a logical of if the intercept should be fitted (default to TRUE)
# Z is a vector of covariates to be controlled for (non-penalized effect estimates, n x q, can be sparse).
# NOTE: Z should NOT include an intercept column
# estimate_prior_variance is a logical if prior variance, V, should be estimated. If "TRUE", supplied value of V is the initial value
# tol is the convergence criterior measuring the change in the ELBO
# maxit is the maximum number of iterations
logistic.susie = function(
  X, y, L = 10, V = 1, prior_weights = NULL,
  init.intercept = NULL, intercept = TRUE, Z = NULL,
  estimate_prior_variance = TRUE, share_prior_variance = FALSE,
  standardize=TRUE,
  tol = 1e-3, maxit = 1000, verbose=0) {

  res <- logistic.susie.init(
    X, y, L, V, prior_weights, init.intercept, intercept, Z,
    estimate_prior_variance, share_prior_variance, standardize)
  if(verbose > 0){
    pb <- progress::progress_bar$new(
      format = "[:bar] :current/:total (:percent)", total=maxit
    )
    pb$tick(0)
  }
  converged <- FALSE
  for (i in 1:maxit){
    if(verbose>0){pb$tick()}
    res <- logistic.susie.iteration(res)
    if(diff(tail(res$elbo, 2)) < tol){
      converged <- TRUE
      if(verbose){message('\nconverged')}
      break
    }
  }
  if(!converged){
    warning(paste0('did not converge after ', maxit, ' iterations'))
  }
  res <- logistic.susie.wrapup(res)
  res$converged <- TRUE
  return(res)
}

logistic.susie.init = function(
  X, y, L = 10, V = 1, prior_weights = NULL,
  init.intercept = NULL, intercept = TRUE, Z = NULL,
  estimate_prior_variance = TRUE, share_prior_variance = FALSE,
  standardize=TRUE) {

  p = ncol(X)
  n = nrow(X)
  se = rep(1, p)
  if(standardize){
    se <- sqrt(sparseMatrixStats::colVars(X))
    se[is.na(se)] <- 1
    se[se < 1e-5] <- 1e-5  # clip small values for stability
    X <- wordspace::scaleMargins(X, cols = 1/se)  # TODO: weird dependency, reimpliment
  }
  else{
    se <- rep(1, p)
  }

  if (is.null(Z)) {
    Z = matrix(1, nrow = n, ncol = 1)
  } else {
    col_variances = apply(Z, MARGIN = 2, var)
    if (any(col_variances == 0)) { # is constant column in Z matrix
      stop("Matrix 'Z' cannot have a constant column")
    }
    Z = cbind(matrix(1, nrow = n, ncol = 1), Z) # add intercept column
  }
  if (is.null(prior_weights)) {
    prior_weights = rep(1 / p, p)
  }

  # place to store posterior info for each l = 1, ..., L
  # initialize: could think of something better
  #delta = glm(y ~ Z - 1, family = "binomial")$coef # initialize to regular GLM solution
  # TODO: what if we have fixed intercept but also covariates?
  if (!is.null(init.intercept)) { # if init.intercept is specified
    delta = init.intercept
  } else { # initialize to regular GLM solution w/ just Z (no X)
    delta = glm(as.numeric(y) ~ Z - 1, family = "binomial")$coef
  }
  Alpha = matrix(prior_weights, nrow = p, ncol = L)
  Mu = matrix(0, nrow = p, ncol = L)
  Sigma2 = matrix(V, nrow = p, ncol = L, byrow = T)
  xi = update_xi(X, Sigma2, Mu, Alpha, Z, delta)

  post_info = list(Sigma2 = Sigma2, Mu = Mu, Alpha = Alpha, delta = delta, xi = xi, V = V)
  dat <- list(
    X=X, X.scale=se, y=y, Z=Z, prior_weights=prior_weights,
    intercept=intercept,
    estimate_prior_variance=estimate_prior_variance,
    share_prior_variance=share_prior_variance
  )
  elbo <- calc_ELBO(
    dat$y, dat$X, post_info$Alpha, post_info$Mu, post_info$Sigma2, post_info$V,
    dat$prior_weights, dat$Z, post_info$delta, post_info$xi
  )
  res = list(dat=dat, post_info=post_info, elbo=elbo)
  return(res)
}

logistic.susie.iteration = function(res, compute_elbo=T){
  dat <- res$dat
  post_info <- res$post_info
  post_info  <- update_all(
    dat$X, dat$y, post_info$xi, dat$prior_weights,
    post_info$V, post_info$Sigma2, post_info$Mu,
    post_info$Alpha, dat$Z, post_info$delta,
    dat$estimate_prior_variance, dat$share_prior_variance, dat$intercept)
  res$post_info <- post_info
  if(compute_elbo){
    elbo <- calc_ELBO(
      dat$y, dat$X, post_info$Alpha, post_info$Mu, post_info$Sigma2, post_info$V,
      dat$prior_weights, dat$Z, post_info$delta, post_info$xi
    )
    res$elbo <- c(res$elbo, elbo)
  }
  return(res)
}

logistic.susie.wrapup = function(res){
  # TODO put back into original X scale
  X <- res$dat$X
  post_info <- res$post_info
  res = list(
    alpha = t(post_info$Alpha),
    mu = t(post_info$Mu),
    mu2 = t(post_info$Sigma2 + post_info$Mu^2),
    elbo = res$elbo,
    dat = res$dat,
    post_info = res$post_info
  )
  class(res) <- 'susie'
  colnames(res$alpha) <- colnames(X)
  colnames(res$mu) <- colnames(X)
  res$pip <- susieR::susie_get_pip(res)
  names(res$pip) <- colnames(X)
  res$sets <- susieR::susie_get_cs(res, X=X)
  res$intercept <- res$post_info$delta[1]
  return(res)
}

# calculate the variational lower bound
# CAREFUL: Not sure what to do when Alpha[j, l] = 0 (we get 0*ln(0)). I will set this to 0
# Note: This expression is only valid when xi has been updated to be sqrt(Q), where Q is the
# expectation of the square of the linear predictor under the variational posterior (what we update xi to nomrally)
calc_ELBO = function(y, X, Alpha, Mu, Sigma2, V, prior_weights, Z, delta, xi) {
  p = nrow(Mu)
  L = ncol(Mu)
  P = matrix(prior_weights, nrow = p, ncol = L)
  b_post = rowSums(Alpha * Mu)

  V = matrix(V, nrow = p, ncol = L, byrow = T) # for proper arithmetic, either if V is a scalar or a vector of length L

  expected_log_lik = sum(log(g(xi)) + (y - .5) * as.numeric(((X %*% b_post) + (Z %*% delta))) - (xi / 2))
  KL_div_vb_prior = Alpha * (log(Alpha) - log(P) + (log(V) / 2) - (log(Sigma2) / 2) - .5 + ((Sigma2 + Mu^2) / (2 * V)))
  KL_div_vb_prior[Alpha == 0] = 0
  KL_div_vb_prior = sum(KL_div_vb_prior)

  if (KL_div_vb_prior < 0) { # to diagnose any issues
    warning("KL Divergence calculated to be < 0")
  }

  # TERMS BELOW WRONG B/C FORGOT MULTIPLICIATIVE FACTOR OF alpha_jl
  #KL_div_vb_prior = sum(log(Alpha)) - L*sum(log(prior_weights)) + (L*p*log(V) / 2) - sum(log(Sigma2) / 2) - (L*p / 2) + (sum(Sigma2 + Mu^2) / (2*V))
  #KL_div_vb_prior = sum(log(Alpha)) - sum(log(Sigma2) / 2) + (sum(Sigma2 + Mu^2) / (2*V)) # up to a constant

  ELBO = expected_log_lik - KL_div_vb_prior
  return(ELBO)
}
