library(Matrix)
source('code/utils.R')

simulate.constant.sim <- function(X, L=1, background.logit=-10, active.logit=10, rep=1, seed=NULL){
  # X is a gene x gene_set matrix indicating membership in gene set  n x p
  # L is the number of gene sets to have active
  # background_logit, active_logit logodds of observing gene in our gene list
  # reps number of replicates to simulate
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- dim(X)[1]
  p <- dim(X)[2]

  logit <- rep(background.logit, n)
  active_gene_sets <- sample(seq(1:p), L, replace = F)
  for (a in active_gene_sets){
    logit[X[, a] == 1] <- active.logit
  }
  active = rep(0, p)
  active[active_gene_sets] <- 1

  Y <- matrix(rbinom(n * rep, 1, logit2prob(logit)), n)
  Y <- Matrix(Y, sparse = T)

  names(active) <- colnames(X)
  names(logit) <- rownames(X)
  rownames(Y) <- rownames(X)

  sim <- list(active=active, logit=logit, Y=Y)
  return(sim)
}

simulate.constant.additive.sim <- function(X, L=1, background.logit=-10, active.logit=1.0, rep=1, seed=NULL){
  # X is a gene x gene_set matrix indicating membership in gene set  n x p
  # L is the number of gene sets to have active
  # background_logit, active_logit logodds of observing gene in our gene list
  # reps number of replicates to simulate
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- dim(X)[1]
  p <- dim(X)[2]

  active_gene_sets <- sample(seq(1:p), L, replace = F)
  effect_sizes <- rep(1, L) * active.logit

  if (L > 1){
    logit <- as.vector(background.logit + X[, active_gene_sets] %*% effect_sizes)
  } else{
    logit <- as.vector(background.logit + X[, active_gene_sets] * effect_sizes)
  }

  active = rep(0, p)
  active[active_gene_sets] <- effect_sizes

  Y <- matrix(rbinom(n * rep, 1, logit2prob(logit)), n)
  Y <- Matrix(Y, sparse = T)

  names(active) <- colnames(X)
  names(logit) <- rownames(X)
  rownames(Y) <- rownames(X)

  sim <- list(active=active, logit=logit, Y=Y)
  return(sim)
}

simulate.normal.additive.sim <- function(X, L=1, background.logit=-10, sigma2=1.0, rep=1, seed=NULL){
  # X is a gene x gene_set matrix indicating membership in gene set  n x p
  # L is the number of gene sets to have active
  # background_logit, active_logit logodds of observing gene in our gene list
  # reps number of replicates to simulate

  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- dim(X)[1]
  p <- dim(X)[2]

  active_gene_sets <- sample(seq(1:p), L, replace = F)
  effect_sizes <- abs(rnorm(L, sd = sqrt(sigma2)))

  print(paste(active_gene_sets, effect_sizes))
  logit <- as.vector(background.logit + X[, active_gene_sets] %*% effect_sizes)

  active = rep(0, p)
  active[active_gene_sets] <- effect_sizes

  Y <- matrix(rbinom(n * rep, 1, logit2prob(logit)), n)
  Y <- Matrix(Y, sparse = T)

  names(active) <- colnames(X)
  names(logit) <- rownames(X)
  rownames(Y) <- rownames(X)

  sim <- list(active=active, logit=logit, Y=Y)
  return(sim)
}

simulate.sim.data.rep <- function(X, L=1, background.logit=-10, active.logit=10, inner.rep=2, outer.rep=10, annotations=list()){
  sim.rep <- t(replicate(outer.rep, simulate.sim.data(
    X, L=L, background.logit=background.logit, active.logit=active.logit,
    rep=inner.rep, annotations =annotations)))
}

simulate.sim.data <- function(X, L=1, background.logit=-10, active.logit=10, rep=1, annotations=list()){
  # X is a gene x gene_set matrix indicating membership in gene set  n x p
  # L is the number of gene sets to have active
  # background_logit, active_logit logodds of observing gene in our gene list
  # reps number of replicates to simulate

  n <- dim(X)[1]
  p <- dim(X)[2]

  logit <- rep(background.logit, n)
  active_gene_sets <- sample(seq(1:p), L, replace = F)
  for (a in active_gene_sets){
    logit[X[, a] == 1] <- active.logit
  }

  Y <- matrix(rbinom(n*rep, 1, logit2prob(logit)), nrow = n, byrow = FALSE)
  sim <- list(active=active_gene_sets, logit=logit, Y=Y, X=X, rep=rep, background.logit=background.logit, active.logit=active.logit)
  sim <- c(sim, annotations)
  return(sim)
}

simulate.sim.data.rep <- function(X, L=1, background.logit=-10, effect.var=1, inner.rep=2, outer.rep=10, annotations=list()){
  sim.rep <- t(replicate(outer.rep, simulate.sim.data(
    X, L=L, background.logit=background.logit, active.logit=active.logit,
    rep=inner.rep, annotations =annotations)))
}
