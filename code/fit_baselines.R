## fit baseline regression models & Fisher's exact test
## to include in pipeline register @ target_components/methods.R
library(glmnet)
library(clustermq)

cv.glmnet.summary <- function(X, y, ...) {
  cv.fit <- cv.glmnet(X, y, ...)
  beta <- coef(cv.fit, s='lambda.min')[, 1]
  selected <- as.numeric(beta != 0)
  names(selected) <- names(beta)
  return(tibble(model=list(cv.fit), beta=list(tail(beta, -1)), selected=list(tail(selected, -1))))
}

fit.lasso <- function(X, y){
  cv.glmnet.summary(X, y, family='gaussian', alpha=1, lower.limits = 0)
}

fit.lasso.binomial <- function(X, y){
  cv.glmnet.summary(X, y, family='binomial', alpha=1, lower.limits = 0)
}

fit.elastic.net.binomial <- function(X, y){
  cv.glmnet.summary(X, y, family='binomial', alpha=0.5, lower.limits = 0)
}

fit.elastic.net <- function(X, y){
  cv.glmnet.summary(X, y, family='gaussian', alpha=0.5, lower.limits = 0)
}

fit.fishers.exact.test <- function(X, y, ...) {
  gene.list <- names(which(y == 1))
  all.genes <- rownames(X)
  
  overlap <- (y %*% X)[1,]
  gene.list.size <- length(gene.list)
  background.size <- length(all.genes) - gene.list.size
  gene.set.size <- (rep(1, dim(X)[1]) %*% X)[1,]
  
  result <- tibble(
    gene.set.name = colnames(X),
    overlap = overlap,
    gene.list.size = gene.list.size,
    background.size = background.size,
    gene.set.size = gene.set.size,
    
  ) %>% mutate(
    odds.ratio = (overlap / (gene.list.size - overlap)) / ((gene.set.size - overlap) / (background.size - gene.set.size + overlap)),
    pValue = phyper(overlap-1, gene.list.size, background.size, gene.set.size, lower.tail= FALSE),
    FDR = p.adjust(pValue, method = 'BH'),
    one_minus_bh_p = 1 - FDR,
    nl10p = -log10(pValue)) %>% chop(everything())
  return(result)
}