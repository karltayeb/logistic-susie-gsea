#' compute TP, FP, TN, FN for decision (prediction > thresh)
#' if flip = T, decision on (prediction < thresh)
score_thresh_single <- function(prediction, label, thresh=0.9, flip = F) {
  if (flip) {
    prediction = 1 - prediction
  }
  unlist(list(
    TP = sum((prediction >= thresh) & (label != 0)),
    FP = sum((prediction >= thresh) & (label == 0)),
    TN = sum((prediction < thresh) & (label == 0)),
    FN = sum((prediction < thresh) & (label != 0)),
    thresh=thresh
  ))
}

THRESH = c(seq(0, 0.98, 0.01), 1 - 10^(-(2:20)))

score_phyper = function(fit, sim){
  map_dfr(THRESH, ~score_thresh_single(fit$pHypergeometric, sim$active, thresh=.x, flip=T))
}

score_pip = function(fit, sim){
  map_dfr(THRESH, ~score_thresh_single(fit$pip, sim$active, thresh=.x, flip=F))
}

score_cvglmnet = function(fit, sim){
  nonzero <- tail(as.integer(coef(fit, s = "lambda.min")[,1] != 0), -1)
  map_dfr(THRESH, ~score_thresh_single(nonzero, sim$active, thresh=.x, flip=F))
}


glm2list = function(fit, s){
  tail(as.integer(coef(fit, s=s)[,1] != 0), -1)
}

score_glm_single = function(fit, sim, s){
  score_thresh_single(glm2list(fit, s=s), sim$active, s, flip=F)
}

score_glm = function(fit, sim){
  map_dfr(fit$lambda, ~score_glm_single(fit, sim, .x))
}

score_glmnet = function(fit, sim){

  x <- which(sim10$active == 1)
  lambdas <- fit$lambda
  map_dfr(lambdas, ~score_thresh_single(
    score_thresh_single(ddd(fit, s=.x), sim10$active, .x, flip=T)))

  coef(fit, lambdas[3])[, 1]
  plot(fit$lambda)
}
