#' compute TP, FP, TN, FN for decision (prediction > thresh)
#' if flip = T, decision on (prediction < thresh)
score_thresh_single <- function(prediction, label, thresh=0.9, flip = F) {
  if (flip) {
    prediction = 1 - prediction
  }
  unlist(list(
    TP = sum((prediction > thresh) & (label != 0)),
    FP = sum((prediction > thresh) & (label == 0)),
    TN = sum((prediction <= thresh) & (label == 0)),
    FN = sum((prediction <= thresh) & (label != 0)),
    thresh=thresh
  ))
}

THRESH = seq(0, 1, 0.01)

score_phyper = function(fit, sim){
  map_dfr(THRESH, ~score_thresh_single(fit$pHypergeometric, sim$active, thresh=.x, flip=T))
}

score_pip = function(fit, sim){
  map_dfr(THRESH, ~score_thresh_single(fit$pip, sim$active, thresh=.x, flip=F))
}
