## score simulations

score.thresh.single <- function(prediction, label, thresh=0.9, flip = F) {
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

count.enriched.thresh <- function(prediction, thresh=0.9) {
  unlist(list(
    positive = sum((prediction > thresh)),
    negative = sum((prediction > thresh)),
    thresh=thresh
  ))
}

get.score <- function(prediction, label, thresholds=NULL, flip=F){
  if (is.null(thresholds)){
    thresholds <- unique(c(seq(0, 0.99, 0.01), c(0.99, 0.999, 0.9999, 0.99999, 1.0)))
  }
  scores <-do.call(
    'rbind',
    map(thresholds, ~ score.thresh.single(prediction, label, .x, flip=flip)))
  return(as_tibble(scores))
}


count.enriched <- function(prediction, thresholds=NULL){
  if (is.null(thresholds)){
    thresholds <- unique(c(seq(0, 0.99, 0.9), c(0.99, 0.999, 0.9999, 0.99999, 1.0)))
  }
  scores <-do.call(
    'rbind',
    map(thresholds, ~ count.enriched.thresh(prediction, .x, flip=flip)))
  return(as_tibble(scores))
}

score.credible.set <- function(fit, active){
  fit %>%
    rowwise() %>%
    mutate(
      cs = list(alpha2cs(alpha)),
      active = list(names(which(active != 0)))
    ) %>%
    unnest(cs) %>% rowwise() %>% mutate(
      cs_in_causal = list(intersect(active, cs)),
      n_causal = length(cs_in_causal),
      cs_size = length(cs)
    ) %>% ungroup() %>%
    select(component, cs, cs_size, cs_in_causal, n_causal, coverage) %>%
    chop(c(component, cs, cs_size, cs_in_causal, n_causal, coverage))
}

count.credible.set <- function(fit, vars, target_coverage=0.95){
  fit %>%
    rowwise() %>%
    mutate(
      cs = list(alpha2cs(alpha, target_coverage))
    ) %>%
    unnest(cs) %>% rowwise() %>% mutate(
      cs_size = length(cs)
    )  %>%
    ungroup() %>%
    group_by_at(all_of(c(vars))) %>%
    select(component, cs, cs_size, coverage) %>%
    chop(c(cs_size, coverage, cs, component)) %>%
    mutate(target_coverage = target_coverage)
}

active.pip.batch <- function(sim, fit, method='susie', col='pip'){
  flip  = (method %in% c('hypergeometric'))
  score <- left_join(fit, sim) %>%
    mutate(value = !!sym(col), metric = col, method=method) %>%
    select(sim.name, rep, batch, active, value, metric, batch, rep, method)
  return(score)
}


