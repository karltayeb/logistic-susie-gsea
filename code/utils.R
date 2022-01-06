# load and row_bind all targets matching a pattern
tar_agg <- function(pattern, selector=dplyr::matches){
  load.env <- new.env()
  tar_load(selector(pattern), envir = load.env)
  result <- grep(pattern , names(load.env), value=TRUE)
  result <- do.call("list", mget(result, envir = load.env))
  result <- bind_rows(result)
}

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

prob2logit <- function(logit){
  return(log(p) - log(1-p))
}

alpha2cs <- function(alpha, target_coverage=0.95){
  as_tibble(as.data.frame(t(alpha))) %>%
    mutate(geneset = colnames(alpha)) %>%
    pivot_longer(-geneset, names_to = 'component', values_to = 'alpha') %>%
    group_by(component) %>%
    arrange(desc(alpha)) %>%
    mutate(coverage = lag(cumsum(alpha), default=0)) %>%
    filter(cumall(!(coverage >= target_coverage))) %>%
    summarise(cs = list(geneset), cs_size = length(geneset), coverage=sum(alpha))
}
