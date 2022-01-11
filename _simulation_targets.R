library(tidyverse)
library(glue)
library(targets)
library(tarchetypes)
library(clustermq)

# callr cant find these when used in 'tar_option_set'
# but these imports would be better placed there
# Both installed via bioconductor so maybe that's a hint
require(GSEABenchmarkeR, quietly = T)
require(EnrichmentBrowser, quietly = T)

source('code/load_gene_sets.R')
source('code/simulate_gene_lists.R')
source('code/score_simulations.R')
#source('code/utils.R')
#source('code/plots.R')

options(clustermq.scheduler = "multiprocess")
tar_option_set(workspace_on_error = TRUE)

tar_option_set(packages = c(
  "tidyverse", "WebGestaltR",
  "Matrix", "susieR", "glmnet",
  "mr.ash.alpha", "VEB.Boost", "fastglm"))
tar_option_set(memory='transient', garbage_collection = T, error='workspace')

load.tar <- list(
  tar_target(X.gonr, load.webGestalt.X()),
  tar_target(X.gobp, load.gobp()),
  #tar_target(X.c2, load.msigdb.X()),
  tar_target(X, list(go_bp = X.gobp, go_no_redundant=X.gonr)),
  tar_target(genesets, tibble(name = c('go_bp', 'go_no_redundant')))
)


# Methods ======================================================================
# to add a method add a line to methods
# method_function must accept
#   - X a gene x pathway matrix
#   - y a binary gene list membership

source('target_components/methods.R')

# Simulations ==================================================================

# L1 sim =======================================================================

get_l1.sim.methods <- function(){
  method_bank %>% dplyr::filter(
    method.name %in% c(
      'fet',
      'susie.10',
      'logistic.susie.veb.boost.10',
      'lasso', 'elastic.net'
    )
  )
}

get_l1.sim.params <- function(){
  return(crossing(
    active.logit=c(-6, -5, -4, -3, -2),
    background.logit=c(-2, -1, 0, 1, 2)))
}

l1.sim.methods <- get_l1.sim.methods() %>%
  mutate(method_arg_names = map(method_args, names))

args2list <- function(X, y, args){
  args <- c(list(X=X, y=y), args)
  return(args)
}

list2namedlist <- function(l1, l2){
  names(l1) <- l2
  return(l1)
}

l1.sim.tar <- list(
  # get simulation parameters
  tar_target(
    l1.sim.spec,
    get_l1.sim.params()
  ),
  # replicate simulation
  tar_rep(
    l1.sim.rep,
    l1.sim.spec %>% rowwise() %>%
      mutate(sim = list(simulate.constant.sim(
        X.gonr, L=1, background.logit = -1, active.logit = 1))) %>%
      unnest_wider(sim),
    batches=1,
    reps=20
  ),
  tar_group_size(l1.sim, l1.sim.rep, 20),
  tar_map(
    values=l1.sim.methods, names=method.name,
    tar_target(
      l1.sim.fit,
      l1.sim %>%
        rowwise() %>%
        mutate(
          fit = list(do.call(get(method_function), args2list(X.gonr, Y[, 1], list2namedlist(method_args, method_arg_names))))
        ), pattern=map(l1.sim)
    ),
    tar_target(
      l1.sim.score,
      l1.sim.fit %>%
        rowwise() %>%
        mutate(
          scores = list(get.score(fit[[prediction_col]][[1]], active)),
          score.cs = if(grepl('susie', method.name)) {list(score.credible.set(fit, active))} else {list(NULL)}
        ),
      pattern=map(l1.sim.fit)
    )
  )
)

# Targets list =================================================================
# final list of targets

list(
  load.tar, l1.sim.tar
)


