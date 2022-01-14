## old
sim_factory_single <- function(name, FUN, params, batches=1, reps=1){
  # replicate simulation
  list(
    tar_rep_raw(
      name,
      expr(lift_dl(match.fun(!!FUN))(params)),
      batches = batches,
      reps = reps,
      iteration='group'
    ) %>%
      tar_hook_outer(
        matches(name),
        hook = as.tibble(.x)
      )
  )
}

## Simulation Factory
sim_factory_old <- function(prefix, X, params, batches=1, reps=1){
  X <- enexpr(X)
  target_name <- paste0(prefix, '.sim')
  tar_rep_raw(
    target_name,
    expr(
      !!params %>%
        rowwise() %>%
        mutate(sim = list(lift_dl(match.fun(sim_function))(c(list(X=!!X), params))))),
    batches = batches,
    reps = reps,
    iteration='group') %>%
    tar_hook_outer(
      matches(target_name),
      hook = as.tibble(.x)
    )
}

## Simulation Factory
# group by param setting so we don't need to repeat computation everytime
# we add to the simulation
make_sim_spec <- function(prefix, params, batches=1, reps=1){
  target_spec_name <- paste0(prefix, '.spec')
  target_spec_sym <- sym(target_spec_name)

  spec <- tidyr::crossing(params, expand.grid(tar_rep=c(1:reps), tar_batch=c(1:batches)))
  spec <- enexpr(spec)

  target_expr <- expr(tar_group_by(!!target_spec_sym, as.data.frame(!!spec), background.logit, active.logit, sim_function, tar_batch))
  eval(target_expr)
}

sim_factory <- function(prefix, X, spec){
  X <- enexpr(X)
  target_spec_name <- paste0(prefix, '.spec')
  target_spec_sym <- sym(target_spec_name)
  target_spec <- expr(target_spec_sym)
  spec <- enexpr(spec)
  target_name <- paste0(prefix, '.sim')

  tar_target_raw(
    target_name,
    expr(!!spec %>%rowwise() %>% mutate(
      sim = list(lift_dl(match.fun(sim_function))(c(list(X=!!X), params))))),
    pattern = expr(map(!!spec))
  )
}

list2namedlist <- function(l1, l2){
  names(l1) <- l2
  return(l1)
}

## Model fitting factory
fit_factory <- function(prefix, X, target_sim, method_spec) {
  ts <- enexpr(target_sim)
  XX <- enexpr(X)
  f <- expr(
    lift(match.fun(method_function))(sym(X), sim$Y[,1])
  )
  target_name <- paste0(prefix, '.fit')
  tar_map(
    values = method_spec, names = method.name,
    tar_target_raw(
      target_name,
      expr(
        !!ts %>% rowwise() %>% mutate(
          method = method.name,
          pf = list(force(purrr::partial(match.fun(method_function), X=!!XX, y=sim$Y[,1]))),
          fit = list(lift_dl(pf)(list2namedlist(method_args, ma_names)))
        ) %>% select(!pf)
      ),
      pattern = expr(map(!!ts))
    )
  )
}

## Model scoring factory
generate_score_target <- function(prefix, method.name, prediction_col){
  target_score_name <- paste0(prefix, '.score_', method.name)
  target_score_name <- expr(!!target_score_name)

  target_fit <-sym(paste0(prefix, '.fit_', method.name))
  target_fit <- expr(!!target_fit)

  prediction_col <- expr(!!prediction_col)
  method_name <- expr(!!method.name)

  target_pattern = expr(map(!!target_fit))

  score_expression <- expr(
    tar_target_raw(
      !!target_score_name,
      expr(!!target_fit %>% rowwise() %>% mutate(
        scores = list(get.score(fit[[!!prediction_col]][[1]], sim$active)),
        score.cs = if(grepl('susie', !!method_name)) {list(score.credible.set(fit, sim$active))} else {list(NULL)},
        method = !!method_name
      )), pattern=expr(!!target_pattern)
    )
  )
  return(score_expression)
}
score_factory <- function(prefix, methods){
  return(lapply(tar_eval(generate_score_target(prefix, method.name, prediction_col), values = methods), eval))
}
