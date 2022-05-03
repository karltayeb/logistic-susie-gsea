## Simulation Factory

# Creates a target with name `{prefix}.spec`
# which returns a data frame to be fed into sim_factory
make_sim_spec <- function(prefix, params, batches=1, reps=1){
  target_spec_name <- paste0(prefix, '.spec')
  target_spec_sym <- sym(target_spec_name)

  spec <- tidyr::crossing(params, expand.grid(tar_rep=c(1:reps), tar_batch=c(1:batches)))
  spec <- enexpr(spec)

  target_expr <- expr(tar_group_by(
    !!target_spec_sym,
    as.data.frame(!!spec),
    background.logit, active.logit, sim_function, tar_batch))
  eval(target_expr)
}

#' take a simulation spec dataframe `spec` and generate the simulations using
#' gene set matrix `X`
#' returns a target with name `{prefix}.sim`
#' which is a tibble with spec info + simulation
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

#' Fit each simulation to each method in the method_spec
#' @param prefix  a prefix for the final target `{prefix}.fit`
#' @param X a gene set matrix
#' @param target_sim output of `sim_factory`
#' @param method_spec a method spec tibble
fit_factory_w_args <- function(prefix, X, target_sim, method_spec) {
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

#' Fit each simulation to each method in the method_spec
#' NOTE: method_args is just a label! The method_function only accepts X, y
#' @param prefix  a prefix for the final target `{prefix}.fit`
#' @param X a gene set matrix
#' @param target_sim output of `sim_factory`
#' @param method_spec a method spec tibble
fit_factory = function(prefix, X, target_sim, method_spec) {
  ts <- enexpr(target_sim)
  XX <- enexpr(X)

  method_spec <-
    method_spec %>%
    mutate(
      method_arg_names = map(method_args, names),
      idx = row_number()
    )

  MS <- method_spec

  target_name <- paste0(prefix, '.fit')
  tar_map(
    values = method_spec, names = method_name,
    tar_target_raw(
      target_name,
      expr(
        !!ts %>% rowwise() %>% mutate(
          method_name = method_name,
          method_args = if_else(
            length(method_args) > 0,
            list(list2namedlist(method_args, method_arg_names)),
            list(NULL)),
          fit = list(method_function(!!XX, sim$Y[,1]))
        )
      ),
      pattern = expr(map(!!ts))
    )
  )
}


## Model scoring factory

#' generates an expression for target `{prefix}.score_{method_name}`
#' which calls `score_function(fit, sim)`
generate_score_target <- function(prefix, method_name, score_function){
  # name the target
  target_score_name <- paste0(prefix, '.score_', method_name)
  target_score_name <- expr(!!target_score_name)

  score_function <- expr(!!score_function)

  # get the fit target
  target_fit <- sym(paste0(prefix, '.fit_', method_name))
  target_fit <- expr(!!target_fit)
  target_pattern = expr(map(!!target_fit))

  score_expression <- expr(
    tar_target_raw(
      !!target_score_name,
      expr(
        !!target_fit %>%
          rowwise() %>%
          mutate(scores = list((!!score_function)(fit, sim)))
      ), pattern=expr(!!target_pattern))
  )
  return(score_expression)
}

#' generate targets for all the methods in `methods`
score_factory <- function(prefix, methods){
  return(lapply(tar_eval(generate_score_target(prefix, method_name, score_function), values = methods), eval))
}


foo = function(a){
  a <- quote(!!a)
  return(a)
}
foo(score_ora)
