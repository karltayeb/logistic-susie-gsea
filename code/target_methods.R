## Methods =====================================================================
# prepare a "methods" tibble with all the parameters for each method
# to add a method add a line to methods
# method_function must accept
#   - X a gene x pathway matrix
#   - y a binary gene list membership

source('code/fit_susie.R')
source('code/fit_mr_ash.R')
source('code/fit_baselines.R')

## mr ash args
mr.ash.sa2 <- (2^(0.2*(seq(1, 20) - 1)) - 1)^2
mr.ash.fixed.args <- list(
  method_q='sigma_indep_q',
  method_g='caisa',
  sa2=mr.ash.sa2,
  update.sigma2 = F
)

mr.ash.args <- list(
  method_q='sigma_indep_q',
  method_g='caisa',
  sa2=mr.ash.sa2,
  update.sigma2 = T
)

## make methods tibble
method_bank <- tribble(
  ~method.name, ~method_function, ~method_args, ~prediction_col,
  'fet', 'fit.fishers.exact.test', NULL, 'one_minus_bh_p',
  'susie.1', 'fit.susie', list(L=1, estimate_residual_variance = T), 'pip',
  'susie.10', 'fit.susie', list(L=10, estimate_residual_variance = T), 'pip',
  'susie.20', 'fit.susie', list(L=20, estimate_residual_variance = T), 'pip',
  'susie.30', 'fit.susie', list(L=30, estimate_residual_variance = T), 'pip',

  'logistic.susie.10', 'fit.logistic.susie', list(L=10, estimate_prior_variance = T), 'pip',


  'logistic.susie.veb.boost.1', 'fit.logistic.susie.veb.boost', list(k=1), 'pip',
  'logistic.susie.veb.boost.10', 'fit.logistic.susie.veb.boost', list(k=10), 'pip',
  'logistic.susie.veb.boost.20', 'fit.logistic.susie.veb.boost', list(k=20), 'pip',
  'logistic.susie.veb.boost.30', 'fit.logistic.susie.veb.boost', list(k=30), 'pip',

  'susie.veb.boost.1', 'fit.susie.veb.boost', list(k=1), 'pip',
  'susie.veb.boost.10', 'fit.susie.veb.boost', list(k=10), 'pip',
  'susie.veb.boost.20', 'fit.susie.veb.boost', list(k=20), 'pip',

  'mr.ash.lasso', 'fit.mr.ash.lasso', mr.ash.args, 'pip',

  'lasso.gaussian', 'fit.lasso', NULL, 'selected',
  'elastic.net.gaussian', 'fit.elastic.net', NULL, 'selected',

  'lasso.binomial', 'fit.lasso.binomial', NULL, 'selected',
  'elastic.net.binomial', 'fit.elastic.net.binomial', NULL, 'selected')
