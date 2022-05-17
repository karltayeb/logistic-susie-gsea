run_enrichment = function(db, experiment, ptop, genesets, data){
  message(paste0('Fitting logistic susie...',
                 '\n\tExperiment = ', experiment,
                 '\n\tDatabase = ', db,
                 '\n\tptop = ', ptop))

  bin.data <- gseasusie::prep_binary_data(genesets[[db]], data[[experiment]], ptop=ptop)
  X <- bin.data$X
  y <- bin.data$y

  logistic.fit <- gseasusie::fit_logistic_susie_veb_boost(X, y, L=20)
  ora <- gseasusie::fit_ora(X, y)

  # linear.fit <- susieR::susie(X, y)
  marginal_regression <- gseasusie::fit_marginal_regression_jax(X, y)
  residual_regression <- gseasusie::fit_residual_regression_jax(X, y, logistic.fit)

  res <- tibble(
    fit=list(logistic.fit),
    # linear.fit = list(linear.fit),
    ora=list(ora),
    marginal_reg = list(marginal_regression),
    residual_reg = list(residual_regression)
  )
  return(res)
}

summarise_enrichment = function(fit, ora, marginal_reg, residual_reg, geneSetDes){
  # put everything into one table
  gssum <- gseasusie:::get_gene_set_summary(fit)
  cssum <- gseasusie:::get_cs_summary_condensed(fit)

  reg <- left_join(
    marginal_reg,
    residual_reg,
    by="geneSet", suffix=c('_marginal', '_residual'))

  gene_set_summary <-
    gssum %>%
    dplyr::left_join(cssum) %>%
    dplyr::left_join(ora) %>%
    left_join(reg) %>%
    left_join(geneSetDes)
  return(gene_set_summary)
}
