fit_affinity_prop_ora = function(ora, X, y){
  sig.ora  <- ora %>% filter(pHypergeometric < 1e-2)
  X.sig <- X[, sig.ora$geneSet]
  X.y <- X.sig[(y == 1),]

  idsInSet <- X.y %>% gseasusie:::X2geneSet() %>%
    split(., .[,'geneSet']) %>%
    map(~.x$gene)

  minusLogP <- -log(ora[match(names(idsInSet), ora$geneSet),]$pHypergeometric)
  minusLogP[minusLogP == Inf] <- -log(.Machine$double.eps)
  apRes <- affinityPropagation(idsInSet, minusLogP)
  return(apRes)
}

fit_weighted_set_cover_ora = function(ora, X, y){
  sig.ora  <- ora %>% filter(pHypergeometric < 1e-2)
  X.sig <- X[, sig.ora$geneSet]
  X.y <- X.sig[(y == 1),]

  idsInSet <- X.y %>% gseasusie:::X2geneSet() %>%
    split(., .[,'geneSet']) %>%
    map(~.x$gene)

  minusLogP <- -log(ora[match(names(idsInSet), ora$geneSet),]$pHypergeometric)
  minusLogP[minusLogP == Inf] <- -log(.Machine$double.eps)
  wscRes <- weightedSetCover(idsInSet, 1 / minusLogP, 10, 1)
  return(wscRes)
}

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

  # other methods for summarizing ora results
  # apres <- fit_affinity_prop_ora(ora, X, y)
  # wscres <- fit_weighted_set_cover_ora(ora, X, y)

  #linear_fit <- susieR::susie(X, y)
  marginal_regression <- gseasusie::fit_marginal_regression_jax(X, y)
  residual_regression <- gseasusie::fit_residual_regression_jax(X, y, logistic.fit)

  res <- tibble(
    y = list(y),
    fit=list(logistic.fit),
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

enrichment_target_factory <- function(prefix, data_loader, .genesets, .propgenes, .experiments=NULL){
  # load data target
  data_target_name <- paste0(prefix, '_data')
  data_target_sym <- rlang::sym(data_target_name)
  data_target <- tar_target_raw(data_target_name, substitute(data_loader()))

  # meta target
  meta_target_name <- paste0(prefix, '_meta')
  meta_target_sym <- rlang::sym(meta_target_name)
  by <- c('db', 'ptop')
  tidy_eval <- targets::tar_option_get("tidy_eval")
  command <- substitute(eval(rlang::expr(crossing(
    db=.genesets,
    experiment= if(is.null(.experiments)) names(data_target_sym) else intersect(names(data_target_sym), .experiments) ,
    ptop = .propgenes)
  )))
  meta_target_command <- tarchetypes:::tar_group_by_command(command, by, tidy_eval)
  meta_target <- tar_target_raw(
    name = meta_target_name,
    command = meta_target_command,
    iteration='group'
  )

  # fit target
  fit_target_name <- paste0(prefix, '_enrichment_fit')
  fit_target_sym <- rlang::sym(fit_target_name)
  fit_target_command <- substitute(eval(rlang::expr(
    meta_target_sym %>%
      rowwise() %>%
      mutate(
        enrichment = list(run_enrichment(db, experiment, ptop, genesets, sc_pbmc_deseq_data))
      ) %>%
      ungroup() %>%
      unnest(enrichment)
  )))
  fit_target_pattern = substitute(map(meta_target_sym))
  fit_target <- tar_target_raw(
    name = fit_target_name,
    command = fit_target_command,
    pattern = fit_target_pattern
  )

  # summary target
  summary_target_name <- paste0(prefix, '_enrichment_summary')
  summary_target_command <- substitute(rlang::expr(
    fit_target_sym %>%
      rowwise() %>%
      mutate(
        enrichment_summary = list(summarise_enrichment(
          fit, ora, marginal_reg, residual_reg, genesets[[db]]$geneSetDes))) %>%
      select(db, experiment, ptop, enrichment_summary)
  ))
  summary_target_pattern = substitute(map(fit_target_sym))
  summary_target <- tar_target_raw(
    name = summary_target_name,
    command = summary_target_command,
    pattern = summary_target_pattern
  )
  return(list(data_target, meta_target, fit_target, summary_target))
}
