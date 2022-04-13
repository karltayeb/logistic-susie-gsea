#' Pepare gene set matrix and binary gene list
#' @param gs a list containing gene set matrix `X`
#'    with genes in `rownames` and gene set names in `colnames
#' @param dat a vector of statistics to threshold hold,
#'    names must match columns of `X`
#' @param threshold to binarize `dat`
prep_binary_data = function(gs, dat, thresh, .sign=c(1, -1)) {
  # get common genes as background
  gs.genes <- rownames(gs$X)
  
  # filter gene list
  y <- dat %>%
    filter(ENTREZID %in% gs.genes, !is.na(threshold.on)) %>%
    mutate(threshold.on = threshold.on * as.integer(sign(beta) == .sign))
    mutate(geneList = as.integer(threshold.on < thresh)) %>%
    select(ENTREZID, geneList) %>% 
    distinct(across(ENTREZID), .keep_all = T) %>%
    tibble2namedlist()

  # filter gene sets
  y.genes <- names(y)
  test.genes <- intersect(gs.genes, y.genes)
  X <- gs$X[test.genes,]
  bad.cols <- BiocGenerics::colSums(X)
  bad.cols <- (bad.cols == 0) | bad.cols == length(test.genes)
  X <- X[, which(!bad.cols)]
  
  # reorder genes in y to match X order
  y <- y[rownames(X)]
  return(list(y=y, X=X))
}

do_logistic_susie = function(experiment,
                             db,
                             thresh,
                             genesets,
                             data,
                             susie.args=NULL,
                             .sign=c(1, -1)) {
  cat(paste0('Fitting logistic susie...\n\tExperiment = ', experiment, '\n\tDatabase = ', db, '\n\tthresh = ', thresh))
  gs <- genesets[[db]]
  dat <- data[[experiment]]
  u <- prep_binary_data(gs, dat, thresh, .sign)  # subset to common genes
  
  if(is.null(susie.args)){  # default SuSiE args
    susie.args = list(
      L=10, init.intercept=0, verbose=1, maxit=500, standardize=FALSE)
  }
  vb.fit <- exec(logistic.susie, u$X, u$y, !!!susie.args)
  res = tibble(
    experiment=experiment,
    db=db,
    thresh=thresh,
    fit=list(vb.fit),
    susie.args = list(susie.args)
  )
  return(res)
}

do_logistic_susie_veb = function(experiment,
                                 db,
                                 thresh,
                                 genesets,
                                 data,
                                 susie.args=NULL,
                                 .sign=c(1,-1)){
  gs <- genesets[[db]]
  dat <- data[[experiment]]
  u <- prep_binary_data(gs, dat, thresh, .sign)  # subset to common genes
  
  if(is.null(susie.args)){  # default SuSiE args
    susie.args = list(L=10)
  }
  vb.fit <- exec(logistic.susie.veb.boost, u$X, u$y, !!!susie.args)
  res = tibble(
    experiment=experiment,
    db=db,
    thresh=thresh,
    fit=list(vb.fit)
  )
  return(res)
}

do_ora = function(experiment,
                  db,
                  thresh,
                  genesets,
                  data,
                  .sign=c(1,-1){
  
  gs <- genesets[[db]]
  dat <- data[[experiment]]
  u <- prep_binary_data(gs, dat, thresh, .sign)  # subset to common genes
  
  do_fet = function(overlap, geneListSize, geneSetSize, nGenes){
    ct <- matrix(c(
      overlap,
      geneListSize-overlap,
      geneSetSize-overlap,
      nGenes - geneSetSize - geneListSize + overlap), nr=2)
    return(fisher.test(ct)$p.value)
  }
  
  do_hyper = function(overlap, geneListSize, geneSetSize, nGenes){
    # genes in list = white balls, not in list = black balls
    # draw k balls w/o replacement where k is the size of the gene set
    # is the overlap we see with our gene set extreme/unlikely?
    return(phyper(
      max(0, overlap - 1),  # p(X >= overlap) = p(X > overlap - 1)
      geneListSize,
      nGenes - geneListSize,
      geneSetSize,
      lower.tail = FALSE))
  }
  
  ora <- tibble(
      geneSet = colnames(u$X),
      geneListSize = sum(u$y),
      geneSetSize = BiocGenerics::colSums(u$X),
      overlap = (u$y %*% u$X)[1,],
      nGenes = length(u$y),
      propInList = overlap / geneListSize,
      propInSet = overlap / geneSetSize,
      oddsRatio = (overlap / (geneListSize - overlap)) / (
        (geneSetSize - overlap) / (nGenes - geneSetSize + overlap))
    ) %>%
    rowwise() %>%
    mutate(
      pHypergeometric = do_hyper(overlap, geneListSize, geneSetSize, nGenes),
      pFishersExact = do_fet(overlap, geneListSize, geneSetSize, nGenes)
    ) %>% 
    ungroup()
  
  res = tibble(
    experiment=experiment,
    db=db,
    thresh=thresh,
    ora=list(ora)
  )
  return(res)
}

get_credible_set_summary = function(res){
  #' report top 50 elements in cs
  beta <- t(res$mu) %>%
    data.frame() %>%
    rownames_to_column(var='geneSet') %>%
    rename_with(~str_replace(., 'X', 'L')) %>%
    dplyr::rename(L1 = 2) %>%  # rename deals with L=1 case
    pivot_longer(starts_with('L'), names_to='component', values_to = 'beta') 
  se <- t(sqrt(res$mu2 - res$mu^2)) %>%
     data.frame() %>%
      rownames_to_column(var='geneSet') %>%
      rename_with(~str_replace(., 'X', 'L')) %>%
      dplyr::rename(L1 = 2) %>%  # rename deals with L=1 case
      pivot_longer(starts_with('L'), names_to='component', values_to = 'beta.se')
  
  credible.set.summary <- t(res$alpha) %>%
    data.frame() %>%
    rownames_to_column(var='geneSet') %>%
    rename_with(~str_replace(., 'X', 'L')) %>%
    dplyr::rename(L1 = 2) %>%  # rename deals with L=1 case
    pivot_longer(starts_with('L'), names_to='component', values_to = 'alpha') %>%
    left_join(beta) %>%
    left_join(se) %>%
    arrange(component, desc(alpha)) %>%
    dplyr::group_by(component) %>%
    filter(row_number() < 50) %>%
    mutate(alpha_rank = row_number(), cumalpha = c(0, head(cumsum(alpha), -1))) %>%
    mutate(in_cs = cumalpha < 0.95) %>%
    mutate(active_cs = component %in% names(res$sets$cs))
  return(credible.set.summary)
}

get_gene_set_summary = function(res){
  #' map each gene set to the component with top alpha
  #' report pip
  res$pip %>% 
    as_tibble(rownames='geneSet') %>%
    dplyr::rename(pip=value) %>%
    mutate(beta=colSums(res$alpha * res$mu))
}

pack_group = function(tbl){
    components <- tbl$component
    unique.components <- unique(components)
    start <- match(unique.components, components)
    end <- c(tail(start, -1) - 1, length(components))
    res <- tbl %>% select(-c(component)) %>% kbl()
    for(i in 1:length(unique.components)){
      res <- pack_rows(res, unique.components[i], start[i], end[i])
    }
    return(res)
}

#' Report credible set based summary of SuSiE
#' @param tbl output of \alias{get_credible_set_summary} to be formatted
#' @param target_coverage the coverage of the credible sets to be reported
#' @param max_coverage report SNPs that are not in the target_coverage c.s. up this value
#' @param max_sets maximum number of gene sets to report for a single credible set
#'    this is useful for large credible sets
#' @return A styled kable table suitable for rendering in an HTML document
report_susie_credible_sets = function(tbl,
                                      target_coverage=0.95,
                                      max_coverage=0.99,
                                      max_sets=10){
  require(kableExtra)
  tbl_filtered <-
    tbl %>%
    group_by(component) %>%
    arrange(component, desc(alpha)) %>%
    filter(cumalpha <= max_coverage, alpha_rank <= max_sets) %>%
    mutate(in_cs = cumalpha <= target_coverage) %>% ungroup()

  tbl_filtered %>%
    select(component, geneSet, description, alpha, beta, beta.se, pHypergeometric, pFishersExact, overlap, geneSetSize, oddsRatio) %>%
    dplyr::mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
    pack_group %>%
    column_spec(c(4), color=ifelse(tbl_filtered$beta >0, 'green', 'red')) %>%
    kableExtra::kable_styling()
}


make_susie_html_table = function(fits, db, experiment, thresh){
  require(htmltools)
  fits %>%
    pluck(db) %>%
    pluck(topic) %>%
    get_credible_set_summary() %>%
    filter(active_cs, in_cs) %>%
    report_susie_credible_sets(target_coverage = 0.95, max_coverage = 0.95)
}
