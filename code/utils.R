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

prob2logit <- function(p){
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

convert_labels <- function(y){
  hs <- org.Hs.eg.db::org.Hs.eg.db
  gene_symbols <- names(y)
  symbol2entrez <- AnnotationDbi::select(hs, keys=gene_symbols, columns=c('ENTREZID', 'SYMBOL'), keytype = 'SYMBOL')
  symbol2entrez <- symbol2entrez[!duplicated(symbol2entrez$SYMBOL),]
  rownames(symbol2entrez) <- symbol2entrez$SYMBOL
  ysub <- y[names(y) %in% symbol2entrez$SYMBOL]
  names(ysub) <- symbol2entrez[names(ysub),]$ENTREZID
  return(ysub)
}

procrustes <- function(X, y){
  idx <- intersect(rownames(X), names(y))
  subX <- X[idx,]
  subX <- subX[,(colSums(subX) > 0 & colSums(subX) < dim(subX)[1])]
  return(list(X=subX, y=y[idx]))
}

make_gene_list <- function(z, z_threshold=2){
  y <- as.integer(abs(z) > z_threshold)
  names(y) <- names(z)
  return(y)
}

process_input <- function(X, y){
  if(length(intersect(rownames(X), colnames(y)))){
    message('No interseccting genes, converting to ENTREZ IDs')
    y <- convert_labels(y)
  }
  return(procrustes(X, y))
}


susie_plot2 <- function(fit, max_set_size=200, ...){
  gs <- colnames(fit$alpha[[1]])
  plot(fit$pip[[1]], cex=0.001, ylab='PIP', xlab='Gene Set', ...)

  to_plot <- fit$cs[[1]] %>%
    mutate(plot = cs_size < 20) %>%
    dplyr::select(plot) %>% purrr::pluck(1)
  n_cs <- sum(to_plot)
  cols <- RColorBrewer::brewer.pal(n_cs,'Set1')

  print(n_cs)
  color <- 1
  for (i in which(to_plot)){
    idx = (gs %in% fit$cs[[1]]$cs[[i]])
    print(which(idx))
    points(x=which(idx), y=fit$pip[[1]][which(idx)], col=cols[color], cex=2, pch=16)
    color <- color+1
  }
  points(fit$pip[[1]], col='black', pch=16, cex=0.5)
}
