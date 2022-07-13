#' a collection of functions for loading
#' standardized versions of the data
library(tidyverse)

generate_map2entrez = function(genes, from, to='ENTREZID'){
  hs <- org.Hs.eg.db::org.Hs.eg.db
  genes <- unique(genes)
  map2entrez <- AnnotationDbi::select(
    hs, keys=genes,
    columns=c(to, from),
    keytype = from)
  return(map2entrez)
}

#' add names to list, useful for piping
add_names = function(l, n){
  names(l) <- n
  return(l)
}


load_droplet = function() {
  load("data/de-droplet/de-droplet-noshrink.RData")
  gene_symbols <- toupper(unique(rownames(de_merged$postmean)))
  map2enrtrez <- generate_map2entrez(gene_symbols, 'SYMBOL')

  # Load the results of the DE analysis.
  extract_celltype = function(k){
    b  <- de_merged$postmean[, k]
    se <- with(de_merged,postmean/z)[, k]
    p <- 2 * pnorm(abs(de_merged$z[, k]), lower.tail=F)
    data.frame(
        beta = b,
        se = se,
        threshold.on = p) %>%
      rownames_to_column(var='SYMBOL') %>%
      mutate(SYMBOL=toupper(SYMBOL)) %>%
      left_join(map2enrtrez, by='SYMBOL') %>%
      relocate(ENTREZID, .after=SYMBOL)
  }

  data <- map(1:6, extract_celltype)
  names(data) <- colnames(de_merged$postmean)
  return(data)
}

# load different versions of the pbmc data
load_sc_pbmc_deseq2 = function() {
  load('data/pbmc-purified/deseq2-pbmc-purified.RData')
  genes <- unique(rownames(deseq[[1]]))
  map2entrez <- generate_map2entrez(genes, 'ENSEMBL')

  standardize_data = function(dat){
    dat %>%
      as.data.frame() %>%
      rownames_to_column(var='ENSEMBL') %>%
      left_join(map2entrez, by='ENSEMBL') %>%
      relocate(ENTREZID, .after=ENSEMBL) %>%
      mutate(  # set default columns
        beta = log2FoldChange,
        se = lfcSE,
        threshold.on = pvalue
      )
  }
  data <- map(deseq, standardize_data)
  return(data)
}

load_sc_pbmc_noshrink = function() {
  load("~/Research/logistic-susie-gsea/data/pbmc-purified/de-pbmc-purified-noshrink.RData")

  hs <- org.Hs.eg.db::org.Hs.eg.db
  names(de$postmean)
  gene_symbols <- unique(rownames(de$postmean))
  map2enrtrez <- generate_map2entrez(gene_symbols, 'ENSEMBL')
  add_names = function(l, n){
    names(l) <- n
    return(l)
  }

  # Load the results of the DE analysis.
  extract_celltype = function(k){
    b  <- de$postmean[, k]
    se <- with(de,postmean/z)[, k]
    p <- 2 * pnorm(abs(de$z[, k]), lower.tail=F)
    data.frame(
        beta = b,
        se = se,
        threshold.on = p) %>%
      rownames_to_column(var='ENSEMBL') %>%
      left_join(map2enrtrez, by='ENSEMBL') %>%
      relocate(ENTREZID, .after=ENSEMBL)
  }

  data <- map(1:6, extract_celltype)
  names(data) <- paste0('k', 1:6)
  return(data)
}

load_sc_pbmc_purified = function() {
  load("~/Research/logistic-susie-gsea/data/pbmc-purified/de-pbmc-purified-seed=1.RData")

  hs <- org.Hs.eg.db::org.Hs.eg.db
  names(de$postmean)
  gene_symbols <- unique(rownames(de$postmean))
  map2enrtrez <- generate_map2entrez(gene_symbols, 'ENSEMBL')
  add_names = function(l, n){
    names(l) <- n
    return(l)
  }

  # Load the results of the DE analysis.
  extract_celltype = function(k){
    b  <- de$postmean[, k]
    se <- with(de,postmean/z)[, k]
    p <- 2 * pnorm(abs(de$z[, k]), lower.tail=F)
    data.frame(
      beta = b,
      se = se,
      threshold.on = p) %>%
      rownames_to_column(var='ENSEMBL') %>%
      left_join(map2enrtrez, by='ENSEMBL') %>%
      relocate(ENTREZID, .after=ENSEMBL)
  }

  data <- map(1:6, extract_celltype)
  names(data) <- paste0('k', 1:6)
  return(data)
}

# deng topic model fits
load_deng_topics = function(){
  nmf <- readRDS('data/deng/nmf.rds')
  snmf <- readRDS('data/deng/snmf.rds')

  hs <- org.Hs.eg.db::org.Hs.eg.db
  gene_symbols <- unique(rownames(nmf$fl$L.pm))
  gene_symbols <- map_chr(gene_symbols, toupper)
  map2entrez <- generate_map2entrez(gene_symbols, 'SYMBOL')

  extract_component = function(fl, k){
    dt <- do.call('cbind', map(c('L.pm', 'L.psd', 'L.lfsr'), ~ pluck(fl, .x)[, k])) %>%
    as.data.frame()
    names(dt) <- c('postMean', 'postMeanSE', 'lfsr')
    dt <-
      dt %>%
      rownames_to_column(var='SYMBOL') %>%
      mutate(SYMBOL = toupper(SYMBOL)) %>%
      left_join(map2entrez) %>%
      relocate(ENTREZID, .after=SYMBOL) %>%
      filter(!is.na(ENTREZID)) %>%
      mutate(
        threshold.on = lfsr,
        beta = postMean,
        se = postMeanSE
      )
    return(dt)
  }

  nmf.standardized = map(2:20, ~extract_component(nmf$fl, .x))
  names(nmf.standardized) <- paste0('NMF_Topic', 2:20)

  snmf.standardized = map(2:20, ~extract_component(snmf$fl, .x))
  names(snmf.standardized) <- paste0('SNMF_Topic', 2:20)

  data <- c(nmf.standardized, snmf.standardized)
  return(data)
}

# human chimp eb for many celltypes
load_human_chimp_eb  = function(){
  de <- readRDS('data/human_chimp_eb/big.df.rds')
  gene_symbols <- unique(de$gene)
  map2enrtrez <- generate_map2entrez(gene_symbols, 'SYMBOL')
  data <- de %>%
    dplyr::rename(SYMBOL=gene) %>%
    left_join(map2enrtrez, by='SYMBOL') %>%
    relocate(ENTREZID, .after=SYMBOL) %>%
    mutate(  # set default columns
      beta = dream.logFC,
      se = dream.SE,
      threshold.on = dream.p.val
    ) %>%
    group_by(celltype) %>%
    group_map(~ .x, .keep = T) %>%
    add_names(map_chr(., ~pluck(.x, 'celltype')[1]))
  return(data)
}


load_baboon_diet = function(){
  de.adipose <- read.table('data/wenhe_baboon_diet/DE_lrt_adipose.txt')
  de.liver <- read.table('data/wenhe_baboon_diet/DE_lrt_liver.txt')
  de.muscle <- read.table('data/wenhe_baboon_diet/DE_lrt_muscle.txt')

  de <- bind_rows(
    rownames_to_column(de.adipose) %>% mutate(tissue='Adipose'),
    rownames_to_column(de.liver) %>% mutate(tissue='Liver'),
    rownames_to_column(de.muscle) %>% mutate(tissue='Muscle')) %>%
    rename(gene=rowname)

  hs <- org.Hs.eg.db::org.Hs.eg.db
  gene_symbols <- unique(de$gene)
  symbol2entrez <- AnnotationDbi::select(
    hs, keys=gene_symbols,
    columns=c('ENTREZID', 'SYMBOL'),
    keytype = 'SYMBOL')

  add_names = function(l, n){
    names(l) <- n
    return(l)
  }

  data <- de %>%
    rename(SYMBOL=gene) %>%
    left_join(symbol2entrez, by='SYMBOL') %>%
    relocate(ENTREZID, .after=SYMBOL) %>%
    mutate(  # set default columns
      beta = logFC,
      se = 0,
      threshold.on = PValue
    ) %>%
    group_by(tissue) %>%
    group_map(~ .x, .keep = T) %>%
    add_names(map_chr(., ~pluck(.x, 'tissue')[1]))

  return(data)
}


load_chondrocyte_data<- function(){
  data <- read.table('data/anthony/DE_table_midSize.csv', sep=',', header=T)
  map2entrez <- generate_map2entrez(data$Row.names, 'ENSEMBL')
  data <- data %>%
    mutate(ENSEMBL = Row.names, beta=1., se=1.) %>%
    select(ENSEMBL, adjPval_CTS, adjPval_il1b, beta, se) %>%
    left_join(map2entrez)

  cts <- data %>%
    select(ENTREZID, beta, se, adjPval_CTS) %>%
    mutate(threshold.on = adjPval_CTS)

  il1b <- data %>%
    select(ENTREZID, beta, se, adjPval_il1b) %>%
    mutate(threshold.on = adjPval_il1b)

  data <- list(
    CTS=cts, IL1B=il1b
  )
  return(data)
}


# NOTE: data has logfc etc, but we're just loading p-values here
load_chondrocyte_data2 <- function(){
  data <- read.table('data/anthony/DE_merged.csv', sep=',', header=T)
  map2entrez <- generate_map2entrez(data$rowname, 'ENSEMBL')
  data <- data %>%
    mutate(ENSEMBL = rowname, beta=1., se=1.) %>%
    select(ENSEMBL, adj.P.Val_CTS, adj.P.Val_il1b, beta, se) %>%
    left_join(map2entrez)

  cts <- data %>%
    select(ENTREZID, beta, se, adj.P.Val_CTS) %>%
    mutate(threshold.on = adj.P.Val_CTS)

  il1b <- data %>%
    select(ENTREZID, beta, se, adj.P.Val_il1b) %>%
    mutate(threshold.on = adj.P.Val_il1b)

  data <- list(
    CTS=cts, IL1B=il1b
  )
  return(data)
}


# yusha's singcle cell tumor data
load_sc_tumor_hnscc = function(){
  # load data
  genes <- read.table('data/yusha_sc_tumor/hnscc/gene_list.txt')$V1

  files <- Sys.glob('data/yusha_sc_tumor/hnscc/factor*', dirmark = FALSE)
  map2entrez <- generate_map2entrez(genes, 'SYMBOL')

  f <- files[1]

  load_factor = function(f){
    gene_list <- read.table(f)$V1

    map2entrez %>%
      mutate(
        beta = as.integer(SYMBOL %in% gene_list),
        threshold.on= 1 - beta
      )
  }

  data <- purrr::map(files, load_factor)

  # get names from filepaths
  names <- stringr::str_split(files, '/') %>%
    map_chr(~pluck(.x, 4)) %>%
    stringr::str_split('.txt') %>%
    map_chr(~pluck(.x, 1))
  names(data) <- names

  return(data)
}


load_sc_tumor_pdac = function(){
  # load data
  genes <- read.table('data/yusha_sc_tumor/pdac/gene_list.txt')$V1

  files <- Sys.glob('data/yusha_sc_tumor/pdac/factor*', dirmark = FALSE)
  map2entrez <- generate_map2entrez(genes, 'SYMBOL')

  f <- files[1]

  load_factor = function(f){
    gene_list <- read.table(f)$V1

    map2entrez %>%
      mutate(
        beta = as.integer(SYMBOL %in% gene_list),
        threshold.on= 1 - beta
      )
  }

  data <- purrr::map(files, load_factor)

  # get names from filepaths
  names <- stringr::str_split(files, '/') %>%
    map_chr(~pluck(.x, 4)) %>%
    stringr::str_split('.txt') %>%
    map_chr(~pluck(.x, 1))
  names(data) <- names

  return(data)
}
