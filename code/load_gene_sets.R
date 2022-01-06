make.gene.set <- function(X){
  referenceGene <- rownames(X)
  geneSets <- colnames(X)
  a <- which(X!=0, arr.ind = T)
  geneSet <- cbind(referenceGene[a[,1]], geneSets[a[,2]])
  colnames(geneSet) <- c('gene', 'geneSet')
  geneSet <- as_tibble(geneSet)
  return(geneSet)
}

make.gene.set.matrix <- function(geneSet){
  X <- geneSet %>%
    add_column(in_geneSet = 1) %>%
    select(geneSet, in_geneSet, gene) %>%
    pivot_wider(names_from = geneSet, values_from = in_geneSet) %>%
    mutate(across(everything(), ~replace_na(.x, 0)))

  backgroundGene <- X$gene
  X_mat <- as.matrix(X[,-1])
  rownames(X_mat) <- X$gene
  return(Matrix(X_mat, sparse = T))
}

load.msigdb.X <- function(){
  msigdb.tb <- msigdbr::msigdbr(species="Homo sapiens", category = c("C2"))

  msigdb.geneSet <- msigdb.tb %>%
    select(gs_id, human_entrez_gene) %>%
    transmute(geneSet = gs_id, gene = human_entrez_gene)

  X <- make.gene.set.matrix(msigdb.geneSet)
  X <- X[, sample(colnames(X), 500)]
  X <- X[rowSums(X) > 0, ]
}

load.webGestalt.X <- function(db='geneontology_Biological_Process_noRedundant'){
  organism <- 'hsapiens'
  interestGeneType = "ensembl_gene_id"
  referenceGeneType = "ensembl_gene_id"
  outputDirectory = './data/WebGestalt/results'
  hostName = 'http://www.webgestalt.org/'

  enrichDatabase <- c(db)
  geneSet <- WebGestaltR::loadGeneSet(
    organism=organism, enrichDatabase=enrichDatabase, hostName=hostName)$geneSet

  X <- make.gene.set.matrix(geneSet)
  return(X)
}


load.gobp <- function(){
  X <- load.webGestalt.X(db='geneontology_Biological_Process')
  sizes <- colSums(X)
  X <- X[, (sizes > 5) & (sizes < 500)]
}
