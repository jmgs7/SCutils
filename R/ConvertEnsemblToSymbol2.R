ConvertEnsembleToSymbol2 <- function(
  mat,
  mirror = NULL,
  species = c('human', 'mouse')
) {
  species <- match.arg(arg = species)
  if (species == 'human') {
    database <- 'hsapiens_gene_ensembl'
    symbol <- 'hgnc_symbol'
  } else if (species == 'mouse') {
    database <- 'mmusculus_gene_ensembl'
    symbol <- 'mgi_symbol'
  } else {
    stop('species name not found')
  }

  library("biomaRt")
  library("dplyr")

  name_df <- data.frame(gene_id = c(rownames(mat)))
  name_df$orig.id <- name_df$gene_id
  #make this a character, otherwise it will throw errors with left_join
  name_df$gene_id <- as.character(name_df$gene_id)
  # in case it's gencode, this mostly works
  #if ensembl, will leave it alone
  name_df$gene_id <- sub("[.][0-9]*", "", name_df$gene_id)
  ensembl <- useEnsembl(
    biomart = "ensembl",
    dataset = database,
    mirror = mirror
  )
  mart <- useDataset(dataset = database, ensembl, )
  genes <- name_df$gene_id
  gene_IDs <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", symbol),
    values = genes,
    mart = mart
  )
  gene.df <- left_join(name_df, gene_IDs, by = c("gene_id" = "ensembl_gene_id"))
  rownames(gene.df) <- make.unique(gene.df$orig.id)
  gene.df <- gene.df[rownames(mat), ]
  gene.df <- gene.df[gene.df[, symbol] != '', ]
  gene.df <- gene.df[!is.na(gene.df$orig.id), ]
  mat.filter <- mat[gene.df$orig.id, ]
  rownames(mat.filter) <- make.unique(gene.df[, symbol])
  return(mat.filter)
}
