#' reduce the full expression matrix to only the highly variable genes, change the readcounts to binary factors
#' prepare the input data for random forest
#' @param  mat: expression matrix with row names as the gene names and each column for different cell
#' @param meta: cluster and cell type info for each cell
#' @param vargene: a list of gene names for highly variable genes
#' @export

Reduce_matrix <- function(mat, meta, vargene){
  
  Expmat <- as.matrix(mat)
  # dim(Expmat)
  # dim(Vargenes)
  # colnames(Vargenes) <- c("Genes")
  
  #Reduce the expression mat to only highly variable genes, and transpose
  Expmat <- t(Expmat[rownames(Expmat) %in% vargene, ])
  # dim(Expmat)
  
  #Convert >0 values to 1
  Expmat[Expmat > 0] <- 1
  Expmat <- data.frame(Cells = row.names(Expmat), Expmat)
  
  # Meta <- data.frame(Cells = row.names(meta), meta)
  Meta <- as.data.frame(meta)
  names(Meta)[names(Meta) == "rn"] <- "Cells"
  Merge <- merge(Meta, Expmat, by="Cells")
  data1 <- as.data.frame(lapply(Merge[ncol(Meta):ncol(Merge)], as.factor))
  return(data1)
}
