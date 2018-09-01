#' Reduce expression matrix to variable genes and binarize.
#'
#' @description Reduce the full expression expr_matrix to only the
#' highly variable genes,
#' change the readcounts to binary factors,
#' prepare the input data for random forest
#'
#' @param expr_mat expression expr_matrix with row names as the gene
#' names and each column for different cell
#' @param metadata cluster and cell type info for each cell
#' @param vargene a list of gene names for highly variable genes
#' @export
reduce_expr_matrix <- function(expr_mat, metadata, vargene) {
  Expexpr_mat <- as.matrix(expr_mat)
  # dim(Expexpr_mat)
  # dim(Vargenes)
  # colnames(Vargenes) <- c("Genes")

  # Reduce the expression expr_mat to only highly variable genes, and transpose
  Expexpr_mat <- t(Expexpr_mat[rownames(Expexpr_mat) %in% vargene, ])
  # dim(Expexpr_mat)

  # Convert >0 values to 1
  Expexpr_mat[Expexpr_mat > 0] <- 1
  Expexpr_mat <- data.frame(Cells = row.names(Expexpr_mat), Expexpr_mat)

  # metadata <- data.frame(Cells = row.names(metadata), metadata)
  metadata <- as.data.frame(metadata)
  names(metadata)[names(metadata) == "rn"] <- "Cells"
  Merge <- merge(metadata, Expexpr_mat, by = "Cells")
  data1 <- as.data.frame(lapply(Merge[ncol(metadata):ncol(Merge)], as.factor))
  return(data1)
}
