#' Returns a list of variable genes based on PCA
#'
#'  @description  Extract genes, i.e. "features", based on the top
#' loadings of principal components
#' formed from the bulk expression data set
#'
#' @param sc_avg_expr Expression data frame. Rownames are genes,
#' colnames are single cell cluster name, and
#' values are average single cell expression (log transformed).
#' @param bulk_expr Bulk RNA expression data frame. Rownames are genes,
#' colnames are bulk sample names, and values are expression counts.
#' @param nr_pcs Number of PCs to selected gene loadings from.
#' See the explore_PCA_corr.Rmd vignette for details.
#' @param percentile Select the percentile of absolute values of
#' PCA loadings to select genes from. E.g. 0.999 would select the
#' top point 1 percent of genes with the largest loadings.
#' @return The list of genes to use as features.
#'
#' @export
getPCAGenes <- function(
                        sc_avg_expr,
                        bulk_expr,
                        nr_pcs,
                        percentile) {

  # Get overlapping genes
  shared.genes <- rownames(sc_avg_expr)[rownames(sc_avg_expr) %in% rownames(bulk_expr)]
  bulkrna <- bulk_expr[shared.genes, ]
  bulkrna <- log(bulkrna + 1)

  # Get the PCs
  pca <- prcomp(t(bulkrna))

  # For the given number PCs, select the genes with the largest loadings
  genes <- c()
  for (i in 1:nr_pcs) {
    cutoff <- quantile(abs(pca$rotation[, i]), probs = percentile)
    genes <- c(genes, rownames(pca$rotation[abs(pca$rotation[, i]) >= cutoff, ]))
  }

  return(genes)
}
