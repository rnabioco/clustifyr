#' Main function to compare scRNA-seq data to bulk RNA-seq data.
#'
#' @param expr_mat single-cell expression matrix
#' @param metadata clustering info of single-cell data
#' @param bulk_mat bulk expression matrix
#' @param query_genes A vector of genes of interest.
#' @param cluster_col column in metadata that contains cluster ids per cell. Will default to first
#' column of metadata if not supplied. Not required if running correlation per cell.
#' @param cell_col column in metadata that contains the cell ids that match to columns
#' in the single-cell expression matrix. Will default to rownames of the metadata data.frame
#' if not supplied. Not required if running correlation per cluster.
#' @param per_cell if true run per cell, otherwise per cluster.
#' @param num_perm number of permutations, set to 0 by default
#' @param return_full Return full results includeing scores instead
#'  of correlation coefficient only.
#' @param compute_method method(s) for computing similarity scores
#' @param ... additional arguments to pass to compute_method function
#' @export
clustify <- function(expr_mat, metadata, bulk_mat, query_genes, cluster_col = NULL,
                     cell_col = NULL,
                     per_cell = F,
                     num_perm = 0, return_full = F, compute_method = "spearman", ...) {
  # select gene subsets
  gene_constraints <- list(rownames(expr_mat), rownames(bulk_mat), query_genes)
  expr_mat <- select_gene_subset(expr_mat, gene_constraints)
  bulk_mat <- select_gene_subset(bulk_mat, gene_constraints)

  if(!compute_method %in% clustifyr_methods){
    stop(paste(compute_method, "correlation method not implemented"))
  }

  if(per_cell){
    if(is.null(cell_col)){
      cluster_ids <- rownames(metadata)
    } else {
      cluster_ids <- metadata[[cell_col]]
    }
  } else { # per cluster
    if(is.null(cluster_col)){
      cluster_ids <- metadata[[colnames(metadata)[1]]]
    } else {
      cluster_ids <- metadata[[cluster_col]]
    }
  }

  # run permutation
  res <- permutation_similarity(
    expr_mat, bulk_mat, cluster_ids = cluster_ids,
    num_perm = num_perm, per_cell = per_cell,
    compute_method = compute_method, ...
  )

  # extract score only by default
  if (!return_full) {
    res <- res$score
  }

  return(res)
}

#' Correlation functions available in clustifyR
#'@export
clustifyr_methods <- c(
  "pearson",
  "spearman",
  "cosine",
  "kl_divergence"
)
