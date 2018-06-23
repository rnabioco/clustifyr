#' run_cor: Main function to compare scRNA-seq data to bulk RNA-seq data.
#' @inheritParams permutation_similarity
#' @param query_gene_list A vector of genes of interest.
#' @param return_full Return full results instead of correlation coefficient only.
#' @export
run_cor <- function(expr_mat, metadata, bulk_mat, query_gene_list, per_cell = F,
                    num_perm = 0, return_full = F, compute_method, ...) {
# select gene subsets
  gene_constraints <- list(rownames(expr_mat), rownames(bulk_mat), query_gene_list)
  expr_mat <- select_gene_subset(expr_mat, gene_constraints)
  bulk_mat <- select_gene_subset(bulk_mat, gene_constraints)

  # run permutation
  res <- permutation_similarity(expr_mat, bulk_mat, metadata,
                                num_permute, per_cell, compute_method, ...)

  # extract score only by default
  if(!return_full) {
    res <- res$score
  }

  return(res)
}
