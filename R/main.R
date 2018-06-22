#' main: global controller function to evaluate corr. coef
#' @param query_gene_list character vector with genes for computing similarity
#' @inheritParams permutation_similarity
#' @export
run_cor <- function(sc_expr, sc_meta, bulk_expr, query_gene_list, per_cell = F,
                    if_permute=TRUE, num_permute=1000, return_full = F,
                    compute_method, ...) {
  # not permute -> num_permute = 0
  if (!if_permute) {
    num_permute <- 0
  }

  # select gene subsets
  gene_constraints <- list(rownames(sc_expr), rownames(bulk_expr), query_gene_list)
  sc_expr <- select_gene_subset(sc_expr, gene_constraints)
  bulk_expr <- select_gene_subset(bulk_expr, gene_constraints)

  # run permutation
  res <- permutation_similarity(sc_expr, bulk_expr, sc_meta,
                                num_permute, per_cell, compute_method, ...)

  # extract score only by default
  if(!return_full) {
    res <- res$score
  }

  return(res)
}
