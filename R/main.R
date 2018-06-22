data("pbmc4k_matrix"); data("pbmc4k_avg"); data("pbmc_bulk_matrix");
# main: global controller function to evaluate and visualize corr. coef
main <- function(sc_expr, sc_meta, bulk_expr, query_gene_list, if_permute=TRUE, num_permute=1000, compute_method, ...) {
  # not permute -> num_permute = 0
  if (!if_permute) {
    num_permute <- 0;
  }

  # select gene subsets
  gene_constraints <- list(rownames(sc_avg), rownames(bulk_expr), query_gene_list);
  sc_expr <- sc_expr[gene_constraints, ]; bulk_expr <- bulk_expr[gene_constraints, ];

  # run permutation
  res <- permutation_similarity(sc_expr, bulk_expr, sc_meta, num_permute, compute_method, ...);

  # extract information
  similarity_score <- res$score;
  p_val <- res$p_val;

  #TODO: PLOT FUNCTION

}
