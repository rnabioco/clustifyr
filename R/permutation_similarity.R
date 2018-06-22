data("pbmc4k_matrix"); data("pbmc4k_meta"); data("pbmc4k_avg"); data("pbmc_bulk_matrix");
gene_constraints <- list(rownames(pbmc4k_avg), rownames(pbmc_bulk_matrix));
sc_expr <- select_gene_subset(pbmc4k_matrix, gene_constraints);
bulk_expr <- select_gene_subset(pbmc_bulk_matrix, gene_constraints);
res1 <- permutation_similarity(sc_expr, bulk_expr, pbmc4k_meta, 1000, corr_coef, method="pearson");
res2 <- permutation_similarity(sc_expr, bulk_expr, pbmc4k_meta, 1000, corr_coef, method="spearman");
res3 <- permutation_similarity(sc_expr, bulk_expr, pbmc4k_meta, 1000, corr_coef, method="cosine");
res4 <- permutation_similarity(sc_expr, bulk_expr, pbmc4k_meta, 1000, corr_coef, method="kl_divergence");

# compute the p-value for data set
# sc_expr, bulk_expr: single-cell expr matrix and bulk epxr matrix.
# sc_meta: clustering info of single-cell data
# assume that genes have ALREADY BEEN filtered
# num_perm: number of permutation
# compute_method, ...: parameters feed in for computing similarity score
permutation_similarity <- function(sc_expr, bulk_expr, sc_meta, num_perm, compute_method, ...) {
  # get cell types
  sc_clust <- sort(unique(sc_meta[,'cluster'])); bulk_clust <- colnames(bulk_expr);
  assigned_score <- run_one_round(compute_mean_expr(sc_expr, sc_meta[,'cluster'], sc_clust), bulk_expr, sc_clust, bulk_clust, compute_method, ...);

  # perform permutation
  sig_counts <- matrix(0L, nrow=length(sc_clust), ncol=length(bulk_clust));
  for (i in 1:num_perm) {
    # permutate assignment
    new_score <- run_one_round(compute_mean_expr(sc_expr, sample(sc_meta[,'cluster'], nrow(sc_meta), replace=FALSE), sc_clust), bulk_expr, sc_clust, bulk_clust, compute_method, ...);
    sig_counts <- sig_counts + as.numeric(new_score>assigned_score);
  }
  rownames(assigned_score) <- sc_clust; colnames(assigned_score) <- bulk_clust;
  rownames(sig_counts) <- sc_clust; colnames(sig_counts) <- bulk_clust;
  return(list(score=assigned_score, p_val=sig_counts/num_perm));
}

compute_mean_expr <- function(sc_expr, sc_assign, sc_clust) {
  sapply(sc_clust, function(x) Matrix::rowMeans(sc_expr[, sc_assign==x]))
}

run_one_round <- function(sc_avg, bulk_expr, sc_clust, bulk_clust, compute_method, ...) {
  num_sc_clust <- length(sc_clust); num_bulk_clust <- length(bulk_clust);
  similarity_score <- matrix(NA, nrow=num_sc_clust, ncol=num_bulk_clust);
  for (i in 1:num_sc_clust) {
    for (j in 1:num_bulk_clust) {
      similarity_score[i,j] <- compute_similarity(sc_avg[, sc_clust[i]], bulk_expr[, bulk_clust[j]], compute_method, ...);
    }
  }
  return(similarity_score);
}

