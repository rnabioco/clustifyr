
#' compute the p-value for similarity using permutation
#'
#' @description Permute cluster labels to calculate empirical p-value
#' @param sc_expr: single-cell expr matrix
#' @param bulk_expr: bulk epxr matrix.
#' @param sc_meta: clustering info of single-cell data assume that genes have ALREADY BEEN filtered
#' @param num_perm: number of permutation
#' @param compute_method parameters feed in for computing similarity score
#' @param metadata metadata column used for clustering
#' @param ...
#' @export
permutation_similarity <- function(sc_expr, bulk_expr, sc_meta, num_perm, compute_method, metadata = "cluster", ...) {
  # get cell types
  sc_clust <- sort(unique(pull(sc_meta, metadata)))
  bulk_clust <- colnames(bulk_expr)
  assigned_score <- run_one_round(compute_mean_expr(sc_expr,
                                                    pull(sc_meta,metadata),
                                                    sc_clust),
                                  bulk_expr,
                                  sc_clust,
                                  bulk_clust,
                                  compute_method,
                                  ...)

  # perform permutation
  sig_counts <- matrix(0L, nrow=length(sc_clust), ncol=length(bulk_clust))

  for (i in 1:num_perm) {
    # permutate assignment
    new_score <- run_one_round(compute_mean_expr(sc_expr,
                                                 sample(pull(sc_meta,metadata),
                                                        nrow(sc_meta),
                                                        replace=FALSE),
                                                 sc_clust),
                               bulk_expr,
                               sc_clust,
                               bulk_clust,
                               compute_method,
                               ...)
    sig_counts <- sig_counts + as.numeric(new_score>assigned_score)
  }

  rownames(assigned_score) <- sc_clust
  colnames(assigned_score) <- bulk_clust
  rownames(sig_counts) <- sc_clust
  colnames(sig_counts) <- bulk_clust
  return(list(score=assigned_score, p_val=sig_counts/num_perm))
}

#' compute mean of clusters
#' @noRd
compute_mean_expr <- function(sc_expr, sc_assign, sc_clust) {
  sapply(sc_clust, function(x) Matrix::rowMeans(sc_expr[, sc_assign==x]))
}

#' compute similarity
#' @noRd
run_one_round <- function(sc_avg, bulk_expr, sc_clust, bulk_clust, compute_method, ...) {
  num_sc_clust <- length(sc_clust)
  num_bulk_clust <- length(bulk_clust)
  similarity_score <- matrix(NA, nrow=num_sc_clust, ncol=num_bulk_clust)
  for (i in 1:num_sc_clust) {
    for (j in 1:num_bulk_clust) {
      similarity_score[i,j] <- compute_similarity(sc_avg[, sc_clust[i]], bulk_expr[, bulk_clust[j]], compute_method, ...)
    }
  }
  return(similarity_score)
}

