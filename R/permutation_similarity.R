
#' Compute the p-value for similarity using permutation
#'
#' @description Permute cluster labels to calculate empirical p-value
#'
#'
#' @param expr_mat single-cell expression matrix
#' @param metadata clustering info of single-cell data assume that
#'  genes have ALREADY BEEN filtered
#' @param bulk_mat bulk expression matrix
#' @param num_perm number of permutations
#' @param per_cell run per cell?
#' @param compute_method method(s) for computing similarity scores
#' @param cluster_col cluster_col column used for clustering
#' @param ... additional paramters to pass to run_one_round
#' @export
permutation_similarity <- function(expr_mat, metadata, bulk_mat, num_perm, per_cell = F,
                                   compute_method, cluster_col = "cluster", ...) {

  # get cluster or cell ids
  clust_ids <- dplyr::pull(metadata, cluster_col)

  # get cell types
  sc_clust <- sort(unique(clust_ids))
  bulk_clust <- colnames(bulk_mat)

  if (!per_cell){
    clust_avg <- compute_mean_expr(expr_mat,
                                   clust_ids,
                                   sc_clust)
  } else {
    clust_avg <- expr_mat
  }

  assigned_score <- run_one_round(clust_avg,
                                  bulk_mat,
                                  sc_clust,
                                  bulk_clust,
                                  compute_method,
                                  ...)

  if(num_perm > 0) {
    # perform permutation
    sig_counts <- matrix(0L, nrow=length(sc_clust), ncol=length(bulk_clust))

    for (i in 1:num_perm) {
      # permutate assignment
      new_score <- run_one_round(compute_mean_expr(expr_mat,
                                                   sample(dplyr::pull(metadata,cluster_col),
                                                          nrow(metadata),
                                                          replace=FALSE),
                                                   sc_clust),
                                 bulk_mat,
                                 sc_clust,
                                 bulk_clust,
                                 compute_method,
                                 ...)
      sig_counts <- sig_counts + as.numeric(new_score>assigned_score)
    }
  } else {
    sig_counts <- matrix(NA, nrow=length(sc_clust), ncol=length(bulk_clust))
  }

  rownames(assigned_score) <- sc_clust
  colnames(assigned_score) <- bulk_clust
  rownames(sig_counts) <- sc_clust
  colnames(sig_counts) <- bulk_clust
  return(list(score=assigned_score, p_val=sig_counts/num_perm))
}

#' compute mean of clusters
#' @noRd
compute_mean_expr <- function(expr_mat, sc_assign, sc_clust) {
  sapply(sc_clust, function(x) Matrix::rowMeans(expr_mat[, sc_assign==x]))
}

#' compute similarity
#' @noRd
run_one_round <- function(sc_avg, bulk_mat, sc_clust, bulk_clust, compute_method, ...) {
  num_sc_clust <- length(sc_clust)
  num_bulk_clust <- length(bulk_clust)
  similarity_score <- matrix(NA, nrow=num_sc_clust, ncol=num_bulk_clust)
  for (i in 1:num_sc_clust) {
    for (j in 1:num_bulk_clust) {
      similarity_score[i,j] <- compute_similarity(sc_avg[, sc_clust[i]], bulk_mat[, bulk_clust[j]], compute_method, ...)
    }
  }
  return(similarity_score)
}

