#' Compute the p-value for similarity using permutation
#'
#' @description Permute cluster labels to calculate empirical p-value
#'
#' @param expr_mat single-cell expression matrix
#' @param ref_mat reference expression matrix
#' @param cluster_ids vector of cluster ids for each cell
#' @param compute_method method(s) for computing similarity scores
#' @param per_cell run per cell?
#' @param ... additional parameters not used yet
#' @export
get_similarity <- function(expr_mat,
                           ref_mat,
                           cluster_ids,
                           compute_method,
                           per_cell = FALSE,
                           ...) {
  ref_clust <- colnames(ref_mat)
  if (sum(is.na(cluster_ids)) > 0) {
    print("reassigning NAs to unknown")
    cluster_ids <- factor(cluster_ids)
    cluster_ids <- factor(cluster_ids, levels = c(levels(cluster_ids), NA), labels = c(levels(cluster_ids), "unknown"), exclude = NULL)
  }
  if (!per_cell) {
    sc_clust <- sort(unique(cluster_ids))
    # print(sc_clust)
    clust_avg <- compute_mean_expr(
      expr_mat,
      cluster_ids,
      sc_clust
    )
  } else {
    sc_clust <- cluster_ids
    clust_avg <- expr_mat
  }

  assigned_score <- calc_similarity(
    clust_avg,
    ref_mat,
    compute_method,
    ...
  )

  rownames(assigned_score) <- sc_clust
  colnames(assigned_score) <- ref_clust

  return(assigned_score)
}

#' Compute a p-value for similarity using permutation
#'
#' @description Permute cluster labels to calculate empirical p-value
#'
#'
#' @param expr_mat single-cell expression matrix
#' @param cluster_ids clustering info of single-cell data assume that
#'  genes have ALREADY BEEN filtered
#' @param ref_mat reference expression matrix
#' @param num_perm number of permutations
#' @param per_cell run per cell?
#' @param compute_method method(s) for computing similarity scores
#' @param ... additional parameters
#' @export
permute_similarity <- function(expr_mat,
                               ref_mat,
                               cluster_ids,
                               num_perm,
                               per_cell = F,
                               compute_method, ...) {
  ref_clust <- colnames(ref_mat)

  if (!per_cell) {
    sc_clust <- sort(unique(cluster_ids))
    clust_avg <- compute_mean_expr(
      expr_mat,
      cluster_ids,
      sc_clust
    )
  } else {
    sc_clust <- colnames(expr_mat)
    clust_avg <- expr_mat
  }

  assigned_score <- calc_similarity(
    clust_avg,
    ref_mat,
    compute_method,
    ...
  )

  # perform permutation
  sig_counts <- matrix(0L, nrow = length(sc_clust), ncol = length(ref_clust))

  for (i in 1:num_perm) {
    resampled <- sample(cluster_ids,
      length(cluster_ids),
      replace = FALSE
    )

    if (!per_cell) {
      permuted_avg <- compute_mean_expr(
        expr_mat,
        resampled,
        sc_clust
      )
    } else {
      permuted_avg <- expr_mat[, resampled, drop = FALSE]
    }

    # permutate assignment
    new_score <- calc_similarity(
      permuted_avg,
      ref_mat,
      compute_method,
      ...
    )
    sig_counts <- sig_counts + as.numeric(new_score > assigned_score)
  }

  rownames(assigned_score) <- sc_clust
  colnames(assigned_score) <- ref_clust
  rownames(sig_counts) <- sc_clust
  colnames(sig_counts) <- ref_clust

  return(list(
    score = assigned_score,
    p_val = sig_counts / num_perm
  ))
}

#' compute mean of clusters
#' @export
compute_mean_expr <- function(expr_mat, sc_assign, sc_clust) {
  sapply(sc_clust, function(x) Matrix::rowMeans(expr_mat[, sc_assign == x, drop = FALSE]))
}

#' compute similarity
#' @export
calc_similarity <- function(sc_avg,
                            ref_mat,
                            compute_method, ...) {

  # use stats::cor matrix method if possible
  if (any(compute_method %in% c("pearson", "spearman"))) {
    similarity_score <- cor(as.matrix(sc_avg),
      ref_mat,
      method = compute_method
    )
    return(similarity_score)
  }

  # remove 0s ?
  # if(any(compute_method %in% c("spearman"))) {
  #   print("considering 0 as missing data")
  #   sc_avg[sc_avg == 0] <- NA
  #   similarity_score <- cor(as.matrix(sc_avg),
  #                           ref_mat, method = "spearman", use = "pairwise.complete.obs")
  #   return(similarity_score)
  # }

  sc_clust <- colnames(sc_avg)
  ref_clust <- colnames(ref_mat)
  features <- intersect(rownames(sc_avg), rownames(ref_mat))
  sc_avg <- sc_avg[features, ]
  ref_mat <- ref_mat[features, ]
  similarity_score <- matrix(NA,
    nrow = length(sc_clust),
    ncol = length(ref_clust)
  )
  for (i in seq_along(sc_clust)) {
    for (j in seq_along(ref_clust)) {
      similarity_score[i, j] <- vector_similarity(
        sc_avg[, sc_clust[i]],
        ref_mat[, ref_clust[j]],
        compute_method, ...
      )
    }
  }
  return(similarity_score)
}

#' Compute similarity between two vectors
#'
#' @description Compute the similarity score between two vectors using a
#' customized scoring function
#' Two vectors may be from either scRNA-seq or bulk RNA-seq data.
#' The lengths of vec1 and vec2 must match, and must be arranged in the
#' same order of genes.
#' Both vectors should be provided to this function after pre-processing,
#' feature selection and dimension reduction.
#'
#' @param vec1 test vector
#' @param vec2 reference vector
#' @param compute_method method to run i.e. corr_coef
#' @param ... arguments to pass to compute_method function
#' @export
vector_similarity <- function(vec1, vec2, compute_method, ...) {
  # examine whether two vectors are of the same size
  if (!is.numeric(vec1) || !is.numeric(vec2) || length(vec1) != length(vec2)) {
    stop("compute_similarity: two input vectors are not numeric or of different sizes.")
  }

  if (!(compute_method %in% c("cosine", "kl_divergence"))) {
    stop(paste(compute_method, "not implemented"))
  }

  if (compute_method == "kl_divergence") {
    res <- kl_divergence(vec1, vec2, ...)
  } else if (compute_method == "cosine") {
    res <- cosine(vec1, vec2, ...)
  }
  # return the similarity score, must be
  return(res)
}

#' Cosine distance
#' @param vec1 test vector
#' @param vec2 reference vector
#' @export
cosine <- function(vec1, vec2) {
  sum(vec1 * vec2) / sqrt(sum(vec1^2) * sum(vec2^2))
}
#' KL divergence
#'
#' @description Use package entropy to compute Kullback-Leibler divergence.
#' The function first converts each vector's reads to pseudo-number of
#' transcripts by normalizing the total reads to total_reads.
#' The normalized read for each gene is then rounded to serve as the
#' pseudo-number of transcripts.
#' Function [entropy::KL.shrink()] is called to compute the KL-divergence between
#' the two vectors, and the maximal allowed divergence is set to max_KL.
#' Finally, a linear transform is performed to convert the KL divergence,
#' which is between 0 and max_KL, to a similarity score between -1 and 1.
#'
#' @param vec1 Test vector
#' @param vec2 Reference vector
#' @param if_logcounts Whether the vectors are log-transformed. If so, the
#' raw count should be computed before computing KL-divergence.
#' @param total_reads Pseudo-library size
#' @param max_KL Maximal allowed value of KL-divergence.
#' @export
kl_divergence <- function(vec1, vec2, if_logcounts = FALSE,
                          total_reads = 1000, max_KL = 1) {
  if (if_logcounts) {
    vec1 <- expm1(vec1)
    vec2 <- expm1(vec2)
  }
  count1 <- round(vec1 * total_reads / sum(vec1))
  count2 <- round(vec2 * total_reads / sum(vec2))
  est_KL <- entropy::KL.shrink(count1, count2,
    unit = "log",
    verbose = FALSE
  )
  return((max_KL - est_KL) / max_KL * 2 - 1)
}
