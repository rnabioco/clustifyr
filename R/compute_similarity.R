
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

  if (!(compute_method %in% c("cosine", "kl_divergence"))){
    stop(paste(compute_method, "not implemented"))
  }

  if(compute_method == "kl_divergence"){
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
cosine <- function(vec1, vec2){
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
kl_divergence <- function(vec1, vec2, if_logcounts = FALSE, total_reads = 1000, max_KL = 1) {
  if (if_logcounts) {
    vec1 <- expm1(vec1)
    vec2 <- expm1(vec2)
  }
  count1 <- round(vec1 * total_reads / sum(vec1))
  count2 <- round(vec2 * total_reads / sum(vec2))
  est_KL <- entropy::KL.shrink(count1, count2, unit = "log",
                               verbose = FALSE)
  return((max_KL - est_KL) / max_KL * 2 - 1)
}
