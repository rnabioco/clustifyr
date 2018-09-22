#' Main function to compare scRNA-seq data to bulk RNA-seq data.
#'
#' @param expr_mat single-cell expression matrix or Seurat object
#' @param metadata cell cluster assignments, supplied as a vector or data.frame. If
#' data.frame is supplied then `cluster_col` and/or `cell_col` need to be set. Not required if running correlation per cell.
#' @param bulk_mat bulk expression matrix
#' @param query_genes A vector of genes of interest to compare. If NULL, then common genes between
#' the expr_mat and bulk_mat will be used for comparision.
#' @param cluster_col column in metadata that contains cluster ids per cell. Will default to first
#' column of metadata if not supplied. Not required if running correlation per cell.
#' @param per_cell if true run per cell, otherwise per cluster.
#' @param num_perm number of permutations, set to 0 by default
#' @param compute_method method(s) for computing similarity scores
#' @param use_var_genes if providing a seurat object, use the variable genes
#'  (stored in seurat_object@var.genes) as the query_genes.
#' @param ... additional arguments to pass to compute_method function
#' @export
clustify <- function(expr_mat,
                     bulk_mat,
                     metadata,
                     query_genes = NULL,
                     cluster_col = NULL,
                     per_cell = FALSE,
                     num_perm = 0,
                     compute_method = "spearman",
                     use_var_genes = FALSE,
                     ...) {
  UseMethod("clustify", expr_mat)
}

#' @rdname clustify
#' @export
clustify.default <- function(expr_mat,
                             bulk_mat,
                             metadata = NULL,
                             query_genes = NULL,
                             cluster_col = NULL,
                             per_cell = FALSE,
                             num_perm = 0,
                             compute_method = "spearman",
                             ...) {

  if(!compute_method %in% clustifyr_methods){
    stop(paste(compute_method, "correlation method not implemented"))
  }

  # select gene subsets
  gene_constraints <- get_common_elements(rownames(expr_mat),
                                           rownames(bulk_mat),
                                           query_genes)

  expr_mat <- expr_mat[gene_constraints, , drop = FALSE]
  bulk_mat <- bulk_mat[gene_constraints, , drop = FALSE]

  if(is.null(metadata) & !per_cell) {
    stop("metadata needed for per cluster analysis")
  }

  if(!per_cell){
    if(is.vector(metadata)){
      cluster_ids <- metadata
    } else if (is.data.frame(metadata) & !is.null(cluster_col)){
      cluster_ids <- metadata[[cluster_col]]
    } else {
      stop("metadata not formatted correctly,
           supply either a character vector or a dataframe")
    }
  }

  if(per_cell){
    cluster_ids <- colnames(expr_mat)
  }

  if (num_perm == 0) {
    res <- get_similarity(
      expr_mat,
      bulk_mat,
      cluster_ids = cluster_ids,
      per_cell = per_cell,
      compute_method = compute_method, ...
    )
  } else {
    # run permutation
    res <- permute_similarity(
      expr_mat,
      bulk_mat,
      cluster_ids = cluster_ids,
      num_perm = num_perm,
      per_cell = per_cell,
      compute_method = compute_method, ...
    )
  }

  return(res)
}

#' @rdname clustify
#' @export
clustify.seurat <- function(seurat_object,
                            bulk_mat,
                            query_genes = NULL,
                            per_cell = FALSE,
                            num_perm = 0,
                            cluster_col = NULL,
                            compute_method = "pearson",
                            use_var_genes = FALSE) {
  expr_mat <- seurat_object@data
  metadata <- seurat_object@meta.data

  if (use_var_genes){
    query_genes <- seurat_object@var.genes
  }

  res <- clustify(expr_mat,
                  bulk_mat,
                  metadata,
                  query_genes,
                  per_cell = per_cell,
                  num_perm = num_perm,
                  cluster_col = cluster_col,
                  compute_method = compute_method
  )

  res
}

#' Correlation functions available in clustifyR
#'@export
clustifyr_methods <- c(
  "pearson",
  "spearman",
  "cosine",
  "kl_divergence"
)
