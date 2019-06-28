#' Run GSEA to compare a gene list(s) to per cell or per cluster expression data
#' @description Use fgsea algorithm to compute normalized enrichment scores and pvalues for gene
#' set ovelap
#' @param expr_mat single-cell expression matrix or Seurat object
#' @param query_genes A vector or named list of vectors of genesets of interest to compare via GSEA. If
#' supplying a named list, then the gene set names will appear in the output.
#' @param cluster_ids vector of cell cluster assignments, supplied as a vector with order that
#' matches columns in `expr_mat`. Not required if running per cell.
#' @param n_perm Number of permutation for fgsea function. Defaults to 1000.
#' @param per_cell if true run per cell, otherwise per cluster.
#' @param scale convert expr_mat into zscores prior to running GSEA?, default = FALSE
#' @param no_warnings suppress warnings from gsea ties
#' @return dataframe of gsea scores (pval, NES), with clusters as rownames
#' @export
run_gsea <- function(expr_mat,
                     query_genes,
                     cluster_ids = NULL,
                     n_perm = 1000,
                     per_cell = FALSE,
                     scale = FALSE,
                     no_warnings = TRUE) {
  if (!is.list(query_genes)) {
    geneset_list <- list("query_genes" = query_genes)
  } else {
    geneset_list <- query_genes
  }

  if (!per_cell & (ncol(expr_mat) != length(cluster_ids))) {
    stop("cluster_ids do not match number of cells (columns) in expr_mat ")
  }

  if (n_perm > 1e4 & per_cell) {
    warning("run_gsea() take a long time if running many permutations and running per cell")
  }

  if (scale) {
    expr_mat <- t(scale(t(as.matrix(expr_mat))))
  }

  if (!per_cell) {
    avg_mat <- average_clusters(expr_mat, cluster_info = cluster_ids)
  } else {
    avg_mat <- expr_mat
  }

  res <- list()
  for (i in seq_along(colnames(avg_mat))) {
    if (!(no_warnings)) {
      gsea_res <- fgsea::fgsea(geneset_list,
        avg_mat[, i],
        minSize = 1,
        maxSize = max(sapply(geneset_list, length)),
        nproc = 1,
        nperm = n_perm
      )
    } else {
      suppressWarnings(gsea_res <- fgsea::fgsea(geneset_list,
        avg_mat[, i],
        minSize = 1,
        maxSize = max(sapply(geneset_list, length)),
        nproc = 1,
        nperm = n_perm
      ))
    }
    res[[i]] <- gsea_res[, c("pathway", "pval", "NES")]
  }
  gsea_res <- dplyr::bind_rows(res)
  gsea_res <- as.data.frame(dplyr::mutate(gsea_res, cell = colnames(avg_mat)))
  gsea_res <- tibble::column_to_rownames(gsea_res, "cell")

  gsea_res
}
