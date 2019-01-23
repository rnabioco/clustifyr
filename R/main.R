#' Main function to compare scRNA-seq data to bulk RNA-seq data.
#'
#'@export
clustify <- function(input, ...) {
  UseMethod("clustify", input)
}

#' @rdname clustify
#' @param input single-cell expression matrix or Seurat object
#' @param metadata cell cluster assignments, supplied as a vector or data.frame. If
#' data.frame is supplied then `cluster_col` needs to be set. Not required if running correlation per cell.
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
#' @param dr stored dimension reduction
#' @param seurat_out output cor matrix or called seurat object
#' @param ... additional arguments to pass to compute_method function
#' @export
clustify.default <- function(input,
                             bulk_mat,
                             metadata = NULL,
                             query_genes = NULL,
                             cluster_col = NULL,
                             per_cell = FALSE,
                             num_perm = 0,
                             compute_method = "spearman",
                             ...) {
  expr_mat <- input
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
#' @param dr stored dimension reduction
#' @param seurat_out output cor matrix or called seurat object
#' @param ... additional arguments to pass to compute_method function

#' @export
clustify.seurat <- function(input,
                            bulk_mat,
                            query_genes = NULL,
                            per_cell = FALSE,
                            num_perm = 0,
                            cluster_col = NULL,
                            compute_method = "spearman",
                            use_var_genes = TRUE,
                            dr = "tsne",
                            seurat_out = TRUE,
                            threshold = 0,
                            ...) {
  s_object <- input
  expr_mat <- s_object@data
  metadata <- use_seurat_meta(s_object, dr = dr)

  if (use_var_genes & is.null(query_genes)){
    query_genes <- s_object@var.genes
  }

  res <- clustify(expr_mat,
                  bulk_mat,
                  metadata,
                  query_genes,
                  per_cell = per_cell,
                  num_perm = num_perm,
                  cluster_col = cluster_col,
                  compute_method = compute_method,
                  ...
  )

  if (seurat_out == F) {
    res
  } else {
    col_meta <- colnames(metadata)
    if("type" %in% col_meta | "type2" %in% col_meta){
      warning('metadata column name clash of "type"/"type2"')
      return()
    }
    if (num_perm != 0) {
      res <- -log(res$p_val+.01,10)
    }
    df_temp <- tibble::as_tibble(res, rownames = "rn")
    df_temp <- tidyr::gather(df_temp, key = type, value = r, -rn)
    df_temp[["type"]][df_temp$r < threshold] <- paste0("r<", threshold,", unassigned")
    df_temp <- dplyr::top_n(dplyr::group_by_at(df_temp, 1), 1, r)
    if (nrow(df_temp) != nrow(res)) {
      clash <- dplyr::group_by(df_temp, rn)
      clash <- dplyr::summarize(clash, n = n())
      clash <- dplyr::filter(clash, n > 1)
      clash <- dplyr::pull(clash, rn)
      df_temp <- dplyr::mutate(df_temp, type = ifelse(rn %in% clash, paste0(type, "-CLASH!"), type))
      df_temp <- dplyr::distinct(df_temp, rn, r, .keep_all = T)
    }
    if (per_cell == F) {
      df_temp <- dplyr::rename(df_temp, !!cluster_col:=rn)
      df_temp_full <- dplyr::left_join(tibble::rownames_to_column(metadata, "rn"), df_temp, by = cluster_col)
      df_temp_full <- tibble::column_to_rownames(df_temp_full, "rn")
      } else {
      df_temp_full <- dplyr::left_join(tibble::rownames_to_column(metadata, "rn"), df_temp, by = "rn")
      df_temp_full <- tibble::column_to_rownames(df_temp_full, "rn")
    }
    tryCatch(s_object@meta.data <- df_temp_full, error=function(e) print("seurat not loaded, returning cor_mat instead"), finally = return(res))
    s_object
  }
}

#' Correlation functions available in clustifyR
#'@export
clustifyr_methods <- c(
  "pearson",
  "spearman",
  "cosine",
  "kl_divergence"
)
