#' Main function to compare scRNA-seq data to bulk RNA-seq data.
#'
#' @export
clustify <- function(input, ...) {
  UseMethod("clustify", input)
}

#' @rdname clustify
#' @param input single-cell expression matrix or Seurat object
#' @param metadata cell cluster assignments, supplied as a vector or data.frame. If
#' data.frame is supplied then `cluster_col` needs to be set. Not required if running correlation per cell.
#' @param ref_mat reference expression matrix
#' @param query_genes A vector of genes of interest to compare. If NULL, then common genes between
#' the expr_mat and ref_mat will be used for comparision.
#' @param cluster_col column in metadata that contains cluster ids per cell. Will default to first
#' column of metadata if not supplied. Not required if running correlation per cell.
#' @param per_cell if true run per cell, otherwise per cluster.
#' @param num_perm number of permutations, set to 0 by default
#' @param compute_method method(s) for computing similarity scores
#' @param use_var_genes if providing a seurat object, use the variable genes
#'  (stored in seurat_object@var.genes) as the query_genes.
#' @param dr stored dimension reduction
#' @param seurat_out output cor matrix or called seurat object
#' @param verbose whether to report certain variables chosen
#' @param ... additional arguments to pass to compute_method function
#' @export
clustify.default <- function(input,
                             ref_mat,
                             metadata = NULL,
                             query_genes = NULL,
                             cluster_col = NULL,
                             per_cell = FALSE,
                             num_perm = 0,
                             compute_method = "spearman",
                             verbose = F,
                             ...) {
  expr_mat <- input
  if (!compute_method %in% clustifyr_methods) {
    stop(paste(compute_method, "correlation method not implemented"))
  }

  # select gene subsets
  gene_constraints <- get_common_elements(
    rownames(expr_mat),
    rownames(ref_mat),
    query_genes
  )

  if (verbose == T) {
    print(paste0("using # of genes: ", length(gene_constraints)))
  }

  expr_mat <- expr_mat[gene_constraints, , drop = FALSE]
  ref_mat <- ref_mat[gene_constraints, , drop = FALSE]

  if (is.null(metadata) & !per_cell) {
    stop("metadata needed for per cluster analysis")
  }

  if (!per_cell) {
    if (is.vector(metadata)) {
      cluster_ids <- metadata
    } else if (is.data.frame(metadata) & !is.null(cluster_col)) {
      cluster_ids <- metadata[[cluster_col]]
    } else {
      stop("metadata not formatted correctly,
           supply either a character vector or a dataframe")
    }
  }

  if (per_cell) {
    cluster_ids <- colnames(expr_mat)
  }

  if (num_perm == 0) {
    res <- get_similarity(
      expr_mat,
      ref_mat,
      cluster_ids = cluster_ids,
      per_cell = per_cell,
      compute_method = compute_method, ...
    )
  } else {
    # run permutation
    res <- permute_similarity(
      expr_mat,
      ref_mat,
      cluster_ids = cluster_ids,
      num_perm = num_perm,
      per_cell = per_cell,
      compute_method = compute_method, ...
    )
  }

  return(res)
}

#' @rdname clustify
#' @param ref_mat reference expression matrix
#' @param query_genes A vector of genes of interest to compare. If NULL, then common genes between
#' the expr_mat and ref_mat will be used for comparision.
#' @param cluster_col column in metadata that contains cluster ids per cell. Will default to first
#' column of metadata if not supplied. Not required if running correlation per cell.
#' @param per_cell if true run per cell, otherwise per cluster.
#' @param num_perm number of permutations, set to 0 by default
#' @param compute_method method(s) for computing similarity scores
#' @param use_var_genes if providing a seurat object, use the variable genes
#'  (stored in seurat_object@var.genes) as the query_genes.
#' @param dr stored dimension reduction
#' @param seurat_out output cor matrix or called seurat object
#' @param verbose whether to report certain variables chosen
#' @param ... additional arguments to pass to compute_method function

#' @export
clustify.seurat <- function(input,
                            ref_mat,
                            query_genes = NULL,
                            per_cell = FALSE,
                            num_perm = 0,
                            cluster_col = NULL,
                            compute_method = "spearman",
                            use_var_genes = TRUE,
                            dr = "tsne",
                            seurat_out = TRUE,
                            threshold = 0,
                            verbose = F,
                            ...) {
  s_object <- input
  # for seurat < 3.0
  expr_mat <- s_object@data
  metadata <- use_seurat_meta(s_object, dr = dr, seurat3 = F)

  if (use_var_genes & is.null(query_genes)) {
    query_genes <- s_object@var.genes
  }

  res <- clustify(expr_mat,
    ref_mat,
    metadata,
    query_genes,
    per_cell = per_cell,
    num_perm = num_perm,
    cluster_col = cluster_col,
    compute_method = compute_method,
    verbose = verbose,
    ...
  )

  if (seurat_out == F) {
    res
  } else {
    col_meta <- colnames(metadata)
    if ("type" %in% col_meta | "type2" %in% col_meta) {
      warning('metadata column name clash of "type"/"type2"')
      return()
    }
    if (num_perm != 0) {
      res <- -log(res$p_val + .01, 10)
    }
    df_temp <- tibble::as_tibble(res, rownames = "rn")
    df_temp <- tidyr::gather(df_temp, key = type, value = r, -rn)
    df_temp[["type"]][df_temp$r < threshold] <- paste0("r<", threshold, ", unassigned")
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
      df_temp <- dplyr::rename(df_temp, !!cluster_col := rn)
      df_temp_full <- dplyr::left_join(tibble::rownames_to_column(metadata, "rn"), df_temp, by = cluster_col)
      df_temp_full <- tibble::column_to_rownames(df_temp_full, "rn")
    } else {
      df_temp_full <- dplyr::left_join(tibble::rownames_to_column(metadata, "rn"), df_temp, by = "rn")
      df_temp_full <- tibble::column_to_rownames(df_temp_full, "rn")
    }
    if ("Seurat" %in% loadedNamespaces()) {
      s_object@meta.data <- df_temp_full
      return(s_object)
    } else {
      print("seurat not loaded, returning cor_mat instead")
      return(res)
    }
    s_object
  }
}

#' @rdname clustify
#' @param ref_mat reference expression matrix
#' @param query_genes A vector of genes of interest to compare. If NULL, then common genes between
#' the expr_mat and ref_mat will be used for comparision.
#' @param cluster_col column in metadata that contains cluster ids per cell. Will default to first
#' column of metadata if not supplied. Not required if running correlation per cell.
#' @param per_cell if true run per cell, otherwise per cluster.
#' @param num_perm number of permutations, set to 0 by default
#' @param compute_method method(s) for computing similarity scores
#' @param use_var_genes if providing a seurat object, use the variable genes
#'  (stored in seurat_object@var.genes) as the query_genes.
#' @param dr stored dimension reduction
#' @param seurat_out output cor matrix or called seurat object
#' @param verbose whether to report certain variables chosen
#' @param ... additional arguments to pass to compute_method function

#' @export
clustify.Seurat <- function(input,
                            ref_mat,
                            query_genes = NULL,
                            per_cell = FALSE,
                            num_perm = 0,
                            cluster_col = NULL,
                            compute_method = "spearman",
                            use_var_genes = TRUE,
                            dr = "tsne",
                            seurat_out = TRUE,
                            threshold = 0,
                            verbose = F,
                            ...) {
  s_object <- input
  # for seurat 3.0 +
  expr_mat <- s_object@assays$RNA@data
  metadata <- use_seurat_meta(s_object, dr = dr, seurat3 = T)

  if (use_var_genes & is.null(query_genes)) {
    query_genes <- s_object@assays$RNA@var.features
  }

  res <- clustify(expr_mat,
    ref_mat,
    metadata,
    query_genes,
    per_cell = per_cell,
    num_perm = num_perm,
    cluster_col = cluster_col,
    compute_method = compute_method,
    verbose = verbose,
    ...
  )

  if (seurat_out == F) {
    res
  } else {
    col_meta <- colnames(metadata)
    if ("type" %in% col_meta | "type2" %in% col_meta) {
      warning('metadata column name clash of "type"/"type2"')
      return()
    }
    if (num_perm != 0) {
      res <- -log(res$p_val + .01, 10)
    }
    df_temp <- tibble::as_tibble(res, rownames = "rn")
    df_temp <- tidyr::gather(df_temp, key = type, value = r, -rn)
    df_temp[["type"]][df_temp$r < threshold] <- paste0("r<", threshold, ", unassigned")
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
      df_temp <- dplyr::rename(df_temp, !!cluster_col := rn)
      df_temp_full <- dplyr::left_join(tibble::rownames_to_column(metadata, "rn"), df_temp, by = cluster_col)
      df_temp_full <- tibble::column_to_rownames(df_temp_full, "rn")
    } else {
      df_temp_full <- dplyr::left_join(tibble::rownames_to_column(metadata, "rn"), df_temp, by = "rn")
      df_temp_full <- tibble::column_to_rownames(df_temp_full, "rn")
    }
    if ("Seurat" %in% loadedNamespaces()) {
      s_object@meta.data <- df_temp_full
      return(s_object)
    } else {
      print("seurat not loaded, returning cor_mat instead")
      return(res)
    }
    s_object
  }
}
#' Correlation functions available in clustifyR
#' @export
clustifyr_methods <- c(
  "pearson",
  "spearman",
  "cosine",
  "kl_divergence"
)

#' Main function to compare scRNA-seq data to gene lists.
#'
#' @export
clustify_lists <- function(input, ...) {
  UseMethod("clustify_lists", input)
}

#' @rdname clustify_lists
#' @param input single-cell expression matrix or Seurat object
#' @param per_cell compare per cell or per cluster
#' @param cluster_info data.frame or vector containing cluster assignments per cell.
#' Order must match column order in supplied matrix. If a data.frame
#' provide the cluster_col parameters.
#' @param log_scale input data is natural log,
#' averaging will be done on unlogged data
#' @param cluster_col column in cluster_info with cluster number
#' @param topn number of top expressing genes to keep from input matrix
#' @param cut expression cut off from input matrix
#' @param marker matrix or dataframe of candidate genes for each cluster
#' @param marker_inmatrix whether markers genes are already in preprocessed matrix form
#' @param genomen number of genes in the genome
#' @param metric adjusted p-value for hypergeometric test, or jaccard index
#' @param output_high if true (by default to fit with rest of package),
#' -log10 transform p-value
#' @param ... passed to matrixize_markers
#'
#' @export

clustify_lists.default <- function(input,
                                   per_cell = F,
                                   cluster_info = NULL,
                                   log_scale = T,
                                   cluster_col = "cluster",
                                   topn = 3000,
                                   cut = 0,
                                   marker,
                                   marker_inmatrix = T,
                                   genomen = 30000,
                                   metric = "hyper",
                                   output_high = TRUE,
                                   ...) {
  if (per_cell == F) {
    input <- average_clusters(input,
      cluster_info,
      log_scale = log_scale,
      cluster_col = cluster_col
    )
  }

  bin_input <- binarize_expr(input, n = topn, cut = cut)

  if (marker_inmatrix != T) {
    marker <- matrixize_markers(
      marker,
      ...
    )
  }

  compare_lists(bin_input,
    marker_m = marker,
    n = genomen,
    metric = metric,
    output_high = output_high
  )
}

#' @rdname clustify_lists
#' @param input seurat object
#' @param per_cell compare per cell or per cluster
#' @param cluster_info data.frame or vector containing cluster assignments per cell.
#' Order must match column order in supplied matrix. If a data.frame
#' provide the cluster_col parameters.
#' @param log_scale input data is natural log,
#' averaging will be done on unlogged data
#' @param cluster_col column in cluster_info with cluster number
#' @param topn number of top expressing genes to keep from input matrix
#' @param cut expression cut off from input matrix
#' @param marker matrix or dataframe of candidate genes for each cluster
#' @param marker_inmatrix whether markers genes are already in preprocessed matrix form
#' @param genomen number of genes in the genome
#' @param metric adjusted p-value for hypergeometric test, or jaccard index
#' @param output_high if true (by default to fit with rest of package),
#' -log10 transform p-value
#' @param dr stored dimension reduction
#' @param seurat_out output cor matrix or called seurat object

#' @param ... passed to matrixize_markers
#'
#' @export
clustify_lists.seurat <- function(input,
                                  per_cell = F,
                                  cluster_info = NULL,
                                  log_scale = T,
                                  cluster_col = "cluster",
                                  topn = 3000,
                                  cut = 0,
                                  marker,
                                  marker_inmatrix = T,
                                  genomen = 30000,
                                  metric = "hyper",
                                  output_high = TRUE,
                                  dr = "tsne",
                                  seurat_out = TRUE,
                                  threshold = 0,
                                  ...) {
  s_object <- input
  # for seurat < 3.0
  input <- s_object@data
  cluster_info <- as.data.frame(use_seurat_meta(s_object, dr = dr, seurat3 = F))

  res <- clustify_lists(input,
    per_cell = per_cell,
    cluster_info = cluster_info,
    log_scale = log_scale,
    cluster_col = cluster_col,
    topn = topn,
    cut = cut,
    marker,
    marker_inmatrix = marker_inmatrix,
    genomen = genomen,
    metric = metric,
    output_high = output_high,
    ...
  )

  if (seurat_out == F) {
    res
  } else {
    if (per_cell == F) {
      df_temp <- cor_to_call(res, cluster_info, col = cluster_col, threshold = threshold)
      df_temp_full <- dplyr::left_join(tibble::rownames_to_column(cluster_info, "rn"), df_temp, by = cluster_col)
      df_temp_full <- tibble::column_to_rownames(df_temp_full, "rn")
    } else {
      df_temp <- cor_to_call(res, cluster_info, col = "rn", threshold = threshold)
      df_temp_full <- dplyr::left_join(tibble::rownames_to_column(cluster_info, "rn"), df_temp, by = "rn")
      df_temp_full <- tibble::column_to_rownames(df_temp_full, "rn")
    }
    if ("Seurat" %in% loadedNamespaces()) {
      s_object@meta.data <- df_temp_full
      return(s_object)
    } else {
      print("seurat not loaded, returning cor_mat instead")
      return(res)
    }
  }
}

#' @rdname clustify_lists
#' @param input seurat object
#' @param per_cell compare per cell or per cluster
#' @param cluster_info data.frame or vector containing cluster assignments per cell.
#' Order must match column order in supplied matrix. If a data.frame
#' provide the cluster_col parameters.
#' @param log_scale input data is natural log,
#' averaging will be done on unlogged data
#' @param cluster_col column in cluster_info with cluster number
#' @param topn number of top expressing genes to keep from input matrix
#' @param cut expression cut off from input matrix
#' @param marker matrix or dataframe of candidate genes for each cluster
#' @param marker_inmatrix whether markers genes are already in preprocessed matrix form
#' @param genomen number of genes in the genome
#' @param metric adjusted p-value for hypergeometric test, or jaccard index
#' @param output_high if true (by default to fit with rest of package),
#' -log10 transform p-value
#' @param dr stored dimension reduction
#' @param seurat_out output cor matrix or called seurat object

#' @param ... passed to matrixize_markers
clustify_lists.Seurat <- function(input,
                                  per_cell = F,
                                  cluster_info = NULL,
                                  log_scale = T,
                                  cluster_col = "cluster",
                                  topn = 3000,
                                  cut = 0,
                                  marker,
                                  marker_inmatrix = T,
                                  genomen = 30000,
                                  metric = "hyper",
                                  output_high = TRUE,
                                  dr = "tsne",
                                  seurat_out = TRUE,
                                  threshold = 0,
                                  ...) {
  s_object <- input
  # for seurat 3.0 +
  input <- s_object@assays$RNA@data
  cluster_info <- as.data.frame(use_seurat_meta(s_object, dr = dr, seurat3 = T))

  res <- clustify_lists(input,
    per_cell = per_cell,
    cluster_info = cluster_info,
    log_scale = log_scale,
    cluster_col = cluster_col,
    topn = topn,
    cut = cut,
    marker,
    marker_inmatrix = marker_inmatrix,
    genomen = genomen,
    metric = metric,
    output_high = output_high,
    ...
  )

  if (seurat_out == F) {
    res
  } else {
    if (per_cell == F) {
      df_temp <- cor_to_call(res, cluster_info, col = cluster_col, threshold = threshold)
      df_temp_full <- dplyr::left_join(tibble::rownames_to_column(cluster_info, "rn"), df_temp, by = cluster_col)
      df_temp_full <- tibble::column_to_rownames(df_temp_full, "rn")
    } else {
      df_temp <- cor_to_call(res, cluster_info, col = "rn", threshold = threshold)
      df_temp_full <- dplyr::left_join(tibble::rownames_to_column(cluster_info, "rn"), df_temp, by = "rn")
      df_temp_full <- tibble::column_to_rownames(df_temp_full, "rn")
    }
    if ("Seurat" %in% loadedNamespaces()) {
      s_object@meta.data <- df_temp_full
      return(s_object)
    } else {
      print("seurat not loaded, returning cor_mat instead")
      return(res)
    }
  }
}
