#' Average expression values per cluster
#'
#' @param mat expression matrix
#' @param cluster_info data.frame or vector containing cluster assignments per cell.
#' Order must match column order in supplied matrix. If a data.frame
#' provide the cluster_col parameters.
#' @param if_log input data is natural log,
#' averaging will be done on unlogged data
#' @param cluster_col column in cluster_info with cluster number
#' @param cell_col if provided, will reorder matrix first
#' @param low_threshold option to remove clusters with too few cells
#' @param method whether to take mean (default) or median
#' @param output_log whether to report log results
#'
#' @export
average_clusters <- function(mat, cluster_info,
                             if_log = TRUE,
                             cluster_col = "cluster",
                             cell_col = NULL,
                             low_threshold = 0,
                             method = "mean",
                             output_log = TRUE) {
  if (!(is.null(cell_col))) {
    if (!(all(colnames(mat) == cluster_info[[cell_col]]))) {
      mat <- mat[, cluster_info[[cell_col]]]
    }
  }
  if (is.vector(cluster_info)) {
    cluster_ids <- split(colnames(mat), cluster_info)
  } else if (is.data.frame(cluster_info) & !is.null(cluster_col)) {
    cluster_ids <- split(colnames(mat), cluster_info[[cluster_col]])
  } else if (class(cluster_info) == "factor") {
    cluster_info <- as.character(cluster_info)
    cluster_ids <- split(colnames(mat), cluster_info)
  } else {
    stop("cluster_info not formatted correctly,
         supply either a  vector or a dataframe")
  }

  if (method == "mean") {
    out <- lapply(
      cluster_ids,
      function(cell_ids) {
        if (!all(cell_ids %in% colnames(mat))) {
          stop("cell ids not found in input matrix")
        }
        if (if_log) {
          mat_data <- expm1(mat[, cell_ids, drop = FALSE])
        } else {
          mat_data <- mat[, cell_ids, drop = FALSE]
        }
        res <- Matrix::rowMeans(mat_data, na.rm = TRUE)
        if (output_log) {
          res <- log1p(res)
        }
        res
      }
    )
  } else {
    out <- lapply(
      cluster_ids,
      function(cell_ids) {
        if (!all(cell_ids %in% colnames(mat))) {
          stop("cell ids not found in input matrix")
        }
        mat_data <- mat[, cell_ids, drop = FALSE]
        res <- apply(mat_data, 1, function(x) {
          median(x[x > 0])
        })
        res[is.na(res)] <- 0
        res
      }
    )
  }

  out <- do.call(cbind, out)
  if (low_threshold > 0) {
    fil <- sapply(cluster_ids, FUN = length) > low_threshold
    out <- out[, as.vector(fil)]
  }
  return(out)
}

#' Percentage detected per cluster
#'
#' @param mat expression matrix
#' @param cluster_info data.frame with cells
#' @param cluster_col column in cluster_info with cluster number
#' @param cut_num binary cutoff for detection
#'
#' @export
percent_clusters <- function(mat, cluster_info,
                             cluster_col = "cluster",
                             cut_num = 0.5) {
  mat[mat >= cut_num] <- 1
  mat[mat <= cut_num] <- 0

  average_clusters(mat, cluster_info,
    if_log = FALSE,
    cluster_col = cluster_col
  )
}

#' Function to make best call from correlation matrix
#'
#' @param cor_mat correlation matrix
#'
#' @export
get_best_match_matrix <- function(cor_mat) {
  best_mat <- as.data.frame(t(apply(cor_mat, 1, function(x) x - max(x))))
  best_mat[best_mat == 0] <- "1"
  best_mat[best_mat != "1"] <- "0"

  return(best_mat)
}

#' Function to make call and attach score
#'
#' @param name name of row to query
#' @param best_mat binarized call matrix
#' @param cor_mat correlation matrix
#' @param carry_cor whether the correlation score gets reported
#'
#' @export
get_best_str <- function(name,
                         best_mat,
                         cor_mat,
                         carry_cor = TRUE) {
  if (sum(as.numeric(best_mat[name, ])) > 0) {
    best.names <- colnames(best_mat)[which(best_mat[name, ] == 1)]
    best.cor <- round(cor_mat[name, which(best_mat[name, ] == 1)], 2)
    for (i in 1:length(best.cor)) {
      if (i == 1) {
        str <- paste0(best.names[i], " (", best.cor[i], ")")
      } else {
        str <- paste0(str, "; ", best.names[i], " (", best.cor[i], ")")
      }
    }
  } else {
    str <- "?"
  }

  if (carry_cor == FALSE) {
    str <- gsub(" \\(.*\\)", "", str)
  }
  return(str)
}

#' Find entries shared in all vectors
#' @description return entries found in all supplied vectors. If the vector supplied
#' is NULL or NA, then it will be excluded from the comparision.
#' @param ... vectors
#' @examples
#' a <- rep(1:5)
#' b <- rep(4:10)
#' c <- rep(4:6)
#' get_common_elements(a, b, c)
#' @export
get_common_elements <- function(...) {
  vecs <- list(...)
  # drop NULL elements of list
  vecs <- vecs[!sapply(vecs, is.null)]
  # drop NA elements of list (NA values OK in a vector)
  vecs <- vecs[!is.na(vecs)]

  Reduce(intersect, vecs)
}

#' Intra-experiment cluster projection for one sample/set to the rest
#'
#' @param expr_mat single-cell expression matrix or Seurat object
#' @param metadata cell cluster assignments, supplied as a vector or data.frame. If
#' data.frame is supplied then `cluster_col` needs to be set. Not required if running correlation per cell.
#' @param query_genes A vector of genes of interest to compare. If NULL, then common genes between
#' the expr_mat and ref_mat will be used for comparision.
#' @param cluster_col column in metadata that contains cluster ids per cell. Will default to first
#' column of metadata if not supplied. Not required if running correlation per cell.
#' @param sample_col column in metadata that contains sample/subset info
#' @param sample_id ids in column to serve as reference
#' @param per_cell if true run per cell, otherwise per cluster.
#' @param compute_method method(s) for computing similarity scores
#' @param ... additional arguments to pass to compute_method function
#'
#' @export
clustify_intra <- function(expr_mat,
                           metadata,
                           query_genes,
                           cluster_col,
                           sample_col,
                           sample_id,
                           per_cell = FALSE,
                           compute_method = "spearman",
                           ...) {
  row_ref <- (metadata[[sample_col]] == sample_id)
  expr_mat_ref <- expr_mat[, row_ref]
  expr_mat_tar <- expr_mat[, !row_ref]
  meta_ref <- metadata[row_ref, ]
  meta_tar <- metadata[!row_ref, ]

  avg_clusters_ref <- average_clusters(expr_mat_ref, meta_ref,
    if_log = FALSE,
    cluster_col = cluster_col
  )

  r2 <- clustify(expr_mat_tar, avg_clusters_ref, meta_tar,
    query_genes = query_genes,
    cluster_col = cluster_col,
    per_cell = per_cell,
    num_perm = 0,
    compute_method = compute_method,
    use_var_genes = FALSE
  )

  r2
}

#' Average expression values per cluster, filtered by set parameter, defaults to calculating background
#'
#' @param mat expression matrix
#' @param cluster_info data.frame or vector containing cluster assignments per cell, and attribute to filter on.
#' Order must match column order in supplied matrix. If a data.frame
#' provide the cluster_col parameters.
#' @param if_log input data is natural log,
#' averaging will be done on unlogged data
#' @param filter_on column in cluster_info to filter on
#' @param group_by column name to use for cluster identity
#' @param filter_method "<", "==", ">" compared to filter_value
#' @param filter_value baseline minimum as background cutoff
#'
#' @export
average_clusters_filter <- function(mat, cluster_info,
                                    if_log = TRUE,
                                    filter_on = "nGene",
                                    group_by = NULL,
                                    filter_method = "<=",
                                    filter_value = 300) {
  cell_ids <- 0
  eval(parse(text = paste0("cell_ids <- cluster_info[[filter_on]] ", filter_method, "filter_value")))
  if (sum(cell_ids) == 0) {
    stop("no cells kept after filtering")
  }

  if (!is.null(group_by)) {
    res <- average_clusters(
      mat, 
      cluster_info,
      if_log = if_log,
      cluster_col = group_by
    )
  } else {
    if (if_log) {
      mat_data <- expm1(mat[, cell_ids, drop = FALSE])
    } else {
      mat_data <- mat[, cell_ids, drop = FALSE]
    }
    res <- Matrix::rowMeans(mat_data)
    if (if_log) {
      res <- log1p(res)
    }
    res
  }
}

#' Remove high background expression genes from matrix
#'
#' @param mat expression matrix
#' @param background vector or dataframe or matrix of high expression genes in background
#' @param n the number of top genes to exclude, 0 defaults to all
#'
#' @export

remove_background <- function(mat, background, n = 0) {
  if (n == 0) {
    n <- length(background)
  }

  if (!is.vector(background)) {
    background <- background[order(background[, 1], decreasing = TRUE), , drop = FALSE]
    background <- rownames(background)[1:n]
  } else if (!is.null(names(background))) {
    background <- names(sort(background, decreasing = TRUE)[1:n])
  }

  mat[!(rownames(mat) %in% background), ]
}

#' Convert expression matrix to GSEA pathway scores (would take a similar place in workflow before average_clusters/binarize)
#'
#' @param mat expression matrix
#' @param pathway_list a list of vectors, each named for a specific pathway, or dataframe
#' @param n_perm Number of permutation for fgsea function. Defaults to 1000.
#' @param scale convert expr_mat into zscores prior to running GSEA?, default = FALSE
#' @param no_warnings suppress warnings from gsea ties

#' @export

calculate_pathway_gsea <- function(mat,
                                   pathway_list,
                                   n_perm = 1000,
                                   scale = TRUE,
                                   no_warnings = TRUE) {
  # pathway_list can be user defined or `my_pathways <- fgsea::reactomePathways(rownames(pbmc4k_matrix))`
  out <- lapply(
    names(pathway_list),
    function(y) {
      marker_list <- list()
      marker_list[[1]] <- pathway_list[[y]]
      names(marker_list) <- y
      v1 <- marker_list
      temp <- run_gsea(mat, v1,
        n_perm = n_perm,
        scale = scale,
        per_cell = TRUE,
        no_warnings = no_warnings
      )
      temp <- temp[, 3, drop = FALSE]
    }
  )
  res <- do.call(cbind, out)
  colnames(res) <- names(pathway_list)
  res
}

#' get best calls for each cluster
#'
#' @param correlation_matrix input similarity matrix
#' @param metadata input metadata with tsne coordinates and cluster ids
#' @param col metadata column, can be cluster or cellid
#' @param collapse_to_cluster if a column name is provided, takes the most frequent call of entire cluster to color in plot
#' @param threshold minimum correlation coefficent cutoff for calling clusters
#' @param rename_suff suffix to add to type and r column names
#'
#' @export
cor_to_call <- function(correlation_matrix,
                        metadata = NULL,
                        cluster_col = "cluster",
                        collapse_to_cluster = FALSE,
                        threshold = 0,
                        rename_suff = NULL) {
  df_temp <- tibble::as_tibble(correlation_matrix, rownames = cluster_col)
  df_temp <- tidyr::gather(df_temp, key = type, value = r, -!!cluster_col)
  df_temp[["type"]][df_temp$r < threshold] <- paste0("r<", threshold, ", unassigned")
  df_temp <- dplyr::top_n(dplyr::group_by_at(df_temp, 1), 1, r)
  if (nrow(df_temp) != nrow(correlation_matrix)) {
    clash <- dplyr::summarize(dplyr::group_by_at(df_temp, 1), n = n())
    clash <- dplyr::filter(clash, n > 1)
    clash <- dplyr::pull(clash, 1)
    df_temp[lapply(df_temp[, 1], FUN = function(x) x %in% clash)[[1]], 2] <- paste0(df_temp[["type"]][lapply(df_temp[, 1], FUN = function(x) x %in% clash)[[1]]], "-CLASH!")
    df_temp <- dplyr::distinct(df_temp, exclude = "type", .keep_all = TRUE)
    df_temp_full <- dplyr::select(df_temp, -exclude)
  } else {
    df_temp_full <- df_temp
  }

  if (collapse_to_cluster != FALSE) {
    if (!(cluster_col %in% colnames(metadata))) {
      metadata <- tibble::as_tibble(metadata, rownames = cluster_col)
    }
    df_temp_full <- dplyr::left_join(df_temp_full, metadata, by = cluster_col)
    df_temp_full[, "type2"] <- df_temp_full[[collapse_to_cluster]]
    df_temp_full2 <- dplyr::group_by(df_temp_full, type, type2)
    df_temp_full2 <- dplyr::summarize(df_temp_full2, sum = sum(r), n = n())
    df_temp_full2 <- dplyr::group_by(df_temp_full2, type2)
    df_temp_full2 <- dplyr::arrange(df_temp_full2, desc(n), desc(sum))
    df_temp_full2 <- dplyr::filter(df_temp_full2, type != paste0("r<", threshold, ", unassigned"))
    df_temp_full2 <- dplyr::slice(df_temp_full2, 1)
    df_temp_full2 <- dplyr::right_join(df_temp_full2, dplyr::select(df_temp_full, -type), by = stats::setNames(collapse_to_cluster, "type2"))
    df_temp_full <- dplyr::mutate(df_temp_full2, type = tidyr::replace_na(type, paste0("r<", threshold, ", unassigned")))
    eval(parse(text = paste0("df_temp_full <- dplyr::rename(df_temp_full,", collapse_to_cluster, " = type2.y)")))
    df_temp_full <- dplyr::select(dplyr::ungroup(df_temp_full), -type2)
  }

  if (!is.null(rename_suff)) {
    eval(parse(text = paste0("df_temp_full <- dplyr::rename(df_temp_full, ", paste0(rename_suff, "_type"), " = type, ", paste0(rename_suff, "_r"), " = r)")))
  }
  df_temp_full
}

#' manually change idents as needed
#'
#' @param metadata column of ident
#' @param cluster_col column in metadata containing cluster info
#' @param ident_col column in metadata containing identity assignment
#' @param clusters names of clusters to change, string or vector of strings
#' @param idents new idents to assign, must be length of 1 or same as clusters

#' @export
assign_ident <- function(metadata,
                         cluster_col = "cluster",
                         ident_col = "type",
                         clusters,
                         idents) {
  if (!is.vector(clusters) | !is.vector(idents)) {
    stop("unsupported clusters or idents")
  } else {
    if (length(idents) == 1) {
      idents <- rep(idents, length(clusters))
    } else if (length(idents) != length(clusters)) {
      stop("unsupported lengths pairs of clusters and idents")
    }
  }

  for (n in 1:length(clusters)) {
    mindex <- metadata[[cluster_col]] == clusters[n]
    metadata[mindex, ident_col] <- idents[n]
  }
  metadata
}

#' get top calls for each cluster
#'
#' @param correlation_matrix input similarity matrix
#' @param metadata input metadata with tsne coordinates and cluster ids
#' @param col metadata column, can be cluster or cellid
#' @param collapse_to_cluster if a column name is provided, takes the most frequent call of entire cluster to color in plot
#' @param threshold minimum correlation coefficent cutoff for calling clusters
#' @param topn number of calls for each cluster
#'
#' @export
cor_to_call_topn <- function(correlation_matrix,
                             metadata = NULL,
                             col = "cluster",
                             collapse_to_cluster = FALSE,
                             threshold = 0,
                             topn = 2) {
  df_temp <- tibble::as_tibble(correlation_matrix, rownames = col)
  df_temp <- tidyr::gather(df_temp, key = type, value = r, -!!col)
  df_temp[["type"]][df_temp$r < threshold] <- paste0("r<", threshold, ", unassigned")
  df_temp <- dplyr::top_n(dplyr::group_by_at(df_temp, 1), topn, r)
  df_temp_full <- df_temp

  if (collapse_to_cluster != FALSE) {
    if (!(col %in% colnames(metadata))) {
      metadata <- tibble::as_tibble(metadata, rownames = col)
    }
    df_temp_full <- dplyr::left_join(df_temp_full, metadata, by = col)
    df_temp_full[, "type2"] <- df_temp_full[[collapse_to_cluster]]
    df_temp_full2 <- dplyr::group_by(df_temp_full, type, type2)
    df_temp_full2 <- dplyr::summarize(df_temp_full2, sum = sum(r), n = n())
    df_temp_full2 <- dplyr::group_by(df_temp_full2, type2)
    df_temp_full2 <- dplyr::arrange(df_temp_full2, desc(n), desc(sum))
    df_temp_full2 <- dplyr::filter(df_temp_full2, type != paste0("r<", threshold, ", unassigned"))
    df_temp_full2 <- dplyr::slice(df_temp_full2, 1:topn)
    df_temp_full2 <- dplyr::right_join(df_temp_full2, dplyr::select(df_temp_full, -c(type, r)), by = stats::setNames(collapse_to_cluster, "type2"))
    df_temp_full <- dplyr::mutate(df_temp_full2, type = tidyr::replace_na(type, paste0("r<", threshold, ", unassigned")))
    df_temp_full <- dplyr::group_by_(df_temp_full, .dots = col)
    df_temp_full <- dplyr::distinct(df_temp_full, type, type2, .keep_all = TRUE)
    dplyr::arrange(df_temp_full, desc(n), desc(sum), .by_group = TRUE)
  } else {
    df_temp_full <- dplyr::group_by_(df_temp_full, .dots = col)
    dplyr::arrange(df_temp_full, desc(r), .by_group = TRUE)
  }
}

#' pct of cells in each cluster that express genelist
#'
#' @param matrix expression matrix
#' @param genelist vector of marker genes for one identity
#' @param clusters vector of cluster identities
#' @param returning whether to return mean, min, or max of the gene pct in the gene list
#'
#' @export
gene_pct <- function(matrix,
                     genelist,
                     clusters,
                     returning = "mean") {
  genelist <- intersect(genelist, rownames(matrix))
  clusters[is.na(clusters)] <- "orig.NA"
  unique_clusters <- unique(clusters)


  if (returning == "mean") {
    sapply(unique_clusters, function(x) {
      celllist <- clusters == x
      tmp <- matrix[genelist, celllist, drop = FALSE]
      tmp[tmp > 0] <- 1
      mean(Matrix::rowSums(tmp) / ncol(tmp))
    })
  } else if (returning == "min") {
    sapply(unique_clusters, function(x) {
      celllist <- clusters == x
      tmp <- matrix[genelist, celllist, drop = FALSE]
      tmp[tmp > 0] <- 1
      min(Matrix::rowSums(tmp) / ncol(tmp))
    })
  } else if (returning == "max") {
    sapply(unique_clusters, function(x) {
      celllist <- clusters == x
      tmp <- matrix[genelist, celllist, drop = FALSE]
      tmp[tmp > 0] <- 1
      max(Matrix::rowSums(tmp) / ncol(tmp))
    })
  }
}

#' pct of cells in every cluster that express a series of genelists
#'
#' @param matrix expression matrix
#' @param marker_m matrixized markers
#' @param cluster_info data.frame or vector containing cluster assignments per cell.
#' Order must match column order in supplied matrix. If a data.frame
#' provide the cluster_col parameters.
#' @param cluster_col column in cluster_info with cluster number
#' @param norm whether and how the results are normalized
#'
#' @export
gene_pct_markerm <- function(matrix,
                             marker_m,
                             cluster_info,
                             cluster_col = NULL,
                             norm = NULL) {
  if (is.vector(cluster_info)) {
  } else if (is.data.frame(cluster_info) & !is.null(cluster_col)) {
    cluster_info <- cluster_info[[cluster_col]]
  } else {
    stop("cluster_info not formatted correctly,
         supply either a  vector or a dataframe")
  }

  # coerce factors in character
  if (is.factor(cluster_info)) {
    cluster_info <- as.character(cluster_info)
  }

  if (class(marker_m) != "data.frame") {
    marker_m <- as.data.frame(marker_m)
  }

  out <- sapply(colnames(marker_m), function(x) {
    gene_pct(
      matrix,
      marker_m[[x]],
      cluster_info
    )
  })

  if (!(is.null(norm))) {
    if (norm == "divide") {
      out <- sweep(out, 2, apply(out, 2, max), "/")
    } else if (norm == "diff") {
      out <- sweep(out, 2, apply(out, 2, max), "-")
    } else {
      out <- sweep(out, 2, apply(out, 2, max) * norm)
      out[out < 0] <- 0
      out[out > 0] <- 1
    }
  }

  # edge cases where all markers can't be found for a cluster
  out[is.na(out)] <- 0
  out
}

#' Combined function to compare scRNA-seq data to bulk RNA-seq data and marker list
#'
#' @export
clustify_nudge <- function(input, ...) {
  UseMethod("clustify_nudge", input)
}

#' @rdname clustify_nudge
#' @param input seurat 2 object
#' @param ref_mat reference expression matrix
#' @param marker matrix of markers
#' @param query_genes A vector of genes of interest to compare. If NULL, then common genes between
#' the expr_mat and ref_mat will be used for comparision.
#' @param cluster_col column in metadata that contains cluster ids per cell. Will default to first
#' column of metadata if not supplied. Not required if running correlation per cell.
#' @param compute_method method(s) for computing similarity scores
#' @param weight relative weight for the gene list scores, when added to correlation score
#' @param dr stored dimension reduction
#' @param seurat_out output cor matrix or called seurat object
#' @param threshold identity calling minimum score threshold
#' @param norm whether and how the results are normalized
#' @param marker_inmatrix whether markers genes are already in preprocessed matrix form
#' @param mode use marker expression pct or ranked cor score for nudging
#' @param ... passed to matrixize_markers

#' @export
clustify_nudge.seurat <- function(input,
                                  ref_mat,
                                  marker,
                                  cluster_col = NULL,
                                  query_genes = NULL,
                                  compute_method = "spearman",
                                  weight = 1,
                                  seurat_out = TRUE,
                                  threshold = -Inf,
                                  dr = "tsne",
                                  norm = "diff",
                                  marker_inmatrix = TRUE,
                                  mode = "rank",
                                  ...) {
  if (marker_inmatrix != TRUE) {
    marker <- matrixize_markers(
      marker,
      ...
    )
  }
  resa <- clustify(
    input = input,
    ref_mat = ref_mat,
    cluster_col = cluster_col,
    query_genes = query_genes,
    seurat_out = FALSE,
    per_cell = FALSE
  )

  if (mode == "pct") {
    resb <- gene_pct_markerm(input@data,
      marker,
      input@meta.data,
      cluster_col = cluster_col,
      norm = norm
    )
  } else if (mode == "rank") {
    resb <- pos_neg_select(input@data,
      marker,
      metadata = input@meta.data,
      cluster_col = cluster_col,
      cutoff_score = norm
    )
    empty_vec <- setdiff(colnames(resa), colnames(resb))
    empty_mat <- matrix(0, nrow = nrow(resb), ncol = length(empty_vec), dimnames = list(rownames(resb), empty_vec))
    resb <- cbind(resb, empty_mat)
  }

  df_temp <- cor_to_call(resa[order(rownames(resa)), order(colnames(resa))] +
    resb[order(rownames(resb)), order(colnames(resb))] * weight,
  threshold = threshold
  )
  colnames(df_temp) <- c(cluster_col, "type", "score")

  if (seurat_out == FALSE) {
    df_temp
  } else {
    cluster_info <- as.data.frame(use_seurat_meta(input, dr = dr, seurat3 = FALSE))
    df_temp_full <- dplyr::left_join(tibble::rownames_to_column(cluster_info, "rn"), df_temp, by = cluster_col)
    df_temp_full <- tibble::column_to_rownames(df_temp_full, "rn")
    if ("Seurat" %in% loadedNamespaces()) {
      input@meta.data <- df_temp_full
      return(input)
    } else {
      message("seurat not loaded, returning cor_mat instead")
      return(df_temp)
    }
  }
}

#' @rdname clustify_nudge
#' @param input seurat 2 object
#' @param ref_mat reference expression matrix
#' @param metadata cell cluster assignments, supplied as a vector or data.frame. If
#' data.frame is supplied then `cluster_col` needs to be set.
#' @param marker matrix of markers
#' @param query_genes A vector of genes of interest to compare. If NULL, then common genes between
#' the expr_mat and ref_mat will be used for comparision.
#' @param cluster_col column in metadata that contains cluster ids per cell. Will default to first
#' column of metadata if not supplied. Not required if running correlation per cell.
#' @param compute_method method(s) for computing similarity scores
#' @param ... passed to matrixize_markers
#' @param norm whether and how the results are normalized
#' @param call make call or just return score matrix
#' @param marker_inmatrix whether markers genes are already in preprocessed matrix form
#' @param mode use marker expression pct or ranked cor score for nudging

#' @export
clustify_nudge.default <- function(input,
                                   ref_mat,
                                   marker,
                                   metadata = NULL,
                                   cluster_col = NULL,
                                   query_genes = NULL,
                                   compute_method = "spearman",
                                   weight = 1,
                                   seurat_out = TRUE,
                                   threshold = -Inf,
                                   dr = "tsne",
                                   norm = "diff",
                                   call = TRUE,
                                   marker_inmatrix = TRUE,
                                   mode = "rank",
                                   ...) {
  if (marker_inmatrix != TRUE) {
    marker <- matrixize_markers(
      marker,
      ...
    )
  }

  if (!(stringr::str_detect(class(input), "atrix"))) {
    input_original <- input
    temp <- parse_loc_object(input,
      type = class(input),
      expr_loc = NULL,
      meta_loc = NULL,
      var_loc = NULL,
      cluster_col = cluster_col
    )
    input <- temp[["expr"]]
    metadata <- temp[["meta"]]
    if (is.null(query_genes)) {
      query_genes <- temp[["var"]]
    }
    if (is.null(cluster_col)) {
      cluster_col <- temp[["col"]]
    }
  }

  resa <- clustify(
    input = input,
    ref_mat = ref_mat,
    metadata = metadata,
    cluster_col = cluster_col,
    query_genes = query_genes,
    seurat_out = FALSE,
    per_cell = FALSE
  )

  if (mode == "pct") {
    resb <- gene_pct_markerm(input,
      marker,
      metadata,
      cluster_col = cluster_col,
      norm = norm
    )
  } else if (mode == "rank") {
    resb <- pos_neg_select(input,
      marker,
      metadata,
      cluster_col = cluster_col,
      cutoff_score = norm
    )
    empty_vec <- setdiff(colnames(resa), colnames(resb))
    empty_mat <- matrix(0, nrow = nrow(resb), ncol = length(empty_vec), dimnames = list(rownames(resb), empty_vec))
    resb <- cbind(resb, empty_mat)
  }

  if (call == TRUE) {
    df_temp <- cor_to_call(resa[order(rownames(resa)), order(colnames(resa))] +
      resb[order(rownames(resb)), order(colnames(resb))] * weight,
    threshold = threshold
    )
    colnames(df_temp) <- c(cluster_col, "type", "score")
    df_temp
  } else {
    resa[order(rownames(resa)), order(colnames(resa))] +
      resb[order(rownames(resb)), order(colnames(resb))] * weight
  }
}

#' more flexible parsing of single cell objects
#'
#' @param input input object
#' @param type look up predefined slots/loc
#' @param expr_loc expression matrix location
#' @param meta_loc metadata location
#' @param var_loc variable genes location
#' @param cluster_col column of clustering from metadata
#' @param lookuptable if not supplied, will look in built-in table for object parsing

#' @export
parse_loc_object <- function(input,
                             type = class(input),
                             expr_loc = NULL,
                             meta_loc = NULL,
                             var_loc = NULL,
                             cluster_col = NULL,
                             lookuptable = NULL) {
  # if (type == "SingleCellExperiment") {
  #   parsed = list(input@assays$data$logcounts,
  #                 as.data.frame(input@colData)),
  #                 NULL,
  #                 "cell_type1")
  # }
  #
  # if (type == "URD") {
  #   parsed = list(input@logupx.data,
  #                 input@meta,
  #                 input@var.genes,
  #                 "cluster")
  # }
  #
  # if (type == "FunctionalSingleCellExperiment") {
  #   parsed = list(input@ExperimentList$rnaseq@assays$data$logcounts,
  #                 input@ExperimentList$rnaseq@colData,
  #                 NULL,
  #                 "leiden_cluster")
  # }
  if (is.null(lookuptable)) {
    object_loc_lookup1 <- clustifyr::object_loc_lookup
  } else {
    object_loc_lookup1 <- lookuptable
  }

  if (type %in% colnames(object_loc_lookup1)) {
    parsed <- list(
      eval(parse(text = object_loc_lookup1[[type]][1])),
      as.data.frame(eval(parse(text = object_loc_lookup1[[type]][2]))),
      eval(parse(text = object_loc_lookup1[[type]][3])),
      object_loc_lookup1[[type]][4]
    )
  } else {
    parsed <- list(NULL, NULL, NULL, NULL)
  }

  names(parsed) <- c("expr", "meta", "var", "col")

  if (!(is.null(expr_loc))) {
    parsed[["expr"]] <- eval(parse(text = paste0("input", expr_loc)))
  }

  if (!(is.null(meta_loc))) {
    parsed[["meta"]] <- as.data.frame(eval(parse(text = paste0("input", meta_loc))))
  }

  if (!(is.null(var_loc))) {
    parsed[["var"]] <- eval(parse(text = paste0("input", var_loc)))
  }

  if (!(is.null(cluster_col))) {
    parsed[["col"]] <- cluster_col
  }

  parsed
}

#' compare clustering parameters and classification outcomes
#'
#' @param expr expression matrix
#' @param metadata metadata including cluster info and dimension reduction plotting
#' @param ref_mat reference matrix
#' @param cluster_col column of clustering from metadata
#' @param x_col column of metadata for x axis plotting
#' @param y_col column of metadata for y axis plotting
#' @param n expand n-fold for over/under clustering
#' @param ngenes number of genes to use for feature selection, use all genes if NULL
#' @param query_genes vector, otherwise genes with be recalculated
#' @param do_label whether to label each cluster at median center
#' @param seed set seed for kmeans
#' @param newclustering use kmeans if NULL on dr or col name for second column of clustering

#' @export
overcluster_test <- function(expr,
                             metadata,
                             ref_mat,
                             cluster_col,
                             x_col = "tSNE_1",
                             y_col = "tSNE_2",
                             n = 5,
                             ngenes = NULL,
                             query_genes = NULL,
                             do_label = TRUE,
                             seed = 42,
                             newclustering = NULL) {
  if (!(is.null(seed))) {
    set.seed(seed)
  }

  if (is.null(newclustering)) {
    metadata$new_clusters <- as.character(stats::kmeans(metadata[, c(x_col, y_col)],
      centers = n * length(unique(metadata[[cluster_col]]))
    )$clust)
  } else {
    metadata$new_clusters <- metadata[[newclustering]]
    n <- length(unique(metadata[[newclustering]])) / length(unique(metadata[[cluster_col]]))
  }

  if (is.null(query_genes)) {
    if (is.null(ngenes)) {
      genes <- rownames(expr)
    } else {
      genes <- ref_feature_select(expr, ngenes)
    }
  } else {
    genes <- query_genes
  }
  res1 <- clustify(expr,
    ref_mat,
    metadata,
    query_genes = genes,
    cluster_col = cluster_col,
    seurat_out = FALSE
  )
  res2 <- clustify(expr,
    ref_mat,
    metadata,
    query_genes = genes,
    cluster_col = "new_clusters",
    seurat_out = FALSE
  )
  o1 <- plot_tsne(metadata,
    feature = cluster_col,
    x = x_col,
    y = y_col,
    do_label = FALSE,
    do_legend = FALSE
  )
  o2 <- plot_tsne(metadata,
    feature = "new_clusters",
    x = x_col,
    y = y_col,
    do_label = FALSE,
    do_legend = FALSE
  )
  p1 <- plot_best_call(res1,
    metadata,
    cluster_col,
    do_label = do_label,
    x = x_col,
    y = y_col
  )
  p2 <- plot_best_call(res2,
    metadata,
    "new_clusters",
    do_label = do_label,
    x = x_col,
    y = y_col
  )
  g <- cowplot::plot_grid(o1, o2, p1, p2,
    labels = c(
      length(unique(metadata[[cluster_col]])),
      n * length(unique(metadata[[cluster_col]]))
    )
  )
  return(g)
}

#' feature select from reference matrix
#'
#' @param mat reference matrix
#' @param n number of genes to return
#' @param mode the method of selecting features
#' @param rm.lowvar whether to remove lower variation genes first
#'
#' @export
ref_feature_select <- function(mat,
                               n = 3000,
                               mode = "var",
                               rm.lowvar = TRUE) {
  if (rm.lowvar == TRUE) {
    v <- RowVar(mat)
    v2 <- v[order(-v)][1:(length(v) / 2)]
    mat <- mat[names(v2)[!is.na(names(v2))], ]
  }

  if (mode == "cor") {
    cor_mat <- cor(t(as.matrix(mat)), method = "spearman")
    diag(cor_mat) <- rep(0, times = nrow(cor_mat))
    cor_mat <- abs(cor_mat)
    score <- apply(cor_mat, 1, max, na.rm = TRUE)
    score <- score[order(-score)]
    cor_genes <- names(score[1:n])
  } else if (mode == "var") {
    cor_genes <- names(v2[1:n])
  } else if (mode == "hybrid") {

  }
  cor_genes
}

#' Returns a list of variable genes based on PCA
#'
#' @description  Extract genes, i.e. "features", based on the top
#' loadings of principal components
#' formed from the bulk expression data set
#'
#' @param mat Expression matrix. Rownames are genes,
#' colnames are single cell cluster name, and
#' values are average single cell expression (log transformed).
#' @param pcs Precalculated pcs if available, will skip over processing on mat.
#' @param n_pcs Number of PCs to selected gene loadings from.
#' See the explore_PCA_corr.Rmd vignette for details.
#' @param percentile Select the percentile of absolute values of
#' PCA loadings to select genes from. E.g. 0.999 would select the
#' top point 1 percent of genes with the largest loadings.
#' @param if_log whether the data is already log transformed
#' @return The list of genes to use as features.
#'
#' @export
feature_select_PCA <- function(mat = NULL,
                               pcs = NULL,
                               n_pcs = 10,
                               percentile = 0.99,
                               if_log = TRUE) {
  if (if_log == FALSE) {
    mat <- log(mat + 1)
  }

  # Get the PCs
  if (is.null(pcs)) {
    pca <- prcomp(t(as.matrix(mat)))$rotation
  } else {
    pca <- pcs
  }

  # For the given number PCs, select the genes with the largest loadings
  genes <- c()
  for (i in 1:n_pcs) {
    cutoff <- quantile(abs(pca[, i]), probs = percentile)
    genes <- c(genes, rownames(pca[abs(pca[, i]) >= cutoff, ]))
  }

  return(genes)
}

#' convert gmt format of pathways to list of vectors
#'
#' @param path gmt file path
#' @param cutoff remove pathways with less genes than this cutoff
#' @param sep sep used in file to split path and genes

#' @export
gmt_to_list <- function(path,
                        cutoff = 0,
                        sep = "\thttp://www.broadinstitute.org/gsea/msigdb/cards/.*?\t") {
  df <- readr::read_csv(path,
    col_names = FALSE
  )
  df <- tidyr::separate(df,
    X1,
    sep = sep,
    into = c("path", "genes")
  )
  pathways <- stringr::str_split(
    df$genes,
    "\t"
  )
  names(pathways) <- stringr::str_replace(
    df$path,
    "REACTOME_",
    ""
  )
  if (cutoff > 0) {
    ids <- sapply(pathways, function(i) length(i) < cutoff)
    pathways <- pathways[!ids]
  }
  return(pathways)
}

#' plot GSEA pathway scores as heatmap, returns a list containing results and plot.
#'
#' @param mat expression matrix
#' @param pathway_list a list of vectors, each named for a specific pathway, or dataframe
#' @param n_perm Number of permutation for fgsea function. Defaults to 1000.
#' @param scale convert expr_mat into zscores prior to running GSEA?, default = TRUE
#' @param topn number of top pathways to plot
#' @param returning to return "both" list and plot, or either one

#' @export
plot_pathway_gsea <- function(mat,
                              pathway_list,
                              n_perm = 1000,
                              scale = TRUE,
                              topn = 5,
                              returning = "both") {
  res <- calculate_pathway_gsea(mat,
    pathway_list,
    n_perm,
    scale = scale
  )
  coltopn <- unique(cor_to_call_topn(res, topn = topn, threshold = -Inf)$type)
  res[is.na(res)] <- 0
  g <- ComplexHeatmap::Heatmap(res[, coltopn], column_names_gp = grid::gpar(fontsize = 6))

  if (returning == "both") {
    return(list(res, g))
  } else if (returning == "plot") {
    return(g)
  } else {
    return(res)
  }
}

#' get var per row for matrix
#'
#' @param x expression matrix
#' @param na.rm logical. Should missing values (including NaN) be omitted from the calculations?

#' @export
RowVar <- function(x, na.rm = TRUE) {
  rowSums((x - rowMeans(x, na.rm = na.rm))^2, na.rm = na.rm) / (dim(x)[2] - 1)
}

#' downsample matrix by cluster or completely random
#'
#' @param mat expression matrix
#' @param n number per cluster or fraction to keep
#' @param keep_cluster_proportions whether to subsample
#' @param cluster_info data.frame or vector containing cluster assignments per cell.
#' Order must match column order in supplied matrix. If a data.frame
#' provide the cluster_col parameters.
#' @param cluster_col column in cluster_info with cluster number
#' @param set_seed random seed

#' @export
downsample_matrix <- function(mat,
                              n = 1,
                              keep_cluster_proportions = TRUE,
                              cluster_info = NULL,
                              cluster_col = "cluster",
                              set_seed = NULL) {
  if (keep_cluster_proportions == FALSE) {
    cluster_ids <- colnames(mat)
    if (n < 1) {
      n <- as.integer(ncol(mat) * n)
    }
    set.seed(set_seed)
    cluster_ids_new <- sample(cluster_ids, n)
  } else {
    if (is.vector(cluster_info)) {
      cluster_ids <- split(colnames(mat), cluster_info)
    } else if (is.data.frame(cluster_info) & !is.null(cluster_col)) {
      cluster_ids <- split(colnames(mat), cluster_info[[cluster_col]])
    } else if (class(cluster_info) == "factor") {
      cluster_info <- as.character(cluster_info)
      cluster_ids <- split(colnames(mat), cluster_info)
    } else {
      stop("cluster_info not formatted correctly,
         supply either a  vector or a dataframe")
    }
    if (n < 1) {
      n2 <- sapply(cluster_ids, function(x) as.integer(length(x) * n))
      n <- n2
    }
    set.seed(set_seed)
    cluster_ids_new <- mapply(sample, cluster_ids, n, SIMPLIFY = FALSE)
  }
  return(mat[, unlist(cluster_ids_new)])
}

#' make combination ref matrix to assess intermixing
#'
#' @param ref_mat reference expression matrix
#' @param if_log whether input data is natural
#' @param sep separator for name combinations
#'
#' @export
make_comb_ref <- function(ref_mat, if_log = TRUE, sep = "_and_") {
  if (if_log == TRUE) {
    ref_mat <- expm1(ref_mat)
  }
  combs <- utils::combn(x = colnames(ref_mat), m = 2, simplify = FALSE)
  comb_mat <- sapply(combs, FUN = function(x) Matrix::rowMeans(ref_mat[, unlist(x)]))
  colnames(comb_mat) <- sapply(combs, FUN = function(x) stringr::str_c(unlist(x), collapse = sep))
  new_mat <- cbind(ref_mat, comb_mat)
  if (if_log == TRUE) {
    new_mat <- log1p(new_mat)
  }
  new_mat
}

#' marker selection from reference matrix
#'
#' @param mat reference matrix
#' @param cut an expression minimum cutoff
#' @param arrange whether to arrange (lower means better)
#' @param compto compare max expression to the value of next 1 or more
#'
#' @export
ref_marker_select <- function(mat, cut = 0.5, arrange = TRUE, compto = 1) {
  mat <- mat[!is.na(rownames(mat)), ]
  mat <- mat[Matrix::rowSums(mat) != 0, ]
  ref_cols <- colnames(mat)
  res <- apply(mat, 1, marker_select, ref_cols, cut, compto = compto)
  if (class(res) == "list") {
    res <- res[!sapply(res, is.null)]
  }
  resdf <- t(as.data.frame(res, stringsAsFactors = FALSE))
  resdf <- tibble::rownames_to_column(as.data.frame(resdf, stringsAsFactors = FALSE), "gene")
  colnames(resdf) <- c("gene", "cluster", "ratio")
  resdf <- dplyr::mutate(resdf, ratio = as.numeric(ratio))
  if (arrange == TRUE) {
    resdf <- dplyr::group_by(resdf, cluster)
    resdf <- dplyr::arrange(resdf, ratio, .by_group = TRUE)
    resdf <- dplyr::ungroup(resdf)
  }
  resdf
}

#' decide for one gene whether it is a marker for a certain cell type
#' @param row1 a numeric vector of expression values (row)
#' @param cols a vector of cell types (column)
#' @param cut an expression minimum cutoff
#' @param compto compare max expression to the value of next 1 or more
#'
#' @export
marker_select <- function(row1, cols, cut = 1, compto = 1) {
  row_sorted <- sort(row1, decreasing = TRUE)
  col_sorted <- names(row_sorted)
  num_sorted <- unname(row_sorted)
  if (num_sorted[1] >= cut) {
    return(c(col_sorted[1], (num_sorted[1 + compto] / num_sorted[1])))
  }
}

#' adapt clustify to tweak score for pos and neg markers
#' @param input single-cell expression matrix
#' @param metadata cell cluster assignments, supplied as a vector or data.frame. If
#' data.frame is supplied then `cluster_col` needs to be set. Not required if running correlation per cell.
#' @param ref_mat reference expression matrix with positive and negative markers(set expression at 0)
#' @param cluster_col column in metadata that contains cluster ids per cell. Will default to first
#' column of metadata if not supplied. Not required if running correlation per cell.
#' @param cutoff_n expression cutoff where genes ranked below n are considered non-expressing
#' @param cutoff_score positive score lower than this cutoff will be considered as 0 to not influence scores
#' @param ... additional arguments to pass to compute_method function
#' @export
pos_neg_select <- function(input,
                           ref_mat,
                           metadata,
                           cluster_col = "cluster",
                           cutoff_n = 0,
                           cutoff_score = 0.5) {
  suppressWarnings(res <- clustify(rbind(input, "clustifyr0" = 0.01),
    ref_mat,
    metadata,
    cluster_col = cluster_col,
    per_cell = TRUE, verbose = TRUE
  ))
  res[is.na(res)] <- 0
  suppressWarnings(res2 <- average_clusters(t(res),
    metadata,
    cluster_col = cluster_col,
    if_log = FALSE,
    output_log = FALSE
  ))
  res2 <- t(res2)

  if (!(is.null(cutoff_score))) {
    res2 <- apply(res2, 2, function(x) {
      maxr <- max(x)
      if (maxr > 0.1) {
        x[x > 0 & x < cutoff_score * maxr] <- 0
      }
      x
    })
  }

  res2
}

#' generate negative markers from a list of exclusive positive markers
#' @param mat matrix or dataframe of markers
#' @export
reverse_marker_matrix <- function(mat) {
  full_vec <- as.vector(t(mat))
  mat_rev <- apply(mat, 2, function(x) {
    full_vec[!(full_vec %in% x)]
  })
  as.data.frame(mat_rev)
}

#' takes files with positive and negative markers, as described in garnett, and returns list of markers
#' @param filename txt file to load
#' @export
file_marker_parse <- function(filename) {
  lines <- readLines(filename)
  count <- 0
  ident_names <- c()
  ident_pos <- c()
  ident_neg <- c()
  for (line in lines) {
    tag <- substr(line, 1, 1)
    if (tag == ">") {
      count <- count + 1
      ident_names[count] <- substr(line, 2, nchar(line))
    } else if (tag == "e") {
      ident_pos[count] <- strsplit(substr(line, 12, nchar(line)), split = ", ")
    } else if (tag == "n") {
      ident_neg[count] <- strsplit(substr(line, 16, nchar(line)), split = ", ")
    }
  }
  names(ident_neg) <- ident_names
  names(ident_pos) <- ident_names
  list("pos" = ident_pos, "neg" = ident_neg)
}
