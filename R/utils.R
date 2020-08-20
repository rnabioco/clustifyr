#' Overcluster by kmeans per cluster
#'
#' @param mat expression matrix
#' @param cluster_id list of ids per cluster
#' @param power decides the number of clusters for kmeans
#' @return new cluster_id list of more clusters
#' @examples
#' res <- overcluster(
#'     mat = pbmc_matrix_small,
#'     cluster_id = split(colnames(pbmc_matrix_small), pbmc_meta$classified)
#' )
#' length(res)
#' @export
overcluster <- function(mat,
                        cluster_id,
                        power = 0.15) {
  mat <- as.matrix(mat)
  new_ids <- list()
  for (name in names(cluster_id)) {
    ids <- cluster_id[[name]]
    if (length(ids) > 1) {
      new_clusters <-
        stats::kmeans(t(mat[, ids]),
                      centers = as.integer(length(ids)^power))
      new_ids1 <-
        split(names(new_clusters$cluster),
              new_clusters$cluster)
      names(new_ids1) <-
        stringr::str_c(name, names(new_ids1), sep = "_")
      new_ids <- append(new_ids, new_ids1)
    } else {
      new_ids <- append(new_ids, cluster_id[name])
    }
  }
  new_ids
}

#' Average expression values per cluster
#'
#' @param mat expression matrix
#' @param metadata data.frame or vector containing cluster assignments per cell.
#' Order must match column order in supplied matrix. If a data.frame
#' provide the cluster_col parameters.
#' @param if_log input data is natural log,
#' averaging will be done on unlogged data
#' @param cluster_col column in metadata with cluster number
#' @param cell_col if provided, will reorder matrix first
#' @param low_threshold option to remove clusters with too few cells
#' @param method whether to take mean (default) or median
#' @param output_log whether to report log results
#' @param subclusterpower whether to get multiple averages per original cluster
#' @param cut_n set on a limit of genes as expressed, lower ranked genes
#' are set to 0, considered unexpressed
#' @return average expression matrix, with genes for row names, and clusters
#'  for column names
#' @examples
#' mat <- average_clusters(
#'     mat = pbmc_matrix_small,
#'     metadata = pbmc_meta,
#'     cluster_col = "classified",
#'     if_log = FALSE
#' )
#' mat[1:3, 1:3]
#' @importFrom matrixStats rowMaxs rowMedians colRanks
#' @export
average_clusters <- function(mat,
                             metadata,
                             cluster_col = "cluster",
                             if_log = TRUE,
                             cell_col = NULL,
                             low_threshold = 0,
                             method = "mean",
                             output_log = TRUE,
                             subclusterpower = 0,
                             cut_n = NULL) {
  cluster_info <- metadata
  if (!(is.null(cell_col))) {
    if (!(all(colnames(mat) == cluster_info[[cell_col]]))) {
      mat <- mat[, cluster_info[[cell_col]]]
    }
  }
  if (is.vector(cluster_info)) {
    if (ncol(mat) != length(cluster_info)) {
      stop("vector of cluster assignments does not match the number of columns in the matrix",
           call. = FALSE)
    }
    cluster_ids <- split(colnames(mat), cluster_info)
  } else if (is.data.frame(cluster_info) & !is.null(cluster_col)) {
    cluster_info_temp <- cluster_info[[cluster_col]]
    if (is.factor(cluster_info_temp)) {
      cluster_info_temp <- droplevels(cluster_info_temp)
    }
    cluster_ids <- split(colnames(mat), cluster_info_temp)
  } else if (is.factor(cluster_info)) {
    cluster_info <- as.character(cluster_info)
    if (ncol(mat) != length(cluster_info)) {
      stop("vector of cluster assignments does not match the number of columns in the matrix",
           call. = FALSE)
    }
    cluster_ids <- split(colnames(mat), cluster_info)
  } else {
    stop("metadata not formatted correctly,
         supply either a  vector or a dataframe",
         call. = FALSE
    )
  }
  
  if (subclusterpower > 0) {
    cluster_ids <-
      overcluster(mat, cluster_ids, power = subclusterpower)
  }
  
  if (method == "mean") {
    out <- lapply(
      cluster_ids,
      function(cell_ids) {
        if (!all(cell_ids %in% colnames(mat))) {
          stop("cell ids not found in input matrix",
               call. = FALSE)
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
          stop("cell ids not found in input matrix",
               call. = FALSE)
        }
        mat_data <- mat[, cell_ids, drop = FALSE]
        mat_data[mat_data == 0] <- NA
        res <- matrixStats::rowMedians(as.matrix(mat_data),
                                       na.rm = TRUE)
        res[is.na(res)] <- 0
        names(res) <- rownames(mat_data)
        res
      }
    )
  }
  
  out <- do.call(cbind, out)
  if (low_threshold > 0) {
    fil <- vapply(cluster_ids,
                  FUN = length,
                  FUN.VALUE = numeric(1)) >= low_threshold
    if (!all(as.vector(fil))) {
      message("The following clusters have less than ", low_threshold, " cells for this analysis: ",
              paste(colnames(out)[!as.vector(fil)], collapse = ", "),
              ". They are excluded.")
    }
    out <- out[, as.vector(fil)]
  } else {
    fil <- vapply(cluster_ids,
                  FUN = length,
                  FUN.VALUE = numeric(1)) >= 10
    if (!all(as.vector(fil))) {
      message("The following clusters have less than ", 10, " cells for this analysis: ",
              paste(colnames(out)[!as.vector(fil)], collapse = ", "),
              ". Classification is likely inaccurate.")
    }
  }
  if (!(is.null(cut_n))) {
    expr_mat <- out
    expr_df <- as.matrix(expr_mat)
    df_temp <- t(matrixStats::colRanks(-expr_df,
                                       ties.method = "average"))
    rownames(df_temp) <- rownames(expr_mat)
    colnames(df_temp) <- colnames(expr_mat)
    expr_mat[df_temp > cut_n] <- 0
    out <- expr_mat
  }
  
  return(out)
}

#' Percentage detected per cluster
#'
#' @param mat expression matrix
#' @param metadata data.frame with cells
#' @param cluster_col column in metadata with cluster number
#' @param cut_num binary cutoff for detection
#' @return matrix of numeric values, with genes for row names,
#' and clusters for column names
percent_clusters <- function(mat,
                             metadata,
                             cluster_col = "cluster",
                             cut_num = 0.5) {
  cluster_info <- metadata
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
#' @return matrix of 1s and 0s
get_best_match_matrix <- function(cor_mat) {
  cor_mat <- as.matrix(cor_mat)
  best_mat <-
    as.data.frame(cor_mat - matrixStats::rowMaxs(as.matrix(cor_mat)))
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
#' @return string with ident call and possibly cor value
get_best_str <- function(name,
                         best_mat,
                         cor_mat,
                         carry_cor = TRUE) {
  if (sum(as.numeric(best_mat[name, ])) > 0) {
    best.names <- colnames(best_mat)[which(best_mat[name, ] == 1)]
    best.cor <-
      round(cor_mat[name, which(best_mat[name, ] == 1)], 2)
    for (i in seq_len(length(best.cor))) {
      if (i == 1) {
        str <- paste0(best.names[i],
                      " (",
                      best.cor[i],
                      ")")
      } else {
        str <- paste0(str,
                      "; ",
                      best.names[i],
                      " (",
                      best.cor[i],
                      ")")
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
#' @description return entries found in all supplied vectors.
#'  If the vector supplied is NULL or NA, then it will be excluded
#'  from the comparision.
#' @param ... vectors
#' @return vector of shared elements
get_common_elements <- function(...) {
  vecs <- list(...)
  # drop NULL elements of list
  vecs <- vecs[!vapply(vecs, is.null, FUN.VALUE = logical(1))]
  # drop NA elements of list (NA values OK in a vector)
  vecs <- vecs[!is.na(vecs)]
  
  Reduce(intersect, vecs)
}

#' Convert expression matrix to GSEA pathway scores
#' (would take a similar place in workflow before average_clusters/binarize)
#'
#' @param mat expression matrix
#' @param pathway_list a list of vectors, each named for a specific pathway,
#' or dataframe
#' @param n_perm Number of permutation for fgsea function. Defaults to 1000.
#' @param scale convert expr_mat into zscores prior to running GSEA?,
#' default = FALSE
#' @param no_warnings suppress warnings from gsea ties
#' @return matrix of GSEA NES values, cell types as row names,
#' pathways as column names
#' @examples
#' gl <- list(
#'     "n" = c("PPBP", "LYZ", "S100A9"),
#'     "a" = c("IGLL5", "GNLY", "FTL")
#' )
#'
#' pbmc_avg <- average_clusters(
#'     mat = pbmc_matrix_small,
#'     metadata = pbmc_meta,
#'     cluster_col = "classified"
#' )
#'
#' calculate_pathway_gsea(
#'     mat = pbmc_avg,
#'     pathway_list = gl
#' )
#' @export
calculate_pathway_gsea <- function(mat,
                                   pathway_list,
                                   n_perm = 1000,
                                   scale = TRUE,
                                   no_warnings = TRUE) {
  # pathway_list can be user defined or
  #`my_pathways <- fgsea::reactomePathways(rownames(pbmc4k_matrix))`
  out <- lapply(
    names(pathway_list),
    function(y) {
      marker_list <- list()
      marker_list[[1]] <- pathway_list[[y]]
      names(marker_list) <- y
      v1 <- marker_list
      temp <- run_gsea(
        mat,
        v1,
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

#' manually change idents as needed
#'
#' @param metadata column of ident
#' @param cluster_col column in metadata containing cluster info
#' @param ident_col column in metadata containing identity assignment
#' @param clusters names of clusters to change, string or
#'  vector of strings
#' @param idents new idents to assign, must be length of 1 or
#' same as clusters
#' @return new dataframe of metadata
assign_ident <- function(metadata,
                         cluster_col = "cluster",
                         ident_col = "type",
                         clusters,
                         idents) {
  if (!is.vector(clusters) | !is.vector(idents)) {
    stop("unsupported clusters or idents", call. = FALSE)
  } else {
    if (length(idents) == 1) {
      idents <- rep(idents, length(clusters))
    } else if (length(idents) != length(clusters)) {
      stop("unsupported lengths pairs of clusters and idents",
           call. = FALSE)
    }
  }
  
  for (n in seq_len(length(clusters))) {
    mindex <- metadata[[cluster_col]] == clusters[n]
    metadata[mindex, ident_col] <- idents[n]
  }
  metadata
}

#' get top calls for each cluster
#'
#' @param cor_mat input similarity matrix
#' @param metadata input metadata with tsne or umap coordinates
#' and cluster ids
#' @param col metadata column, can be cluster or cellid
#' @param collapse_to_cluster if a column name is provided,
#' takes the most frequent call of entire cluster to color in plot
#' @param threshold minimum correlation coefficent cutoff for calling clusters
#' @param topn number of calls for each cluster
#' @return dataframe of cluster, new potential ident, and r info
#' @examples
#' res <- clustify(
#'     input = pbmc_matrix_small,
#'     metadata = pbmc_meta,
#'     ref_mat = cbmc_ref,
#'     query_genes = pbmc_vargenes,
#'     cluster_col = "classified"
#' )
#'
#' cor_to_call_topn(
#'     cor_mat = res,
#'     metadata = pbmc_meta,
#'     col = "classified",
#'     collapse_to_cluster = FALSE,
#'     threshold = 0.5
#' )
#' @export
cor_to_call_topn <- function(cor_mat,
                             metadata = NULL,
                             col = "cluster",
                             collapse_to_cluster = FALSE,
                             threshold = 0,
                             topn = 2) {
  correlation_matrix <- cor_mat
  df_temp <- tibble::as_tibble(correlation_matrix, rownames = col)
  df_temp <-
    tidyr::gather(
      df_temp,
      key = !!dplyr::sym("type"),
      value = !!dplyr::sym("r"),
      -!!col
    )
  df_temp[["type"]][df_temp$r < threshold] <-
    paste0("r<", threshold, ", unassigned")
  df_temp <-
    dplyr::top_n(dplyr::group_by_at(df_temp, 1),
                 topn,
                 !!dplyr::sym("r"))
  df_temp_full <- df_temp
  
  if (collapse_to_cluster != FALSE) {
    if (!(col %in% colnames(metadata))) {
      metadata <- tibble::as_tibble(metadata, rownames = col)
    }
    df_temp_full <-
      dplyr::left_join(df_temp_full, metadata, by = col)
    df_temp_full[, "type2"] <- df_temp_full[[collapse_to_cluster]]
    df_temp_full2 <-
      dplyr::group_by(
        df_temp_full,
        !!dplyr::sym("type"),
        !!dplyr::sym("type2")
      )
    df_temp_full2 <-
      dplyr::summarize(df_temp_full2,
                       sum = sum(!!dplyr::sym("r")),
                       n = n()
      )
    df_temp_full2 <-
      dplyr::group_by(df_temp_full2, !!dplyr::sym("type2"))
    df_temp_full2 <-
      dplyr::arrange(df_temp_full2, desc(n), desc(sum))
    df_temp_full2 <-
      dplyr::filter(
        df_temp_full2,
        !!dplyr::sym("type") != paste0("r<",
                                       threshold,
                                       ", unassigned")
      )
    df_temp_full2 <- dplyr::slice(df_temp_full2, seq_len(topn))
    df_temp_full2 <-
      dplyr::right_join(
        df_temp_full2,
        dplyr::select(df_temp_full, -c(
          !!dplyr::sym("type"), !!dplyr::sym("r")
        )),
        by = stats::setNames(collapse_to_cluster, "type2")
      )
    df_temp_full <-
      dplyr::mutate(df_temp_full2, type = tidyr::replace_na(
        !!dplyr::sym("type"),
        paste0("r<", threshold, ", unassigned")
      ))
    df_temp_full <- dplyr::group_by(df_temp_full,
                                    !!dplyr::sym(col))
    df_temp_full <-
      dplyr::distinct(df_temp_full,
                      !!dplyr::sym("type"),
                      !!dplyr::sym("type2"),
                      .keep_all = TRUE
      )
    dplyr::arrange(df_temp_full, desc(n), desc(sum),
                   .by_group = TRUE)
  } else {
    df_temp_full <- dplyr::group_by(df_temp_full,
                                    !!dplyr::sym(col))
    dplyr::arrange(df_temp_full, desc(!!dplyr::sym("r")),
                   .by_group = TRUE)
  }
}

#' pct of cells in each cluster that express genelist
#'
#' @param matrix expression matrix
#' @param genelist vector of marker genes for one identity
#' @param clusters vector of cluster identities
#' @param returning whether to return mean, min,
#' or max of the gene pct in the gene list
#' @return vector of numeric values
gene_pct <- function(matrix,
                     genelist,
                     clusters,
                     returning = "mean") {
  genelist <- intersect(genelist, rownames(matrix))
  if (is.factor(clusters)) {
    clusters <-
      factor(clusters, levels = c(levels(clusters), "orig.NA"))
  }
  clusters[is.na(clusters)] <- "orig.NA"
  unique_clusters <- unique(clusters)
  
  if (returning == "mean") {
    vapply(unique_clusters, function(x) {
      celllist <- clusters == x
      tmp <- matrix[genelist, celllist, drop = FALSE]
      tmp[tmp > 0] <- 1
      mean(Matrix::rowSums(tmp) / ncol(tmp))
    }, FUN.VALUE = numeric(1))
  } else if (returning == "min") {
    vapply(unique_clusters, function(x) {
      celllist <- clusters == x
      tmp <- matrix[genelist, celllist, drop = FALSE]
      tmp[tmp > 0] <- 1
      min(Matrix::rowSums(tmp) / ncol(tmp))
    }, FUN.VALUE = numeric(1))
  } else if (returning == "max") {
    vapply(unique_clusters, function(x) {
      celllist <- clusters == x
      tmp <- matrix[genelist, celllist, drop = FALSE]
      tmp[tmp > 0] <- 1
      max(Matrix::rowSums(tmp) / ncol(tmp))
    }, FUN.VALUE = numeric(1))
  }
}

#' pct of cells in every cluster that express a series of genelists
#'
#' @param matrix expression matrix
#' @param marker_m matrixized markers
#' @param metadata data.frame or vector containing cluster
#' assignments per cell.
#' Order must match column order in supplied matrix. If a data.frame
#' provide the cluster_col parameters.
#' @param cluster_col column in metadata with cluster number
#' @param norm whether and how the results are normalized
#' @return matrix of numeric values, clusters from mat as row names,
#'  cell types from marker_m as column names
#' @examples
#' gene_pct_markerm(
#'     matrix = pbmc_matrix_small,
#'     marker_m = cbmc_m,
#'     metadata = pbmc_meta,
#'     cluster_col = "classified"
#' )
#' @export
gene_pct_markerm <- function(matrix,
                             marker_m,
                             metadata,
                             cluster_col = NULL,
                             norm = NULL) {
  cluster_info <- metadata
  if (is.vector(cluster_info)) {
    
  } else if (is.data.frame(cluster_info) & !is.null(cluster_col)) {
    cluster_info <- cluster_info[[cluster_col]]
  } else {
    stop("metadata not formatted correctly,
         supply either a  vector or a dataframe",
         call. = FALSE
    )
  }
  
  # coerce factors in character
  if (is.factor(cluster_info)) {
    cluster_info <- as.character(cluster_info)
  }
  
  if (!is.data.frame(marker_m)) {
    marker_m <- as.data.frame(marker_m)
  }
  
  out <- vapply(colnames(marker_m), function(x) {
    gene_pct(
      matrix,
      marker_m[[x]],
      cluster_info
    )
  }, FUN.VALUE = numeric(length(unique(cluster_info))))
  
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

#' Combined function to compare scRNA-seq data to
#'  bulk RNA-seq data and marker list
#'
#' @examples
#' # Seurat2
#' clustify_nudge(
#'     input = s_small,
#'     ref_mat = cbmc_ref,
#'     marker = cbmc_m,
#'     cluster_col = "res.1",
#'     threshold = 0.8,
#'     seurat_out = FALSE,
#'     mode = "pct",
#'     dr = "tsne"
#' )
#'
#' # Seurat3
#' clustify_nudge(
#'     input = s_small3,
#'     ref_mat = cbmc_ref,
#'     marker = cbmc_m,
#'     cluster_col = "RNA_snn_res.1",
#'     threshold = 0.8,
#'     seurat_out = FALSE,
#'     mode = "pct",
#'     dr = "tsne"
#' )
#'
#' # Matrix
#' clustify_nudge(
#'     input = pbmc_matrix_small,
#'     ref_mat = cbmc_ref,
#'     metadata = pbmc_meta,
#'     marker = as.matrix(cbmc_m),
#'     query_genes = pbmc_vargenes,
#'     cluster_col = "classified",
#'     threshold = 0.8,
#'     call = FALSE,
#'     marker_inmatrix = FALSE,
#'     mode = "pct"
#' )
#' @export
clustify_nudge <- function(input, ...) {
  UseMethod("clustify_nudge", input)
}

#' @rdname clustify_nudge
#' @param input express matrix or object
#' @param ref_mat reference expression matrix
#' @param metadata cell cluster assignments, supplied as a vector
#' or data.frame. If
#' data.frame is supplied then `cluster_col` needs to be set.
#' @param marker matrix of markers
#' @param query_genes A vector of genes of interest to compare.
#' If NULL, then common genes between
#' the expr_mat and ref_mat will be used for comparision.
#' @param cluster_col column in metadata that contains cluster ids per cell.
#'  Will default to first
#' column of metadata if not supplied.
#' Not required if running correlation per cell.
#' @param compute_method method(s) for computing similarity scores
#' @param weight relative weight for the gene list scores,
#' when added to correlation score
#' @param dr stored dimension reduction
#' @param ... passed to matrixize_markers
#' @param norm whether and how the results are normalized
#' @param call make call or just return score matrix
#' @param marker_inmatrix whether markers genes are already
#'  in preprocessed matrix form
#' @param mode use marker expression pct or ranked cor score for nudging
#' @param obj_out whether to output object instead of cor matrix
#' @param seurat_out output cor matrix or called seurat object
#' @param rename_prefix prefix to add to type and r column names
#' @param lookuptable if not supplied, will look in built-in
#' table for object parsing
#' @param threshold identity calling minimum score threshold,
#'  only used when obj_out = T

#' @return single cell object, or matrix of numeric values,
#'  clusters from input as row names, cell types from ref_mat as column names
#' @export
clustify_nudge.default <- function(input,
                                   ref_mat,
                                   marker,
                                   metadata = NULL,
                                   cluster_col = NULL,
                                   query_genes = NULL,
                                   compute_method = "spearman",
                                   weight = 1,
                                   seurat_out = FALSE,
                                   threshold = -Inf,
                                   dr = "umap",
                                   norm = "diff",
                                   call = TRUE,
                                   marker_inmatrix = TRUE,
                                   mode = "rank",
                                   obj_out = FALSE,
                                   rename_prefix = NULL,
                                   lookuptable = NULL,
                                   ...) {
  if (marker_inmatrix != TRUE) {
    marker <- matrixize_markers(
      marker,
      ...
    )
  }
  
  if (!inherits(input, c("matrix", "Matrix", "data.frame"))) {
    input_original <- input
    temp <- parse_loc_object(
      input,
      type = class(input),
      expr_loc = NULL,
      meta_loc = NULL,
      var_loc = NULL,
      cluster_col = cluster_col,
      lookuptable = lookuptable
    )
    
    if (!(is.null(temp[["expr"]]))) {
      message("recognized object type - ", class(input))
    }
    
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
    if (ncol(marker) > 1 && is.character(marker[1, 1])) {
      marker <- pos_neg_marker(marker)
    }
    resb <- pos_neg_select(input,
                           marker,
                           metadata,
                           cluster_col = cluster_col,
                           cutoff_score = NULL
    )
    empty_vec <- setdiff(colnames(resa), colnames(resb))
    empty_mat <-
      matrix(
        0,
        nrow = nrow(resb),
        ncol = length(empty_vec),
        dimnames = list(rownames(resb), empty_vec)
      )
    resb <- cbind(resb, empty_mat)
  }
  
  res <- resa[order(rownames(resa)), order(colnames(resa))] +
    resb[order(rownames(resb)), order(colnames(resb))] * weight
  
  if ((obj_out ||
       seurat_out) &&
      !inherits(input_original, c("matrix", "Matrix", "data.frame"))) {
    df_temp <- cor_to_call(
      res,
      metadata = metadata,
      cluster_col = cluster_col,
      threshold = threshold
    )
    
    df_temp_full <- call_to_metadata(
      df_temp,
      metadata = metadata,
      cluster_col = cluster_col,
      per_cell = FALSE,
      rename_prefix = rename_prefix
    )
    
    out <- insert_meta_object(input_original,
                              df_temp_full,
                              lookuptable = lookuptable
    )
    
    return(out)
  } else {
    if (call == TRUE) {
      df_temp <- cor_to_call(res,
                             threshold = threshold
      )
      colnames(df_temp) <- c(cluster_col, "type", "score")
      return(df_temp)
    } else {
      return(res)
    }
  }
}

#' @rdname clustify_nudge
#' @export
clustify_nudge.seurat <- function(input,
                                  ref_mat,
                                  marker,
                                  cluster_col = NULL,
                                  query_genes = NULL,
                                  compute_method = "spearman",
                                  weight = 1,
                                  seurat_out = TRUE,
                                  obj_out = FALSE,
                                  threshold = -Inf,
                                  dr = "umap",
                                  norm = "diff",
                                  marker_inmatrix = TRUE,
                                  mode = "rank",
                                  rename_prefix = NULL,
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
    per_cell = FALSE,
    dr = dr
  )
  
  if (mode == "pct") {
    resb <- gene_pct_markerm(
      object_data(input, "data"),
      marker,
      object_data(input, "meta.data"),
      cluster_col = cluster_col,
      norm = norm
    )
  } else if (mode == "rank") {
    if (ncol(marker) > 1 && is.character(marker[1, 1])) {
      marker <- pos_neg_marker(marker)
    }
    resb <- pos_neg_select(
      object_data(input, "data"),
      marker,
      metadata = object_data(input, "meta.data"),
      cluster_col = cluster_col,
      cutoff_score = NULL
    )
    empty_vec <- setdiff(colnames(resa), colnames(resb))
    empty_mat <-
      matrix(
        0,
        nrow = nrow(resb),
        ncol = length(empty_vec),
        dimnames = list(rownames(resb), empty_vec)
      )
    resb <- cbind(resb, empty_mat)
  }
  
  res <- resa[order(rownames(resa)), order(colnames(resa))] +
    resb[order(rownames(resb)), order(colnames(resb))] * weight
  
  if (!(seurat_out || obj_out)) {
    res
  } else {
    df_temp <- cor_to_call(
      res,
      object_data(input, "meta.data"),
      cluster_col = cluster_col,
      threshold = threshold
    )
    
    df_temp_full <- call_to_metadata(
      df_temp,
      object_data(input, "meta.data"),
      cluster_col = cluster_col,
      per_cell = FALSE,
      rename_prefix = rename_prefix
    )
    
    if ("Seurat" %in% loadedNamespaces()) {
      input <- write_meta(input, df_temp_full)
      return(input)
    } else {
      message("seurat not loaded, returning cor_mat instead")
      return(res)
    }
    input
  }
}

#' @rdname clustify_nudge
#' @export
clustify_nudge.Seurat <- function(input,
                                  ref_mat,
                                  marker,
                                  cluster_col = NULL,
                                  query_genes = NULL,
                                  compute_method = "spearman",
                                  weight = 1,
                                  seurat_out = TRUE,
                                  obj_out = FALSE,
                                  threshold = -Inf,
                                  dr = "umap",
                                  norm = "diff",
                                  marker_inmatrix = TRUE,
                                  mode = "rank",
                                  rename_prefix = NULL,
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
    per_cell = FALSE,
    dr = dr
  )
  
  if (mode == "pct") {
    resb <- gene_pct_markerm(
      object_data(input, "data"),
      marker,
      object_data(input, "meta.data"),
      cluster_col = cluster_col,
      norm = norm
    )
  } else if (mode == "rank") {
    if (ncol(marker) > 1 && is.character(marker[1, 1])) {
      marker <- pos_neg_marker(marker)
    }
    resb <- pos_neg_select(
      object_data(input, "data"),
      marker,
      object_data(input, "meta.data"),
      cluster_col = cluster_col,
      cutoff_score = NULL
    )
    empty_vec <- setdiff(colnames(resa), colnames(resb))
    empty_mat <-
      matrix(
        0,
        nrow = nrow(resb),
        ncol = length(empty_vec),
        dimnames = list(rownames(resb), empty_vec)
      )
    resb <- cbind(resb, empty_mat)
  }
  
  res <- resa[order(rownames(resa)), order(colnames(resa))] +
    resb[order(rownames(resb)), order(colnames(resb))] * weight
  
  if (!(seurat_out || obj_out)) {
    res
  } else {
    df_temp <- cor_to_call(
      res,
      metadata = object_data(input, "meta.data"),
      cluster_col = cluster_col,
      threshold = threshold
    )
    
    df_temp_full <- call_to_metadata(
      df_temp,
      metadata = object_data(input, "meta.data"),
      cluster_col = cluster_col,
      per_cell = FALSE,
      rename_prefix = rename_prefix
    )
    
    if ("Seurat" %in% loadedNamespaces()) {
      input <- write_meta(input, df_temp_full)
      return(input)
    } else {
      message("seurat not loaded, returning cor_mat instead")
      return(res)
    }
    input
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
#' @return list of expression, metadata, vargenes, cluster_col info from object
#' @examples
#' obj <- parse_loc_object(s_small3)
#' length(obj)
#' @export
parse_loc_object <- function(input,
                             type = class(input),
                             expr_loc = NULL,
                             meta_loc = NULL,
                             var_loc = NULL,
                             cluster_col = NULL,
                             lookuptable = NULL) {
  if (is.null(lookuptable)) {
    object_loc_lookup1 <- clustifyr::object_loc_lookup
  } else {
    object_loc_lookup1 <- lookuptable
  }
  
  if (length(intersect(type, colnames(object_loc_lookup1))) > 0) {
    type <- intersect(type, colnames(object_loc_lookup1))[1]
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
    parsed[["meta"]] <-
      as.data.frame(eval(parse(text = paste0(
        "input", meta_loc
      ))))
  }
  
  if (!(is.null(var_loc))) {
    parsed[["var"]] <- eval(parse(text = paste0("input", var_loc)))
  }
  
  if (!(is.null(cluster_col))) {
    parsed[["col"]] <- cluster_col
  }
  
  parsed
}

#' more flexible metadata update of single cell objects
#'
#' @param input input object
#' @param new_meta new metadata table to insert back into object
#' @param type look up predefined slots/loc
#' @param meta_loc metadata location
#' @param lookuptable if not supplied,
#' will look in built-in table for object parsing
#' @return new object with new metadata inserted
#' @examples
#' \dontrun{
#' insert_meta_object(s_small3, seurat_meta(s_small3, dr = "tsne"))
#' }
#' @export
insert_meta_object <- function(input,
                               new_meta,
                               type = class(input),
                               meta_loc = NULL,
                               lookuptable = NULL) {
  if (is.null(lookuptable)) {
    object_loc_lookup1 <- clustifyr::object_loc_lookup
  } else {
    object_loc_lookup1 <- lookuptable
  }
  
  if (!type %in% colnames(object_loc_lookup1)) {
    stop("unrecognized object type", call. = FALSE)
  } else {
    text1 <- paste0(object_loc_lookup1[[type]][2], " <- ", "new_meta")
    eval(parse(text = text1))
    return(input)
  }
}

#' compare clustering parameters and classification outcomes
#'
#' @param expr expression matrix
#' @param metadata metadata including cluster info and
#' dimension reduction plotting
#' @param ref_mat reference matrix
#' @param cluster_col column of clustering from metadata
#' @param x_col column of metadata for x axis plotting
#' @param y_col column of metadata for y axis plotting
#' @param n expand n-fold for over/under clustering
#' @param ngenes number of genes to use for feature selection,
#' use all genes if NULL
#' @param query_genes vector, otherwise genes with be recalculated
#' @param do_label whether to label each cluster at median center
#' @param do_legend whether to draw legend
#' @param newclustering use kmeans if NULL on dr
#' or col name for second column of clustering
#' @param threshold type calling threshold
#' @param combine if TRUE return a single plot with combined panels, if
#' FALSE return list of plots (default: TRUE)
#' @return faceted ggplot object
#' @examples
#' set.seed(42)
#' overcluster_test(
#'     expr = pbmc_matrix_small,
#'     metadata = pbmc_meta,
#'     ref_mat = cbmc_ref,
#'     cluster_col = "classified",
#'     x_col = "UMAP_1",
#'     y_col = "UMAP_2"
#' )
#' @export
overcluster_test <- function(expr,
                             metadata,
                             ref_mat,
                             cluster_col,
                             x_col = "UMAP_1",
                             y_col = "UMAP_2",
                             n = 5,
                             ngenes = NULL,
                             query_genes = NULL,
                             threshold = 0,
                             do_label = TRUE,
                             do_legend = FALSE,
                             newclustering = NULL,
                             combine = TRUE) {
  if (is.null(newclustering)) {
    metadata$new_clusters <-
      as.character(stats::kmeans(metadata[, c(x_col, y_col)],
                                 centers = n * length(unique(metadata[[cluster_col]]))
      )$clust)
  } else {
    metadata$new_clusters <- metadata[[newclustering]]
    n <-  length(unique(metadata[[newclustering]])) /
      length(unique(metadata[[cluster_col]]))
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
  res1 <- clustify(
    expr,
    ref_mat,
    metadata,
    query_genes = genes,
    cluster_col = cluster_col,
    seurat_out = FALSE
  )
  res2 <- clustify(
    expr,
    ref_mat,
    metadata,
    query_genes = genes,
    cluster_col = "new_clusters",
    seurat_out = FALSE
  )
  o1 <- plot_dims(
    metadata,
    feature = cluster_col,
    x = x_col,
    y = y_col,
    do_label = FALSE,
    do_legend = FALSE
  )
  o2 <- plot_dims(
    metadata,
    feature = "new_clusters",
    x = x_col,
    y = y_col,
    do_label = FALSE,
    do_legend = FALSE
  )
  p1 <- plot_best_call(
    res1,
    metadata,
    cluster_col,
    threshold = threshold,
    do_label = do_label,
    do_legend = do_legend,
    x = x_col,
    y = y_col
  )
  p2 <- plot_best_call(
    res2,
    metadata,
    "new_clusters",
    threshold = threshold,
    do_label = do_label,
    do_legend = do_legend,
    x = x_col,
    y = y_col
  )
  n_orig_clusters <- length(unique(metadata[[cluster_col]]))
  n_new_clusters <- n * length(unique(metadata[[cluster_col]]))
  
  if(combine){
    g <- suppressWarnings(cowplot::plot_grid(o1, o2, p1, p2,
                                             labels = c(
                                               n_orig_clusters,
                                               n_new_clusters
                                             ))
    )
  } else {
    g <- list(original_clusters = o1,
              new_clusters = o2,
              original_cell_types = p1,
              new_cell_types = p2)
  }
  
  return(g)
}

#' feature select from reference matrix
#'
#' @param mat reference matrix
#' @param n number of genes to return
#' @param mode the method of selecting features
#' @param rm.lowvar whether to remove lower variation genes first
#' @return vector of genes
#' @examples
#' pbmc_avg <- average_clusters(
#'     mat = pbmc_matrix_small,
#'     metadata = pbmc_meta,
#'     cluster_col = "classified"
#' )
#'
#' ref_feature_select(
#'     mat = pbmc_avg[1:100, ],
#'     n = 5
#' )
#' @export
ref_feature_select <- function(mat,
                               n = 3000,
                               mode = "var",
                               rm.lowvar = TRUE) {
  if (rm.lowvar == TRUE) {
    if (!(is.matrix(mat))) {
      mat <- as.matrix(mat)
    }
    v <- matrixStats::rowVars(mat)
    names(v) <- rownames(mat)
    v2 <- v[order(-v)][seq_len(length(v) / 2)]
    mat <- mat[names(v2)[!is.na(names(v2))], ]
  }
  
  if (mode == "cor") {
    cor_mat <- cor(t(as.matrix(mat)), method = "spearman")
    diag(cor_mat) <- rep(0, times = nrow(cor_mat))
    cor_mat <- abs(cor_mat)
    score <- matrixStats::rowMaxs(cor_mat, na.rm = TRUE)
    names(score) <- rownames(cor_mat)
    score <- score[order(-score)]
    cor_genes <- names(score[seq_len(n)])
  } else if (mode == "var") {
    cor_genes <- names(v2[seq_len(n)])
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
#' @return vector of genes
#' @examples
#' feature_select_PCA(
#'     cbmc_ref,
#'     if_log = FALSE
#' )
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
  for (i in seq_len(n_pcs)) {
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
#' @return list of genes in each pathway
#' @examples
#' gmt_file <- system.file(
#'     "extdata",
#'     "c2.cp.reactome.v6.2.symbols.gmt.gz",
#'     package = "clustifyr"
#' )
#'
#' gl <- gmt_to_list(path = gmt_file)
#' length(gl)
#' @export
gmt_to_list <- function(path,
                        cutoff = 0,
                        sep = "\thttp://www.broadinstitute.org/gsea/msigdb/cards/.*?\t") {
  df <- readr::read_csv(path,
                        col_names = FALSE
  )
  df <- tidyr::separate(df, !!dplyr::sym("X1"),
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
    ids <- vapply(pathways, function(i) {
      length(i) < cutoff
    }, FUN.VALUE = logical(1))
    pathways <- pathways[!ids]
  }
  return(pathways)
}

#' plot GSEA pathway scores as heatmap,
#'  returns a list containing results and plot.
#'
#' @param mat expression matrix
#' @param pathway_list a list of vectors, each named for a specific pathway,
#' or dataframe
#' @param n_perm Number of permutation for fgsea function. Defaults to 1000.
#' @param scale convert expr_mat into zscores prior to running GSEA?,
#'  default = TRUE
#' @param topn number of top pathways to plot
#' @param returning to return "both" list and plot, or either one
#' @return list of matrix and plot, or just plot, matrix of GSEA NES values,
#' cell types as row names, pathways as column names
#' @examples
#' gl <- list(
#'     "n" = c("PPBP", "LYZ", "S100A9"),
#'     "a" = c("IGLL5", "GNLY", "FTL")
#' )
#'
#' pbmc_avg <- average_clusters(
#'     mat = pbmc_matrix_small,
#'     metadata = pbmc_meta,
#'     cluster_col = "classified"
#' )
#'
#' plot_pathway_gsea(
#'     pbmc_avg,
#'     gl,
#'     5
#' )
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
  coltopn <-
    unique(cor_to_call_topn(res, topn = topn, threshold = -Inf)$type)
  res[is.na(res)] <- 0
  
  g <- suppressWarnings(ComplexHeatmap::Heatmap(res[, coltopn],
                                                column_names_gp = grid::gpar(fontsize = 6)
  ))
  
  if (returning == "both") {
    return(list(res, g))
  } else if (returning == "plot") {
    return(g)
  } else {
    return(res)
  }
}

#' downsample matrix by cluster or completely random
#'
#' @param mat expression matrix
#' @param n number per cluster or fraction to keep
#' @param keep_cluster_proportions whether to subsample
#' @param metadata data.frame or
#' vector containing cluster assignments per cell.
#' Order must match column order in supplied matrix. If a data.frame
#' provide the cluster_col parameters.
#' @param cluster_col column in metadata with cluster number
#' @return new smaller mat with less cell_id columns
#' @examples
#' set.seed(42)
#' mat <- downsample_matrix(
#'     mat = pbmc_matrix_small,
#'     metadata = pbmc_meta$classified,
#'     n = 10,
#'     keep_cluster_proportions = TRUE
#' )
#' mat[1:3, 1:3]
#' @export
downsample_matrix <- function(mat,
                              n = 1,
                              keep_cluster_proportions = TRUE,
                              metadata = NULL,
                              cluster_col = "cluster") {
  cluster_info <- metadata
  if (keep_cluster_proportions == FALSE) {
    cluster_ids <- colnames(mat)
    if (n < 1) {
      n <- as.integer(ncol(mat) * n)
    }
    cluster_ids_new <- sample(cluster_ids, n)
  } else {
    if (is.vector(cluster_info)) {
      cluster_ids <- split(colnames(mat), cluster_info)
    } else if (is.data.frame(cluster_info) &
               !is.null(cluster_col)) {
      cluster_ids <- split(colnames(mat), cluster_info[[cluster_col]])
    } else if (is.factor(cluster_info)) {
      cluster_info <- as.character(cluster_info)
      cluster_ids <- split(colnames(mat), cluster_info)
    } else {
      stop("metadata not formatted correctly,
         supply either a  vector or a dataframe",
           call. = FALSE
      )
    }
    if (n < 1) {
      n2 <- vapply(cluster_ids, function(x) {
        as.integer(length(x) * n)
      }, FUN.VALUE = numeric(1))
      n <- n2
    }
    cluster_ids_new <-
      mapply(sample, cluster_ids, n, SIMPLIFY = FALSE)
  }
  return(mat[, unlist(cluster_ids_new)])
}

#' marker selection from reference matrix
#'
#' @param mat reference matrix
#' @param cut an expression minimum cutoff
#' @param arrange whether to arrange (lower means better)
#' @param compto compare max expression to the value of next 1 or more
#' @return dataframe, with gene, cluster, ratio columns
#' @examples
#' ref_marker_select(
#'     cbmc_ref,
#'     cut = 2
#' )
#' @export
ref_marker_select <-
  function(mat,
           cut = 0.5,
           arrange = TRUE,
           compto = 1) {
    mat <- mat[!is.na(rownames(mat)), ]
    mat <- mat[Matrix::rowSums(mat) != 0, ]
    ref_cols <- colnames(mat)
    res <-
      apply(mat, 1, marker_select, ref_cols, cut, compto = compto)
    if (is.list(res)) {
      res <- res[!vapply(res, is.null, FUN.VALUE = logical(1))]
    }
    resdf <- t(as.data.frame(res, stringsAsFactors = FALSE))
    if (tibble::has_rownames(as.data.frame(resdf,
                                           stringsAsFactors = FALSE))) {
      resdf <- tibble::remove_rownames(as.data.frame(resdf,
                                                     stringsAsFactors = FALSE))
    }
    resdf <- tibble::rownames_to_column(resdf,
                                        "gene")
    colnames(resdf) <- c("gene", "cluster", "ratio")
    resdf <-
      dplyr::mutate(resdf,
                    ratio = as.numeric(!!dplyr::sym("ratio")))
    if (arrange == TRUE) {
      resdf <- dplyr::group_by(resdf, cluster)
      resdf <-
        dplyr::arrange(resdf, !!dplyr::sym("ratio"),
                       .by_group = TRUE)
      resdf <- dplyr::ungroup(resdf)
    }
    resdf
  }

#' decide for one gene whether it is a marker for a certain cell type
#' @param row1 a numeric vector of expression values (row)
#' @param cols a vector of cell types (column)
#' @param cut an expression minimum cutoff
#' @param compto compare max expression to the value of next 1 or more
#' @return vector of cluster name and ratio value
#' @examples
#' pbmc_avg <- average_clusters(
#'     mat = pbmc_matrix_small,
#'     metadata = pbmc_meta,
#'     cluster_col = "classified",
#'     if_log = FALSE
#' )
#'
#' marker_select(
#'     row1 = pbmc_avg["PPBP", ],
#'     cols = names(pbmc_avg["PPBP", ])
#' )
#' @export
marker_select <- function(row1,
                          cols,
                          cut = 1,
                          compto = 1) {
  row_sorted <- sort(row1, decreasing = TRUE)
  col_sorted <- names(row_sorted)
  num_sorted <- unname(row_sorted)
  if (num_sorted[1] >= cut) {
    return(c(col_sorted[1], (num_sorted[1 + compto] / num_sorted[1])))
  }
}

#' adapt clustify to tweak score for pos and neg markers
#' @param input single-cell expression matrix
#' @param metadata cell cluster assignments,
#' supplied as a vector or data.frame. If
#' data.frame is supplied then `cluster_col` needs to be set.
#'  Not required if running correlation per cell.
#' @param ref_mat reference expression matrix with positive and
#' negative markers(set expression at 0)
#' @param cluster_col column in metadata that contains cluster ids per cell.
#' Will default to first
#' column of metadata if not supplied.
#' Not required if running correlation per cell.
#' @param cutoff_n expression cutoff where genes ranked below n are
#'  considered non-expressing
#' @param cutoff_score positive score lower than this cutoff will be
#' considered as 0 to not influence scores
#' @return matrix of numeric values, clusters from input as row names,
#'  cell types from ref_mat as column names
#' @examples
#' pn_ref <- data.frame(
#'     "Myeloid" = c(1, 0.01, 0),
#'     row.names = c("CD74", "clustifyr0", "CD79A")
#' )
#'
#' pos_neg_select(
#'     input = pbmc_matrix_small,
#'     ref_mat = pn_ref,
#'     metadata = pbmc_meta,
#'     cluster_col = "classified",
#'     cutoff_score = 0.8
#' )
#' @export
pos_neg_select <- function(input,
                           ref_mat,
                           metadata,
                           cluster_col = "cluster",
                           cutoff_n = 0,
                           cutoff_score = 0.5) {
  suppressWarnings(
    res <- clustify(
      rbind(input, "clustifyr0" = 0.01),
      ref_mat,
      metadata,
      cluster_col = cluster_col,
      per_cell = TRUE,
      verbose = TRUE,
      query_genes = rownames(ref_mat)
    )
  )
  res[is.na(res)] <- 0
  
  suppressWarnings(
    res2 <- average_clusters(
      t(res),
      metadata,
      cluster_col = cluster_col,
      if_log = FALSE,
      output_log = FALSE
    )
  )
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
#' @return matrix of gene names
#' @examples
#' reverse_marker_matrix(cbmc_m)
#' @export
reverse_marker_matrix <- function(mat) {
  full_vec <- as.vector(t(mat))
  mat_rev <- apply(mat, 2, function(x) {
    full_vec[!(full_vec %in% x)]
  })
  as.data.frame(mat_rev)
}

#' generate pos and negative marker expression matrix from a
#' list/dataframe of positive markers
#' @param mat matrix or dataframe of markers
#' @return matrix of gene expression
#' @examples
#' m1 <- pos_neg_marker(cbmc_m)
#' @export
pos_neg_marker <- function(mat) {
  if (is.data.frame(mat)) {
    mat <- as.list(mat)
  } else if (is.matrix(mat)) {
    mat <- as.list(as.data.frame(mat,
                                 stringsAsFactors = FALSE))
  } else if (!is.list(mat)) {
    stop("unsupported marker format,
             must be dataframe, matrix, or list",
         call. = FALSE
    )
  }
  genelist <- mat
  typenames <- names(genelist)
  
  g2 <- lapply(genelist, function(x) {
    data.frame(gene = x, stringsAsFactors = FALSE)
  })
  
  g2 <- dplyr::bind_rows(g2, .id = "type")
  g2 <- dplyr::mutate(g2, expression = 1)
  g2 <- tidyr::spread(g2, key = "type", value = "expression")
  if (tibble::has_rownames(g2)) {
    g2 <- tibble::remove_rownames(g2)
  }
  g2 <- tibble::column_to_rownames(g2, "gene")
  g2[is.na(g2)] <- 0
  g2
}
#' takes files with positive and negative markers, as described in garnett,
#' and returns list of markers
#' @param filename txt file to load
#' @return list of positive and negative gene markers
#' @examples
#' marker_file <- system.file(
#'     "extdata",
#'     "hsPBMC_markers.txt",
#'     package = "clustifyr"
#' )
#'
#' file_marker_parse(marker_file)
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
      ident_pos[count] <-
        strsplit(substr(line, 12, nchar(line)), split = ", ")
    } else if (tag == "n") {
      ident_neg[count] <-
        strsplit(substr(line, 16, nchar(line)), split = ", ")
    }
  }
  
  if (!(is.null(ident_neg))) {
    names(ident_neg) <- ident_names
  }
  if (!(is.null(ident_pos))) {
    names(ident_pos) <- ident_names
  }
  list("pos" = ident_pos, "neg" = ident_neg)
}

#' Generate a unique column id for a dataframe
#' @param df dataframe with column names
#' @param id desired id if unique
#' @return character
get_unique_column <- function(df, id = NULL) {
  if (!is.null(id)) {
    out_id <- id
  } else {
    out_id <- "x"
  }
  
  res <- ifelse(out_id %in% colnames(df),
                make.unique(c(colnames(df),
                              out_id))[length(c(colnames(df),
                                                out_id))],
                out_id
  )
  
  res
}

#' Find rank bias
#' @param mat original query expression matrix
#' @param metadata metadata after calling types
#' @param type_col column name in metadata that contains ids
#' @param ref_mat reference expression matrix
#' @param query_genes original vector of genes used to clustify
#' @param filter_out whether to only report filtered results
#' @param expr_cut consider all lower expressing genes as off
#' @param threshold diff threshold
#' @param consensus_cut filter out if lower
#' han number of types show large diff
#' @return matrix of rank diff values
#' @examples
#' res <- clustify(
#'     input = pbmc_matrix_small,
#'     metadata = pbmc_meta,
#'     ref_mat = cbmc_ref,
#'     query_genes = pbmc_vargenes,
#'     cluster_col = "classified"
#' )
#' call1 <- cor_to_call(
#'     res,
#'     metadata = pbmc_meta,
#'     cluster_col = "classified",
#'     collapse_to_cluster = FALSE,
#'     threshold = 0.8
#' )
#' pbmc_meta2 <- call_to_metadata(
#'     call1,
#'     pbmc_meta,
#'     "classified"
#' )
#' find_rank_bias(
#'     pbmc_matrix_small,
#'     pbmc_meta2, "type",
#'     cbmc_ref,
#'     query_genes = pbmc_vargenes
#' )
#' @export
find_rank_bias <- function(mat,
                           metadata,
                           type_col,
                           ref_mat,
                           query_genes = NULL,
                           filter_out = TRUE,
                           threshold = 0.33,
                           expr_cut = 3000,
                           consensus_cut = 1) {
  if (is.null(query_genes)) {
    query_genes <- intersect(
      rownames(mat),
      rownames(ref_mat)
    )
  } else {
    query_genes <- intersect(
      query_genes,
      intersect(
        rownames(mat),
        rownames(ref_mat)
      )
    )
  }
  avg2 <- average_clusters(
    mat[, rownames(metadata)],
    metadata[[type_col]]
  )
  r2 <- t(matrixStats::colRanks(-avg2[query_genes, ],
                                ties.method = "average"))
  rownames(r2) <- query_genes
  colnames(r2) <- colnames(avg2)
  r2 <- r2[, colnames(r2)[!stringr::str_detect(colnames(r2),
                                               "unassigned"),
                          drop = FALSE],
           drop = FALSE]
  
  r1 <- t(matrixStats::colRanks(-ref_mat[query_genes, ],
                                ties.method = "average"))
  rownames(r1) <- query_genes
  colnames(r1) <- colnames(ref_mat)
  r1 <- r1[, colnames(r2), drop = FALSE]
  
  if (!(is.null(expr_cut))) {
    r1[r1 > expr_cut] <- expr_cut
    r2[r2 > expr_cut] <- expr_cut
    nthreshold <- expr_cut * threshold
  } else {
    nthreshold <- length(query_genes) * threshold
  }
  rdiff <- r1 - r2
  if (filter_out) {
    rp <- rdiff > nthreshold | rdiff < -nthreshold
    rp[r1 > 0.9 * expr_cut & r2 > 0.9 * expr_cut] <- NA
    v <- rowMeans(rp, na.rm = TRUE) == 1
    v[is.na(v)] <- FALSE
    
    v2 <- Matrix::rowSums(rp, na.rm = TRUE) == 1
    prob <- rdiff[v & !v2, , drop = FALSE]
    
    return(prob)
  } else {
    return(rdiff)
  }
}

#' Given a reference matrix and a list of genes, take the union of 
#' all genes in vector and genes in reference matrix 
#' and insert zero counts for all remaining genes. 
#' @param gene_vector char vector with gene names
#' @param ref_matrix Reference matrix containing cell types vs. 
#' gene expression values
#' @return Reference matrix with union of all genes
#' @examples
#' mat <- append_genes(
#' gene_vector = human_genes_10x,
#' ref_matrix = cbmc_ref 
#' ) 
#' @export
append_genes <- function(gene_vector, ref_matrix)
{
  rownamesGSEMatrix <- rownames(ref_matrix) 
  #Get rownames from GSEMatrix (new GSE file)
  
  rowCountHumanGenes <- length(gene_vector) 
  #Calculate number of rows from list of full human genes
  rowCountNewGSEFile <- nrow(ref_matrix) 
  #Calculate number of rows of GSE matrix
  
  missing_rows <- setdiff(gene_vector, rownamesGSEMatrix) 
  #Use setdiff function to figure out rows which are different/missing 
  #from GSE matrix
  
  zeroExpressionMatrix <- matrix(
    0, 
    nrow = length(missing_rows), 
    ncol = ncol(ref_matrix)) 
  #Create a placeholder matrix with zeroes and missing_rows length
  
  rownames(zeroExpressionMatrix) <- missing_rows 
  #Assign row names
  colnames(zeroExpressionMatrix) <- colnames(ref_matrix) 
  #Assign column names
  
  full_matrix <- rbind(ref_matrix, zeroExpressionMatrix) 
  #Bind GSEMatrix and zeroExpressionMatrix together
  
  #Reorder matrix
  full_matrix <- full_matrix[gene_vector, ] 
  #Reorder fullMatrix to preserve gene order
  return(full_matrix) 
  #Return fullMatrix
}

#' Given a count matrix, determine if the matrix has been either 
#' log-normalized, normalized, or contains raw counts
#' @param counts_matrix Count matrix containing scRNA-seq read data
#' @param max_log_value Static value to determine if a matrix is normalized
#' @return String either raw counts, log-normalized or normalized
#' @examples
#' check_raw_counts(pbmc_matrix_small)
#' @export
check_raw_counts <- function(counts_matrix, max_log_value = 50)
{
  if (is(counts_matrix, 'sparseMatrix')) {
    counts_matrix <- as.matrix(counts_matrix)
  }
  if(!is.matrix(counts_matrix))
  {
    counts_matrix <- as.matrix(counts_matrix)
  }
  if (is.integer(counts_matrix))
  {
    return("raw counts")
  }
  else if (is.double(counts_matrix))
  {
    if (all(counts_matrix == floor(counts_matrix)))
    {
      return("raw counts")
    }
    if(max(counts_matrix) > max_log_value)
    {
      return("normalized")
    }
    else if (min(counts_matrix) < 0)
    {
      stop("negative values detected, likely scaled data")
    }
    else
    {
      return("log-normalized")
    }
  }
  else
  {
    stop("unknown matrix format: ", typeof(counts_matrix))
  }
}

#' Function to combine records into single atlas
#'
#' @param matrix_fns character vector of paths to study matrices stored as .rds files.
#' If a named character vector, then the name will be added as a suffix to the cell type
#' name in the final matrix. If it is not named, then the filename will be used (without .rds)
#' @param genes_fn text file with a single column containing genes and the ordering desired
#' in the output matrix
#' @param matrix_objs Checks to see whether .rds files will be read or R objects in a local environment 
#' @param output_fn output filename for .rds file. If NULL the matrix will be returned instead of
#' saving
#' @return Combined matrix with all genes given
#' @export
build_atlas <- function(matrix_fns = NULL,
                        genes_fn,
                        matrix_objs = NULL,
                        output_fn = NULL)
{
  genesVector <- readr::read_lines(genes_fn)
  if(is.null(matrix_obs) && !is.null(matrix_fns))
  {
    ref_mats <- lapply(matrix_fns, readRDS)
    if(is.null(names(matrix_fns)))
    {
      names(ref_mats) <- basename(ref_matrices_fns) %>% stringr::str_remove(".rds$")
    } 
    else 
    {
      names(ref_mats) <- names(matrix_fns)
    }
  } 
  else if(!is.null(matrix_obs)) 
  {
    ref_mats <- matrix_obs
    if(is.null(names(matrix_obs)))
    {
      names(ref_mats) <- basename(ref_matrices_fns) %>% stringr::str_remove(".rds$")
    } 
    else 
    {
      names(ref_mats) <- names(matrix_obs)
    }
  }
  new_mats <- list()
  for(i in seq_along(ref_mats))
  {
    # standardize genes in matrix
    mat <- appendGenes(geneVector = genesVector,
                       GSEMatrix = as.matrix(ref_mats[[i]]))
    # get study name
    mat_name <- names(ref_mats)[i]
    
    # append study name to cell type names
    new_cols <- paste0(colnames(mat),
                       " (",
                       mat_name,
                       ")")
    colnames(mat) <- new_cols
    
    # assign to list
    new_mats[[i]] <- mat
  }
  
  # cbind a list of matrices
  atlas <- do.call(cbind, new_mats)
  
  if(!is.null(output_fn))
  {
    saveRDS(atlas, output_fn)
  } 
  else 
  {
    return(atlas)
  }
}