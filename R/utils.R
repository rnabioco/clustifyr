#' Average expression values per cluster
#'
#' @param mat expression matrix
#' @param cluster_info data.frame or vector containing cluster assignments per cell.
#' Order must match column order in supplied matrix. If a data.frame
#' provide the cluster_col parameters.
#' @param log_scale input data is natural log,
#' averaging will be done on unlogged data
#' @param cluster_col column in cluster_info with cluster number
#'
#' @export
average_clusters <- function(mat, cluster_info,
                             log_scale = T,
                             cluster_col = "cluster") {

  if(is.vector(cluster_info)){
    cluster_ids <- split(colnames(mat), cluster_info)
  } else if (is.data.frame(cluster_info) & !is.null(cluster_col)){
    cluster_ids <- split(colnames(mat), cluster_info[[cluster_col]])
  } else {
    stop("cluster_info not formatted correctly,
         supply either a  vector or a dataframe")
  }

  out <- lapply(
    cluster_ids,
    function(cell_ids) {
      if (!all(cell_ids %in% colnames(mat))) {
        stop("cell ids not found in input matrix")
      }
      if (log_scale) {
        mat_data <- expm1(mat[, cell_ids, drop = FALSE])
      } else {
        mat_data <- mat[, cell_ids, drop = FALSE]
      }
      res <- Matrix::rowMeans(mat_data)
      if (log_scale) {
        res <- log1p(res)
      }
      res
    }
  )
  return(do.call(cbind, out))
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
    log_scale = F,
    cluster_col = cluster_col
  )
}

#' Function to make best call from correlation matrix
#'
#' @param cor_mat correlation matrix
#'
#' @export
get_best_match_matrix <- function(correlation_matrix) {

  best_mat <- as.data.frame(t(apply(correlation_matrix, 1, function(x) x - max(x))))
  best_mat[best_mat==0]="1"
  best_mat[best_mat!="1"]="0"

  return (best_mat)
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
#' @param num_perm number of permutations, set to 0 by default
#' @param compute_method method(s) for computing similarity scores
#' @param use_var_genes if providing a seurat object, use the variable genes
#'  (stored in seurat_object@var.genes) as the query_genes.
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
                           ...){
  row_ref <- (metadata[[sample_col]] == sample_id)
  expr_mat_ref <- expr_mat[,row_ref]
  expr_mat_tar <- expr_mat[,!row_ref]
  meta_ref <- metadata[row_ref,]
  meta_tar <- metadata[!row_ref,]

  avg_clusters_ref <- average_clusters(expr_mat_ref, meta_ref,
                                       log_scale = F,
                                       cluster_col = cluster_col)

  r2 <- clustify(expr_mat_tar, avg_clusters_ref, meta_tar,
                 query_genes = query_genes,
                 cluster_col = cluster_col,
                 per_cell = per_cell,
                 num_perm = 0,
                 compute_method = compute_method,
                 use_var_genes = FALSE)

  r2
}

#' Average expression values per cluster, filtered by set parameter, defaults to calculating background
#'
#' @param mat expression matrix
#' @param cluster_info data.frame or vector containing cluster assignments per cell, and attribute to filter on.
#' Order must match column order in supplied matrix. If a data.frame
#' provide the cluster_col parameters.
#' @param log_scale input data is natural log,
#' averaging will be done on unlogged data
#' @param filter_on column in cluster_info to filter on
#' @param filter_method "<", "==", ">" compared to filter_value
#' @param filter_value
#'
#' @export
average_clusters_filter <- function(mat, cluster_info,
                                   log_scale = T,
                                   filter_on = "nGene",
                                   group_by = NULL,
                                   filter_method = "<=",
                                   filter_value = 300) {
  eval(parse(text = paste0("cell_ids <- cluster_info[[filter_on]] ", filter_method, "filter_value")))
  if (sum(cell_ids) == 0) {
    stop("no cells kept after filtering")
  }

  if (!is.null(group_by)) {
    res <- average_clusters(mat, cluster_info,
                     log_scale = log_scale,
                     cluster_col = group_by)
  } else {
    if (log_scale) {
      mat_data <- expm1(mat[, cell_ids, drop = FALSE])
    } else {
      mat_data <- mat[, cell_ids, drop = FALSE]
    }
    res <- Matrix::rowMeans(mat_data)
    if (log_scale) {
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

remove_background <- function(mat, background, n = 0){
  if (n == 0) {
    n = length(background)
  }

  if (!is.vector(background)) {
    background <- background[order(background[,1], decreasing = T), , drop = F]
    background <- rownames(t3)[1:n]
  } else if (!is.null(names(background))) {
    background <- names(sort(background, decreasing = T)[1:n])
  }

  mat[!(rownames(mat) %in% background), ]
}

#' Convert expression matrix to GSEA pathway scores (would take a similar place in workflow before average_clusters/binarize)
#'
#' @param mat expression matrix
#' @param pathway_list a list of vectors, each named for a specific pathway, or dataframe
#' @param n_perm Number of permutation for fgsea function. Defaults to 1000.
#' @param scale convert expr_mat into zscores prior to running GSEA?, default = FALSE

#' @export

calculate_pathway_gsea <- function(mat,
                                   pathway_list,
                                   n_perm = 1000,
                                   scale = FALSE) {
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
                       per_cell = T)
      temp <- temp[ , 3, drop = F]
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
#'
#' @export
cor_to_call <- function(correlation_matrix,
                        metadata = NULL,
                        col = "cluster",
                        collapse_to_cluster = FALSE,
                        threshold = 0) {
  df_temp <- tibble::as_tibble(correlation_matrix, rownames = col)
  df_temp <- tidyr::gather(df_temp, key = type, value = r, -!!col)
  df_temp[["type"]][df_temp$r < threshold] <- paste0("r<", threshold,", unassigned")
  df_temp <- dplyr::top_n(dplyr::group_by_at(df_temp, 1), 1, r)
  if (nrow(df_temp) != nrow(correlation_matrix)) {
    clash <- dplyr::summarize(dplyr::group_by_at(df_temp, 1), n = n())
    clash <- dplyr::filter(clash, n > 1)
    clash <- dplyr::pull(clash, 1)
    df_temp[lapply(df_temp[,1], FUN = function(x) x %in% clash)[[1]],2] <- paste0(df_temp[["type"]][lapply(df_temp[,1], FUN = function(x) x %in% clash)[[1]]], "-CLASH!")
    df_temp <- dplyr::distinct(df_temp, exclude = "type", .keep_all = T)
    df_temp_full <- dplyr::select(df_temp, -exclude)
  } else {
    df_temp_full <- df_temp
  }

  if(collapse_to_cluster != FALSE){
    if (!(col %in% colnames(metadata))) {
      metadata <- tibble::as_tibble(metadata, rownames = col)
    }
    df_temp_full <- dplyr::left_join(df_temp_full, metadata, by = col)
    df_temp_full[,"type2"] <- df_temp_full[[collapse_to_cluster]]
    df_temp_full2 <- dplyr::group_by(df_temp_full, type, type2)
    df_temp_full2 <- dplyr::summarize(df_temp_full2, sum = sum(r), n = n())
    df_temp_full2 <- dplyr::group_by(df_temp_full2, type2)
    df_temp_full2 <- dplyr::arrange(df_temp_full2, desc(n), desc(sum))
    df_temp_full2 <- dplyr::filter(df_temp_full2, type != paste0("r<", threshold,", unassigned"))
    df_temp_full2 <- dplyr::slice(df_temp_full2, 1)
    df_temp_full2 <- dplyr::right_join(df_temp_full2, select(df_temp_full, -type), by = stats::setNames(collapse_to_cluster, "type2"))
    df_temp_full <- dplyr::mutate(df_temp_full2, type = tidyr::replace_na(type, paste0("r<", threshold,", unassigned")))
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
  df_temp[["type"]][df_temp$r < threshold] <- paste0("r<", threshold,", unassigned")
  df_temp <- dplyr::top_n(dplyr::group_by_at(df_temp, 1), topn, r)
  df_temp_full <- df_temp

  if(collapse_to_cluster != FALSE){
    if (!(col %in% colnames(metadata))) {
      metadata <- tibble::as_tibble(metadata, rownames = col)
    }
    df_temp_full <- dplyr::left_join(df_temp_full, metadata, by = col)
    df_temp_full[,"type2"] <- df_temp_full[[collapse_to_cluster]]
    df_temp_full2 <- dplyr::group_by(df_temp_full, type, type2)
    df_temp_full2 <- dplyr::summarize(df_temp_full2, sum = sum(r), n = n())
    df_temp_full2 <- dplyr::group_by(df_temp_full2, type2)
    df_temp_full2 <- dplyr::arrange(df_temp_full2, desc(n), desc(sum))
    df_temp_full2 <- dplyr::filter(df_temp_full2, type != paste0("r<", threshold,", unassigned"))
    df_temp_full2 <- dplyr::slice(df_temp_full2, 1:topn)
    df_temp_full2 <- dplyr::right_join(df_temp_full2, select(df_temp_full, -c(type, r)), by = stats::setNames(collapse_to_cluster, "type2"))
    df_temp_full <- dplyr::mutate(df_temp_full2, type = tidyr::replace_na(type, paste0("r<", threshold,", unassigned")))
    df_temp_full <- dplyr::group_by_(df_temp_full, .dots = col)
    df_temp_full <- dplyr::distinct(df_temp_full, type, type2, .keep_all = T)
    dplyr::arrange(df_temp_full, desc(n), desc(sum), .by_group = T)
  } else {
    df_temp_full <- dplyr::group_by_(df_temp_full, .dots = col)
    dplyr::arrange(df_temp_full, desc(r), .by_group = T)
  }

}

#' pct of cells in each cluster that express genelist
#'
#' @param matrix expression matrix
#' @param genelist vector of marker genes for one identity
#' @param clusters vector of cluster identities
#'
#' @export
gene_pct <- function(matrix, genelist, clusters){
  genelist <- intersect(genelist, rownames(matrix))
  unique_clusters <- unique(clusters)
  sapply(unique_clusters, function(x) {
    celllist <- clusters == x
    tmp <- matrix[genelist, celllist, drop = F]
    tmp[tmp > 0] <- 1
    mean(Matrix::rowSums(tmp)/ncol(tmp))
  })
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
gene_pct_markerm <- function(matrix, marker_m, cluster_info, cluster_col = NULL, norm = NULL) {
  if(is.vector(cluster_info)){
  } else if (is.data.frame(cluster_info) & !is.null(cluster_col)){
    cluster_info <- cluster_info[[cluster_col]]
  } else {
    stop("cluster_info not formatted correctly,
         supply either a  vector or a dataframe")
  }

  out <- sapply(colnames(marker_m), function(x) {
    gene_pct(matrix,
             marker_m[[x]],
             cluster_info)
  })

  if (!(is.null(norm))) {
    if (norm == "divide") {
      out <- sweep(out,2,apply(out,2,max),"/")
    } else if (norm == "diff") {
      out <- sweep(out,2,apply(out,2,max),"-")
    } else {
      out <- sweep(out,2,apply(out,2,max) * norm)
      out[out < 0] <- 0
      out[out > 0] <- 1
    }
  }
  out
}

#' Combined function to compare scRNA-seq data to bulk RNA-seq data and marker list
#'
#'@export
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
#' @param dr stored dimension reduction
#' @param seurat_out output cor matrix or called seurat object
#' @param ... passed to matrixize_markers
#' @param norm whether and how the results are normalized
#' @param marker_inmatrix whether markers genes are already in preprocessed matrix form

#' @export
clustify_nudge.seurat <- function(input,
                           ref_mat,
                           marker,
                           cluster_col = NULL,
                           query_genes = NULL,
                           compute_method = "spearman",
                           weight = 1,
                           seurat_out = T,
                           threshold = -Inf,
                           dr = "tsne",
                           norm = "diff",
                           marker_inmatrix = T,
                           ...){
  if (marker_inmatrix != T) {
    marker <- matrixize_markers(marker,
                                ...)
  }

  resb <- gene_pct_markerm(input@data, marker,
                           input@meta.data,
                           cluster_col = cluster_col,
                           norm = norm)

  resa <- clustify(
    input = input,
    ref_mat = ref_mat,
    cluster_col = cluster_col,
    query_genes = query_genes,
    seurat_out = F,
    per_cell = F
  )

  df_temp <- cor_to_call(resa[order(rownames(resa)), order(colnames(resa))] +
                           resb[order(rownames(resb)), order(colnames(resb))] * weight,
                         threshold = threshold)
  colnames(df_temp) <- c(cluster_col, "type", "score")

  if (seurat_out == F) {
    df_temp
  } else {
    cluster_info <- as.data.frame(use_seurat_meta(input, dr = dr, seurat3 = F))
    df_temp_full <- dplyr::left_join(tibble::rownames_to_column(cluster_info, "rn"), df_temp, by = cluster_col)
    df_temp_full <- tibble::column_to_rownames(df_temp_full, "rn")
    if ("Seurat" %in% loadedNamespaces()) {
      input@meta.data <- df_temp_full
      return(input)
    } else {
      print("seurat not loaded, returning cor_mat instead")
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

#' @export
clustify_nudge.default <- function(input,
                                  ref_mat,
                                  marker,
                                  metadata = NULL,
                                  cluster_col = NULL,
                                  query_genes = NULL,
                                  compute_method = "spearman",
                                  weight = 1,
                                  seurat_out = T,
                                  threshold = -Inf,
                                  dr = "tsne",
                                  norm = "diff",
                                  call = T,
                                  marker_inmatrix = T,
                                  ...) {
  if (marker_inmatrix != T) {
    marker <- matrixize_markers(marker,
                                ...)
  }

  if (!(stringr::str_detect(class(input), "atrix"))) {
    input_original <- input
    temp <- parse_loc_object(input,
                             type = class(input),
                             expr_loc = NULL,
                             meta_loc = NULL,
                             var_loc = NULL,
                             cluster_col = cluster_col)
    input <- temp[["expr"]]
    metadata <- temp[["meta"]]
    if (is.null(query_genes)) {
    query_genes <- temp[["var"]]
    }
    if (is.null(cluster_col)) {
      cluster_col <- temp[["col"]]
    }
  }

  resb <- gene_pct_markerm(input, marker,
                           metadata,
                           cluster_col = cluster_col,
                           norm = norm)

  resa <- clustify(
    input = input,
    ref_mat = ref_mat,
    metadata = metadata,
    cluster_col = cluster_col,
    query_genes = query_genes,
    seurat_out = F,
    per_cell = F
  )

  if (call == T) {
    df_temp <- cor_to_call(resa[order(rownames(resa)), order(colnames(resa))] +
                             resb[order(rownames(resb)), order(colnames(resb))] * weight,
                           threshold = threshold)
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
#' @export
parse_loc_object <- function(input,
                             type = class(input),
                             expr_loc = NULL,
                             meta_loc = NULL,
                             var_loc = NULL,
                             cluster_col = NULL) {
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
  if (type %in% colnames(object_loc_lookup)) {
    parsed = list(eval(parse(text=object_loc_lookup[[type]][1])),
                  as.data.frame(eval(parse(text=object_loc_lookup[[type]][2]))),
                  eval(parse(text=object_loc_lookup[[type]][3])),
                  object_loc_lookup[[type]][4])
  } else {
    parsed = list(NULL, NULL, NULL, NULL)
  }

  names(parsed) <- c("expr", "meta", "var", "col")

  if (!(is.null(expr_loc))) {
    parsed[["expr"]] <- eval(parse(text = paste0("input",expr_loc)))
  }

  if (!(is.null(meta_loc))) {
   parsed[["meta"]] <- as.data.frame(eval(parse(text = paste0("input",meta_loc))))
  }

  if (!(is.null(var_loc))) {
    parsed[["var"]] <- eval(parse(text = paste0("input",var_loc)))
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
#' @param cluster_col column of clustering from metadata
#' @param x_col column of metadata for x axis plotting
#' @param y_col column of metadata for y axis plotting
#' @param n expand n-fold for over/under clustering
#' @param do.label whether to label each cluster at median center
#' @param seed set seed for kmeans
#' @param newclustering use kmeans if NULL on dr or col name for second column of clustering

#' @export
overcluster_test <- function(expr,
                             metadata,
                             cluster_col,
                             x_col = "tSNE_1",
                             y_col = "tSNE_2",
                             n = 5,
                             do.label = T,
                             seed = 42,
                             newclustering = NULL) {
  if (!(is.null(seed))) {
    set.seed(seed)
  }

  if (is.null(newclustering)) {
    metadata$new_clusters <- as.character(kmeans(metadata[,c(x_col, y_col)],
                                                 centers = n * length(unique(metadata[[cluster_col]])))$clust)
  } else {
    metadata$new_clusters <- metadata[[newclustering]]
    n <- length(unique(metadata[[newclustering]]))/length(unique(metadata[[cluster_col]]))
  }
  res1 <- clustify_nudge(expr,
                         cbmc_ref,
                         cbmc_m,
                         metadata,
                         query_genes = pbmc4k_markers_M3Drop$Gene,
                         cluster_col = cluster_col,
                         call = F)
  res2 <- clustify_nudge(expr,
                         cbmc_ref,
                         cbmc_m,
                         metadata,
                         query_genes = pbmc4k_markers_M3Drop$Gene,
                         cluster_col = "new_clusters",
                         call = F)
  o1 <- plot_tsne(metadata,
                  feature = cluster_col,
                  x = x_col,
                  y = y_col,
                  do.label = F,
                  do.legend = F)
  o2 <- plot_tsne(metadata,
                  feature = "new_clusters",
                  x = x_col,
                  y = y_col,
                  do.label = F,
                  do.legend = F)
  p1 <- plot_best_call(res1,
                       metadata,
                       cluster_col,
                       do.label = do.label,
                       x = x_col,
                       y = y_col)
  p2 <- plot_best_call(res2,
                       metadata,
                       "new_clusters",
                       do.label = do.label,
                       x = x_col,
                       y = y_col)
  cowplot::plot_grid(o1, o2, p1, p2,
                     labels = c(length(unique(metadata[[cluster_col]])),
                                n * length(unique(metadata[[cluster_col]]))))
}
