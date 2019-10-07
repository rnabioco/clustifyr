#' get best calls for each cluster
#'
#' @param correlation_matrix input similarity matrix
#' @param metadata input metadata with tsne or umap coordinates and cluster ids
#' @param cluster_col metadata column, can be cluster or cellid
#' @param collapse_to_cluster if a column name is provided, takes the most frequent call of entire cluster to color in plot
#' @param threshold minimum correlation coefficent cutoff for calling clusters
#' @param rename_prefix prefix to add to type and r column names
#' @return dataframe of cluster, new ident, and r info
#' @examples
#' res <- clustify(
#'   input = pbmc_matrix_small,
#'   metadata = pbmc_meta,
#'   cluster_info = pbmc_meta,
#'   cluster_col = "classified",
#'   ref_mat = pbmc_bulk_matrix
#' )
#'
#' res2 <- cor_to_call(res)
#' @export
cor_to_call <- function(correlation_matrix,
                        metadata = NULL,
                        cluster_col = "cluster",
                        collapse_to_cluster = FALSE,
                        threshold = 0,
                        rename_prefix = NULL) {
  if (threshold == "auto") {
    threshold = round(0.75 * max(correlation_matrix), 2)
    message(paste0("using threshold of ", threshold))
  }
  df_temp <- tibble::as_tibble(correlation_matrix, rownames = cluster_col)
  df_temp <- tidyr::gather(df_temp, key = !!dplyr::sym("type"), value = !!dplyr::sym("r"), -!!cluster_col)
  df_temp[["type"]][df_temp$r < threshold] <- paste0("r<", threshold, ", unassigned")
  df_temp <- dplyr::top_n(dplyr::group_by_at(df_temp, 1), 1, !!dplyr::sym("r"))
  if (nrow(df_temp) != nrow(correlation_matrix)) {
    clash <- dplyr::summarize(dplyr::group_by_at(df_temp, 1), n = n())
    clash <- dplyr::filter(clash, n > 1)
    clash <- dplyr::pull(clash, 1)
    df_temp[lapply(df_temp[, 1], FUN = function(x) x %in% clash)[[1]], 2] <- paste0(df_temp[["type"]][lapply(df_temp[, 1], FUN = function(x) x %in% clash)[[1]]], "-CLASH!")
    df_temp <- dplyr::distinct(df_temp, exclude = "type", .keep_all = TRUE)
    df_temp_full <- dplyr::select(df_temp, -!!dplyr::sym("exclude"))
  } else {
    df_temp_full <- df_temp
  }

  if (collapse_to_cluster != FALSE) {
    if (!(cluster_col %in% colnames(metadata))) {
      metadata <- tibble::as_tibble(metadata, rownames = cluster_col)
    }
    df_temp_full <- dplyr::left_join(df_temp_full, metadata, by = cluster_col)
    df_temp_full[, "type2"] <- df_temp_full[[collapse_to_cluster]]
    df_temp_full2 <- dplyr::group_by(df_temp_full, !!dplyr::sym("type"), !!dplyr::sym("type2"))
    df_temp_full2 <- dplyr::summarize(df_temp_full2, sum = sum(!!dplyr::sym("r")), n = n())
    df_temp_full2 <- dplyr::group_by(df_temp_full2, !!dplyr::sym("type2"))
    df_temp_full2 <- dplyr::arrange(df_temp_full2, desc(n), desc(sum))
    df_temp_full2 <- dplyr::filter(df_temp_full2, !!dplyr::sym("type") != paste0("r<", threshold, ", unassigned"))
    df_temp_full2 <- dplyr::slice(df_temp_full2, 1)
    df_temp_full2 <- dplyr::right_join(df_temp_full2, dplyr::select(df_temp_full, -!!dplyr::sym("type")), by = stats::setNames(collapse_to_cluster, "type2"))
    df_temp_full <- dplyr::mutate(df_temp_full2, type = tidyr::replace_na(!!dplyr::sym("type"), paste0("r<", threshold, ", unassigned")))
    eval(parse(text = paste0("df_temp_full <- dplyr::rename(df_temp_full,", collapse_to_cluster, " = type2.y)")))
    df_temp_full <- dplyr::select(dplyr::ungroup(df_temp_full), -!!dplyr::sym("type2"))
  }

  if (!is.null(rename_prefix)) {
    eval(parse(text = paste0("df_temp_full <- dplyr::rename(df_temp_full, ", paste0(rename_prefix, "_type"), " = type, ", paste0(rename_prefix, "_r"), " = r)")))
  }
  df_temp_full
}

#' Insert called ident results into metadata
#'
#' @param res dataframe of idents, such as output of cor_to_call
#' @param metadata input metadata with tsne or umap coordinates and cluster ids
#' @param cluster_col metadata column, can be cluster or cellid
#' @param per_cell whether the res dataframe is listed per cell
#' @param rename_prefix prefix to add to type and r column names
#' @return new metadata with added columns
#' @examples
#' \donttest{
#' res <- clustify(
#'   input = pbmc_matrix_small,
#'   metadata = pbmc_meta,
#'   cluster_info = pbmc_meta,
#'   cluster_col = "classified",
#'   ref_mat = pbmc_bulk_matrix
#' )
#'
#' res2 <- cor_to_call(res)
#'
#' pbmc
#'
#' pbmc_meta2 <- call_to_metadata(
#'   res = res2,
#'   metadata = pbmc_meta,
#'   cluster_col = "cluster",
#'   rename_suff = "_assigned"
#' )
#' }
#' @export
call_to_metadata <- function(res,
                             metadata,
                             cluster_col,
                             per_cell = FALSE,
                             rename_prefix = NULL) {
  temp_col_id <- get_unique_column(metadata, "rn")

  df_temp <- res
  if (per_cell == FALSE) {
    if (!(cluster_col %in% colnames(metadata))) {
      stop("cluster_col is not a column of metadata", call. = FALSE)
    }

    if (!(cluster_col %in% colnames(res))) {
      stop("cluster_col is not a column of called cell type dataframe", call. = FALSE)
    }

    if (!(all(unique(df_temp[[cluster_col]]) %in% unique(metadata[[cluster_col]])))) {
      stop("cluster_col from clustify step and joining to metadata step are not the same", call. = FALSE)
    }

    df_temp_full <- suppressWarnings(dplyr::left_join(tibble::rownames_to_column(metadata,
                                                                                 temp_col_id),
                                                      df_temp,
                                                      by = cluster_col,
                                                      suffix = c("", ".clustify")))

    df_temp_full <- tibble::column_to_rownames(df_temp_full,
                                               temp_col_id)
  } else {
    colnames(df_temp)[1] <- cluster_col
    names(cluster_col) <- temp_col_id

    df_temp_full <- suppressWarnings(dplyr::left_join(tibble::rownames_to_column(metadata,
                                                                                 temp_col_id),
                                                      df_temp,
                                                      by = cluster_col,
                                                      suffix = c("", ".clustify")))

    df_temp_full <- tibble::column_to_rownames(df_temp_full, temp_col_id)
  }
  if (!is.null(rename_prefix)) {
    eval(parse(text = paste0("df_temp_full <- dplyr::rename(df_temp_full, ", paste0(rename_prefix, "_type"), " = type, ", paste0(rename_prefix, "_r"), " = r)")))
  }
  df_temp_full
}

#' From per-cell calls, take highest freq call in each cluster
#'
#' @param res dataframe of idents, such as output of cor_to_call
#' @param metadata input metadata with tsne or umap coordinates and cluster ids
#' @param cluster_col metadata column for cluster
#' @param threshold minimum correlation coefficent cutoff for calling clusters
#' @return new metadata with added columns
#' @export
collapse_to_cluster <- function(res,
                                metadata,
                                cluster_col,
                                threshold = 0) {
  df_temp_full <- res
  df_temp_full2 <- dplyr::mutate(df_temp_full, type2 = metadata[[cluster_col]])
  df_temp_full2 <- dplyr::group_by(df_temp_full2, !!dplyr::sym("type"), !!dplyr::sym("type2"))
  df_temp_full2 <- dplyr::summarize(df_temp_full2, sum = sum(!!dplyr::sym("r")), n = n())
  df_temp_full2 <- dplyr::group_by(df_temp_full2, !!dplyr::sym("type2"))
  df_temp_full2 <- dplyr::arrange(df_temp_full2, desc(n), desc(sum))
  df_temp_full2 <- dplyr::filter(df_temp_full2, !!dplyr::sym("type") != paste0("r<", threshold, ", unassigned"))
  df_temp_full2 <- dplyr::slice(df_temp_full2, 1)
  df_temp_full2 <- dplyr::right_join(df_temp_full2, select(df_temp_full, -!!dplyr::sym("type")), by = stats::setNames(cluster_col, "type2"))
  df_temp_full <- dplyr::mutate(df_temp_full2, type = tidyr::replace_na(!!dplyr::sym("type"), paste0("r<", threshold, ", unassigned")))
}
