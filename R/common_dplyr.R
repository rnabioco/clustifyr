#' get best calls for each cluster
#'
#' @param correlation_matrix input similarity matrix
#' @param metadata input metadata with tsne or umap coordinates and cluster ids
#' @param cluster_col metadata column, can be cluster or cellid
#' @param collapse_to_cluster if a column name is provided, takes the most frequent call of entire cluster to color in plot
#' @param threshold minimum correlation coefficent cutoff for calling clusters
#' @param rename_suff suffix to add to type and r column names
#' @return dataframe of cluster, new ident, and r info
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

#' Insert called ident results into metadata
#'
#' @param res dataframe of idents, such as output of cor_to_call
#' @param metadata input metadata with tsne or umap coordinates and cluster ids
#' @param cluster_col metadata column, can be cluster or cellid
#' @param rename_suff suffix to add to type and r column names
#' @return new metadata with added columns
#' @export
call_to_metadata <- function(res,
                             metadata,
                             cluster_col,
                             per_cell = FALSE,
                             rename_suff = NULL) {
  df_temp <- res
  if (per_cell == FALSE) {
    df_temp_full <- dplyr::left_join(tibble::rownames_to_column(metadata, "rn"), df_temp, by = cluster_col)
    df_temp_full <- tibble::column_to_rownames(df_temp_full, "rn")
  } else {
    df_temp <- dplyr::rename(df_temp, "rn" = 1)
    df_temp_full <- dplyr::left_join(tibble::rownames_to_column(metadata, "rn"), df_temp, by = "rn")
    df_temp_full <- tibble::column_to_rownames(df_temp_full, "rn")
  }
  if (!is.null(rename_suff)) {
    eval(parse(text = paste0("df_temp_full <- dplyr::rename(df_temp_full, ", paste0(rename_suff, "_type"), " = type, ", paste0(rename_suff, "_r"), " = r)")))
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
  df_temp_full2 <- dplyr::group_by(df_temp_full2, type, type2)
  df_temp_full2 <- dplyr::summarize(df_temp_full2, sum = sum(r), n = n())
  df_temp_full2 <- dplyr::group_by(df_temp_full2, type2)
  df_temp_full2 <- dplyr::arrange(df_temp_full2, desc(n), desc(sum))
  df_temp_full2 <- dplyr::filter(df_temp_full2, type != paste0("r<", threshold, ", unassigned"))
  df_temp_full2 <- dplyr::slice(df_temp_full2, 1)
  df_temp_full2 <- dplyr::right_join(df_temp_full2, select(df_temp_full, -type), by = stats::setNames(cluster_col, "type2"))
  df_temp_full <- dplyr::mutate(df_temp_full2, type = tidyr::replace_na(type, paste0("r<", threshold, ", unassigned")))
}
