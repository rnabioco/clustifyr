#' Binarize scRNAseq data
#'
#' @param mat single-cell expression matrix
#' @param n number of top expressing genes to keep
#' @param cut cut off to set to 0
#' @return matrix of 1s and 0s
#' @examples
#' pbmc_avg <- average_clusters(
#'   mat = pbmc_matrix_small,
#'   metadata = pbmc_meta,
#'   cluster_col = "classified"
#' )
#'
#' pbmc_avgb <- binarize_expr(pbmc_avg)
#' @export
binarize_expr <- function(mat,
                          n = 1000,
                          cut = 0) {
  expr_mat <- mat
  if (n < nrow(expr_mat)) {
    expr_df <- as.data.frame(as.matrix(expr_mat))
    df_temp <-
      apply(expr_df, 2, function(x) {
        x - x[order(x, decreasing = TRUE)[n + 1]]
      })
    df_temp[df_temp > cut] <- 1
    df_temp[df_temp < cut] <- 0
    as.matrix(df_temp)
  } else {
    expr_mat[expr_mat > cut] <- 1
    expr_mat[expr_mat < cut] <- 0
    expr_mat
  }
}

#' Convert candidate genes list into matrix
#'
#' @param marker_df dataframe of candidate genes, must contain "gene" and "cluster" columns, or a matrix of gene names to convert to ranked
#' @param ranked unranked gene list feeds into hyperp, the ranked
#' gene list feeds into regular corr_coef
#' @param n number of genes to use
#' @param step_weight ranked genes are tranformed into pseudo expression by descending weight
#' @param background_weight ranked genes are tranformed into pseudo expression with
#' added weight
#' @param unique whether to use only unique markers to 1 cluster
#' @param metadata vector or dataframe of cluster names, should have column named cluster
#' @param cluster_col column for cluster names to replace original cluster, if metadata is dataframe
#' @param remove_rp do not include rps, rpl, rp1-9 in markers
#' @return matrix of unranked gene marker names, or matrix of ranked expression
#' @examples
#' pbmc_mm <- matrixize_markers(pbmc_markers)
#' @export
matrixize_markers <- function(marker_df,
                              ranked = FALSE,
                              n = NULL,
                              step_weight = 1,
                              background_weight = 0,
                              unique = FALSE,
                              metadata = NULL,
                              cluster_col = "classified",
                              remove_rp = FALSE) {
  # takes marker in dataframe form
  # equal number of marker genes per known cluster
  marker_df <- dplyr::as_tibble(marker_df)

  # if "gene" not present in column names, assume df is a matrix to be converted to ranked
  if (!("gene" %in% colnames(marker_df))) {
    marker_df <-
      data.frame(lapply(marker_df, as.character), stringsAsFactors = FALSE)
    marker_df <-
      tidyr::gather(marker_df,
        factor_key = TRUE,
        key = "cluster",
        value = "gene"
      )
  }

  if (remove_rp) {
    marker_df <-
      dplyr::filter(marker_df, !(stringr::str_detect(gene, "^RP[0-9,L,S]|^Rp[0-9,l,s]")))
  }

  if (unique) {
    nonunique <- dplyr::group_by(marker_df, gene)
    nonunique <- dplyr::summarize(nonunique, n = dplyr::n())
    nonunique <- dplyr::filter(nonunique, n > 1)
    marker_df <- dplyr::anti_join(marker_df, nonunique, by = "gene")
  }

  cut_temp <- dplyr::group_by(marker_df, cluster)
  cut_temp <- dplyr::summarize(cut_temp, n = n())
  cut_num <- min(cut_temp$n)

  if (!is.null(n)) {
    if (n < cut_num) {
      cut_num <- n
    }
  }

  marker_temp <- dplyr::select(marker_df, gene, cluster)
  marker_temp <- dplyr::group_by(marker_temp, cluster)
  marker_temp <- dplyr::slice(marker_temp, 1:cut_num)
  if (ranked) {
    marker_temp <-
      dplyr::mutate(
        marker_temp,
        n = seq(
          step_weight * cut_num,
          by = -step_weight,
          length.out = cut_num
        ) + background_weight
      )
    marker_temp2 <-
      tidyr::spread(marker_temp, key = "cluster", value = n)
    marker_temp2 <-
      as.data.frame(replace(marker_temp2, is.na(marker_temp2), 0))
    rownames(marker_temp2) <- marker_temp2$gene
    marker_temp2 <- dplyr::select(marker_temp2, -gene)
  } else {
    marker_temp <- dplyr::mutate(marker_temp, n = 1:cut_num)
    marker_temp2 <-
      tidyr::spread(marker_temp, key = "cluster", value = "gene")
    marker_temp2 <- as.data.frame(dplyr::select(marker_temp2, -n))
  }

  # if metadata is vector, adopt names in vector; if metadata is a metadata dataframe, pulls names from cluster_col column
  if (!is.null(metadata)) {
    if (typeof(metadata) != "character") {
      metadata <-
        suppressWarnings(dplyr::left_join(
          tibble::tibble(cluster = colnames(marker_temp2)),
          unique(
            tibble::tibble(
              cluster = metadata$cluster,
              classified = metadata[[cluster_col]]
            )
          ),
          by = "cluster"
        ))
      metadata <- metadata[[cluster_col]]
    }
    colnames(marker_temp2) <- metadata
  }

  marker_temp2
}

#' Generate variable gene list from marker matrix
#'
#' @description Variable gene list is required for `clustify` main function. This
#' function parses variables genes from a matrix input.
#'
#' @param marker_mat matrix or dataframe of candidate genes for each cluster
#' @return vector of marker gene names
#' @examples
#' get_vargenes(cbmc_m)
#' @export
get_vargenes <- function(marker_mat) {
  if (rownames(marker_mat)[1] != "1") {
    unique(rownames(marker_mat))
  } else {
    unique(unlist(marker_mat, use.names = FALSE))
  }
}

#' Calculate adjusted p-values for hypergeometric test of gene lists
#' or jaccard index
#'
#' @param bin_mat binarized single-cell expression matrix, feed in by_cluster mat, if desired
#' @param marker_mat matrix or dataframe of candidate genes for each cluster
#' @param n number of genes in the genome
#' @param metric adjusted p-value for hypergeometric test, or jaccard index
#' @param output_high if true (by default to fit with rest of package),
#' -log10 transform p-value
#' @return matrix of numeric values, clusters from expr_mat as row names, cell types from marker_mat as column names
#' @examples
#' pbmc_mm <- matrixize_markers(pbmc_markers)
#'
#' pbmc_avg <- average_clusters(
#'   pbmc_matrix_small,
#'   pbmc_meta,
#'   cluster_col = "classified"
#' )
#'
#' pbmc_avgb <- binarize_expr(pbmc_avg)
#'
#' compare_lists(
#'   pbmc_avgb,
#'   pbmc_mm,
#'   metric = "spearman"
#' )
#' @export
compare_lists <- function(bin_mat,
                          marker_mat,
                          n = 30000,
                          metric = "hyper",
                          output_high = TRUE) {
  # check if matrix is binarized
  if ((length(unique(bin_mat[, 1])) > 2) & (metric != "gsea")) {
    warning("non-binarized data, running spearman instead")
    metric <- "spearman"
  }

  if (metric == "hyper") {
    out <- lapply(
      colnames(bin_mat),
      function(x) {
        per_col <- lapply(
          colnames(marker_mat),
          function(y) {
            marker_list <- unlist(marker_mat[, y], use.names = FALSE)
            bin_temp <- bin_mat[, x][bin_mat[, x] == 1]
            list_top <- names(bin_temp)

            t <- length(intersect(list_top, marker_list))
            a <- max(length(list_top), length(marker_list))
            b <- min(length(list_top), length(marker_list))
            sum(dhyper(t:b, a, n - a, b))
          }
        )
        do.call(cbind, as.list(p.adjust(per_col)))
      }
    )
    if (any(sapply(out, is.na))) {
      error("NaN produced, possibly due to wrong n")
    }
  } else if (metric == "jaccard") {
    out <- lapply(
      colnames(bin_mat),
      function(x) {
        per_col <- lapply(
          colnames(marker_mat),
          function(y) {
            marker_list <- unlist(marker_mat[, y], use.names = FALSE)
            bin_temp <- bin_mat[, x][bin_mat[, x] == 1]
            list_top <- names(bin_temp)

            I <- length(intersect(list_top, marker_list))
            I / (length(list_top) + length(marker_list) - I)
          }
        )
        do.call(cbind, per_col)
      }
    )
  } else if (metric == "spearman") {
    out <- lapply(
      colnames(bin_mat),
      function(x) {
        per_col <- lapply(
          colnames(marker_mat),
          function(y) {
            marker_list <- unlist(marker_mat[, y], use.names = FALSE)
            v1 <-
              marker_list[marker_list %in% names(as.matrix(bin_mat)[, x])]
            bin_temp <- as.matrix(bin_mat)[, x]
            bin_temp <- bin_temp[order(bin_temp, decreasing = TRUE)]
            list_top <- names(bin_temp)
            v2 <- list_top[list_top %in% v1]
            v1 <- v1
            v2 <- v2
            sum(sapply(seq_along(v1), function(i) {
              abs(i - (
                which(v2 == v1[i])
              ))
            }))
          }
        )
        do.call(cbind, per_col)
      }
    )
  } else if (metric == "gsea") {
    out <- lapply(
      colnames(marker_mat),
      function(y) {
        marker_list <- list()
        marker_list[[1]] <- marker_mat[, y]
        names(marker_list) <- y
        v1 <- marker_list
        run_gsea(bin_mat, v1, n_perm = 1000, per_cell = TRUE)
      }
    )
    res <- do.call(cbind, out)
    n <- ncol(res)
    res2 <- res[, 3 * c(1:(ncol(res) / 3)) - 1, drop = FALSE]
    rownames(res2) <- rownames(res)
    colnames(res2) <- colnames(marker_mat)
    res <- res2
    if (output_high) {
      res <- -log10(res)
    }
  } else {
    stop("unrecognized metric", call. = FALSE)
  }

  if (metric != "gsea") {
    res <- do.call(rbind, out)
    rownames(res) <- colnames(bin_mat)
    colnames(res) <- colnames(marker_mat)
  }

  if (output_high) {
    if (metric == "hyper") {
      res <- -log10(res)
    } else if (metric == "spearman") {
      res <- -res + max(res)
    }
  }

  res
}
