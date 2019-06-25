#' Binarize scRNA seq data
#'
#' @param expr_mat single-cell expression matrix
#' @param n number of top expressing genes to keep
#' @param cut cut off to set to 0
#'
#' @export
binarize_expr <- function(expr_mat,
                          n = 1000,
                          cut = 0) {
  if (n < nrow(expr_mat)) {
    expr_df <- as.data.frame(as.matrix(expr_mat))
    df_temp <- apply(expr_df, 2, function(x) x - x[order(x, decreasing = TRUE)[n + 1]])
    df_temp[df_temp > cut] <- 1
    df_temp[df_temp < cut] <- 0
    as.matrix(df_temp)
  } else {
    expr_mat[expr_mat > cut] <- 1
    expr_mat[expr_mat < cut] <- 0
    expr_mat
  }
}

#' convert candidate genes list into matrix
#'
#' @param marker_df dataframe of candidate genes, must contain "gene" and "cluster" columns, or a matrix of gene names to convert to ranked
#' @param ranked unranked gene list feeds into hyperp, the ranked
#' gene list feeds into regular corr_coef
#' @param n number of genes to use
#' @param step_weight ranked genes are tranformed into pseudo expression by descending weight
#' @param background_weight ranked genes are tranformed into pseudo expression with
#' added weight
#' @param unique whether to use only unique markers to 1 cluster
#' @param labels vector or dataframe of cluster names
#' @param remove_rp do not include rps, rpl, rp1-9 in markers
#'
#' @export
matrixize_markers <- function(marker_df,
                              ranked = FALSE,
                              n = NULL,
                              step_weight = 1,
                              background_weight = 0,
                              unique = FALSE,
                              labels = NULL,
                              remove_rp = FALSE) {
  # takes marker in dataframe form
  # equal number of marker genes per known cluster
  marker_df <- dplyr::as_tibble(marker_df)

  # if "gene" not present in column names, assume df is a matrix to be converted to ranked
  if (!("gene" %in% colnames(marker_df))) {
    marker_df <- data.frame(lapply(marker_df, as.character), stringsAsFactors = FALSE)
    marker_df <- tidyr::gather(marker_df, factor_key = TRUE, key = "cluster", value = "gene")
  }

  if (remove_rp == TRUE) {
    marker_df <- dplyr::filter(marker_df, !(stringr::str_detect(gene, "^RP[0-9,L,S]")))
  }

  if (unique == TRUE) {
    nonunique <- dplyr::group_by(marker_df, gene)
    nonunique <- dplyr::summarize(nonunique, n = n())
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
  if (ranked == TRUE) {
    marker_temp <- dplyr::mutate(marker_temp, n = seq(step_weight * cut_num, by = -step_weight, length.out = cut_num) + background_weight)
    marker_temp2 <- tidyr::spread(marker_temp, key = "cluster", value = n)
    marker_temp2 <- as.data.frame(replace(marker_temp2, is.na(marker_temp2), 0))
    rownames(marker_temp2) <- marker_temp2$gene
    marker_temp2 <- dplyr::select(marker_temp2, -gene)
  } else {
    marker_temp <- dplyr::mutate(marker_temp, n = 1:cut_num)
    marker_temp2 <- tidyr::spread(marker_temp, key = "cluster", value = "gene")
    marker_temp2 <- as.data.frame(dplyr::select(marker_temp2, -n))
  }

  # if labels is vector, adopt names in vector; if labels is a metadata dataframe, pulls names from "classified" column
  if (!is.null(labels)) {
    if (typeof(labels) != "character") {
      label_df <- labels
      labels <- left_join(data_frame(cluster = colnames(marker_temp2)),
        unique(data_frame(
          cluster = labels$cluster,
          classified = labels$classified
        )),
        by = "cluster"
      )
      labels <- dplyr::pull(labels, classified)
    }
    colnames(marker_temp2) <- labels
  }

  marker_temp2
}

#' generate variable gene list from marker matrix
#'
#' @description Variable gene list is required for `clustify` main function. This
#' function parses variables genes from a matrix input.
#'
#' @param marker_m matrix or dataframe of candidate genes for each cluster
#'
#' @export
get_vargenes <- function(marker_m) {
  if (rownames(marker_m)[1] != "1") {
    unique(rownames(marker_m))
  } else {
    unique(unlist(marker_m, use.names = FALSE))
  }
}

#' calculate adjusted p-values for hypergeometric test of gene lists
#' or jaccard index
#'
#' @param bin_mat binarized single-cell expression matrix, feed in by_cluster mat, if desired
#' @param marker_m matrix or dataframe of candidate genes for each cluster
#' @param n number of genes in the genome
#' @param metric adjusted p-value for hypergeometric test, or jaccard index
#' @param output_high if true (by default to fit with rest of package),
#' -log10 transform p-value
#'
#' @export
compare_lists <- function(bin_mat,
                          marker_m,
                          n = 30000,
                          metric = "hyper",
                          output_high = TRUE) {
  # check if matrix is binarized
  if ((length(unique(bin_mat[, 1])) > 2) & (metric != "gsea")) {
    warning("non-binarized data, running spearman instead")
    metric <- "spearman"
  }

  # "expressed" genes per single cell data cluster
  if (metric == "hyper") {
    out <- lapply(
      colnames(bin_mat),
      function(x) {
        per_col <- lapply(
          colnames(marker_m),
          function(y) {
            marker_list <- unlist(marker_m[, y], use.names = FALSE)
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
  }

  if (metric == "jaccard") {
    out <- lapply(
      colnames(bin_mat),
      function(x) {
        per_col <- lapply(
          colnames(marker_m),
          function(y) {
            marker_list <- unlist(marker_m[, y], use.names = FALSE)
            bin_temp <- bin_mat[, x][bin_mat[, x] == 1]
            list_top <- names(bin_temp)

            I <- length(intersect(list_top, marker_list))
            I / (length(list_top) + length(marker_list) - I)
          }
        )
        do.call(cbind, per_col)
      }
    )
  }

  if (metric == "spearman") {
    out <- lapply(
      colnames(bin_mat),
      function(x) {
        per_col <- lapply(
          colnames(marker_m),
          function(y) {
            marker_list <- unlist(marker_m[, y], use.names = FALSE)
            v1 <- marker_list[marker_list %in% names(as.matrix(bin_mat)[, x])]
            bin_temp <- as.matrix(bin_mat)[, x]
            bin_temp <- bin_temp[order(bin_temp, decreasing = TRUE)]
            list_top <- names(bin_temp)
            v2 <- list_top[list_top %in% v1]
            v1 <<- v1
            v2 <<- v2
            sum(sapply(seq_along(v1), function(i) abs(i - (which(v2 == v1[i])))))
          }
        )
        do.call(cbind, per_col)
      }
    )
  }

  if (metric != "gsea") {
    res <- do.call(rbind, out)
    rownames(res) <- colnames(bin_mat)
    colnames(res) <- colnames(marker_m)
  }

  if (metric == "gsea") {
    out <- lapply(
      colnames(marker_m),
      function(y) {
        marker_list <- list()
        marker_list[[1]] <- marker_m[, y]
        names(marker_list) <- y
        v1 <- marker_list
        run_gsea(bin_mat, v1, n_perm = 1000, per_cell = TRUE)
      }
    )
    res <- do.call(cbind, out)
    n <- ncol(res)
    res2 <- res[, 3 * c(1:(ncol(res) / 3)), drop = FALSE]
    rownames(res2) <- rownames(res)
    colnames(res2) <- colnames(marker_m)
    res <- res2
  }

  if (output_high == TRUE) {
    if (metric == "hyper") {
      res <- -log10(res)
    } else if (metric == "spearman") {
      res <- -res + max(res)
    }
  }

  res
}
