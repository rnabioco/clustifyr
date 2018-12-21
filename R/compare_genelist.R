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
  if(n < nrow(expr_mat)){
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
#'
#' @export
matrixize_markers <- function(marker_df,
                              ranked = FALSE,
                              n = NULL,
                              step_weight = 1,
                              background_weight = 0,
                              unique = FALSE,
                              labels = NULL) {
  # takes marker in dataframe form
  # equal number of marker genes per known cluster
  marker_df <- as_tibble(marker_df)

  # if "gene" not present in column names, assume df is a matrix to be converted to ranked
  if (!("gene" %in% colnames(marker_df))) {
    marker_df <- data.frame(lapply(marker_df, as.character), stringsAsFactors = FALSE) %>% tidyr::gather(factor_key = TRUE, key = "cluster", value = "gene")
  }

  if (unique == TRUE) {
    nonunique <- marker_df %>% group_by(gene) %>% summarize(n = n()) %>% filter(n > 1)
    marker_df <- anti_join(marker_df, nonunique, by = "gene")
  }

  cut_num <- min((marker_df %>% group_by(cluster) %>% summarize(n = n()))$n)

  if (!is.null(n)) {
    if (n < cut_num) {
      cut_num <- n
    }
  }

  marker_temp <- marker_df %>% dplyr::select(gene, cluster) %>% group_by(cluster) %>% dplyr::slice(1:cut_num)
  if (ranked == TRUE) {
    marker_temp <- marker_temp %>% mutate(n = seq(step_weight * cut_num, by = -step_weight, length.out = cut_num) + background_weight)
    marker_temp2 <- as.data.frame(tidyr::spread(marker_temp, key = "cluster", value = n) %>% replace(is.na(.), 0))
    rownames(marker_temp2) <- marker_temp2$gene
    marker_temp2 <- marker_temp2 %>% dplyr::select(-gene)
  } else {
    marker_temp <- marker_temp %>% mutate(n = 1:cut_num)
    marker_temp2 <- as.data.frame(tidyr::spread(marker_temp, key = "cluster", value = "gene") %>% dplyr::select(-n))
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
      ) %>%
        pull(classified)
    }
    colnames(marker_temp2) <- labels
  }

  marker_temp2
}

#' generate variable gene list from marker matrix
#'
#' @description Variable gene list is required for run_cor main function. This
#' function parses variables genes from a matrix input.
#'
#' @param marker_m matrix or dataframe of candidate genes for each cluster
#'
#' @export
get_vargenes <- function(marker_m) {
  unique(unlist(marker_m, use.names = FALSE))
}

#' calculate adjusted p-values for hypergeometric test of gene lists
#' or jaccard index
#'
#' @param bin_mat binarized single-cell expression matrix
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
  if (length(unique(bin_mat[, 1])) > 2) {
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
            v1 <- marker_list
            bin_temp <- as.matrix(bin_mat)[, x]
            bin_temp <- bin_temp[order(bin_temp, decreasing = TRUE)]
            list_top <- names(bin_temp)
            v2 <- list_top[list_top %in% v1]
            sum(sapply(seq_along(v1), function(i) abs(i - (which(v2 == v1[i])))))
          }
        )
        do.call(cbind, per_col)
      }
    )
  }

  res <- do.call(rbind, out)
  rownames(res) <- colnames(bin_mat)
  colnames(res) <- colnames(marker_m)

  if (output_high == TRUE) {
    if (metric == "hyper") {
      res <- -log10(res)
    } else if (metric == "spearman") {
      res <- -res
    }
  }

  res
}
