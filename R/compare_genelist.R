#' Binarize scRNA seq data
#'
#' @param expression_matrix single-cell expression matrix
#' @param n number of top expressing genes to keep
#' @export
binarize_expr <- function(expression_matrix,
                          n = 1000){
  expr_df <- as.data.frame(expression_matrix)
  df_temp <- apply(expr_df, 2, function(x) x - x[order(x, decreasing=TRUE)[n + 1]])
  df_temp[df_temp > 0] = 1
  df_temp[df_temp < 0] = 0
  df_temp
}

#' convert candidate genes list into matrix
#'
#' @param marker_df dataframe of candidate genes
#' @param ranked unranked gene list feeds into hyperp, ranked gene list feeds into regular corr_coef
#' @param weight ranked genes are tranformed into pseudo expression with added weight
#' @param marker_df dataframe of candidate genes
#' @export
matrixize_markers <- function(marker_df, ranked = FALSE, weight = 0){
  # takes marker in dataframe form
  # equal number of marker genes per known cluster
  cut_num <- min((marker_df %>% group_by(cluster) %>% summarize(n =n()))$n)
  marker_temp <- marker_df %>% select(gene, cluster) %>% group_by(cluster) %>% slice(1:cut_num)
  if(ranked == TRUE){
    marker_temp <- marker_temp %>% mutate(n = (cut_num + weight) : (1 + weight))
    marker_temp2 <- as.data.frame(tidyr::spread(marker_temp, key = "cluster", value = n) %>% replace(is.na(.), 0))
    rownames(marker_temp2) <- marker_temp2$gene
    marker_temp2 <- marker_temp2 %>% select(-gene)
  } else {
    marker_temp <- marker_temp %>% mutate(n = 1:cut_num)
    marker_temp2 <- as.data.frame(tidyr::spread(marker_temp, key = "cluster", value = "gene") %>% select(-n))
  }
  marker_temp2
}

#' calculate adjusted p-values for hypergeometric test of gene lists
#'
#' @param bin_matrix binarized single-cell expression matrix
#' @param marker_m matrix or dataframe of candidate genes for each cluster
#' @param n number of genes in the genome
#' @export
hyperp <- function(bin_matrix, marker_m, n = 30000){
  # "expressed" genes per single cell data cluster
  out <- lapply(colnames(bin_matrix),
                function(x){
                  per_col <- lapply(colnames(marker_m),
                                    function(y){
                                      marker_list <- unlist(marker_m[,y],use.names = FALSE)
                                      bin_temp <- bin_matrix[,x][bin_matrix[,x] == 1]
                                      list_top <- names(bin_temp)

                                      t <- length(intersect(list_top, marker_list))
                                      a <- max(length(list_top),length(marker_list))
                                      b <- min(length(list_top),length(marker_list))
                                      sum(dhyper(t:b, a, n - a, b))})
                  do.call(cbind, as.list(p.adjust(per_col)))
                })

  res <- do.call(rbind, out)
  rownames(res) <- colnames(bin_matrix)
  colnames(res) <- colnames(marker_m)
  res
}
