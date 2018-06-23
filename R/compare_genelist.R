#' Binarize scRNA seq data
#'
#' @param expr_mat single-cell expression matrix
#' @param n number of top expressing genes to keep
#'
#' @export
binarize_expr <- function(expr_mat,
                          n = 1000){
  expr_df <- as.data.frame(as.matrix(expr_mat))
  df_temp <- apply(expr_df, 2, function(x) x - x[order(x, decreasing=TRUE)[n + 1]])
  df_temp[df_temp > 0] = 1
  df_temp[df_temp < 0] = 0
  df_temp
}

#' convert candidate genes list into matrix
#'
#' @param marker_df dataframe of candidate genes
#' @param ranked unranked gene list feeds into hyperp, the ranked
#' gene list feeds into regular corr_coef
#' @param weight ranked genes are tranformed into pseudo expression with
#' added weight
#' @param marker_df dataframe of candidate genes
#'
#' @export
matrixize_markers <- function(marker_df, ranked = FALSE, weight = 0){
  # takes marker in dataframe form
  # equal number of marker genes per known cluster
  marker_df <- as_tibble(marker_df)
  cut_num <- min((marker_df %>% group_by(cluster) %>% summarize(n =n()))$n)
  marker_temp <- marker_df %>% dplyr::select(gene, cluster) %>% group_by(cluster) %>% dplyr::slice(1:cut_num)
  if(ranked == TRUE){
    marker_temp <- marker_temp %>% mutate(n = (cut_num + weight) : (1 + weight))
    marker_temp2 <- as.data.frame(tidyr::spread(marker_temp, key = "cluster", value = n) %>% replace(is.na(.), 0))
    rownames(marker_temp2) <- marker_temp2$gene
    marker_temp2 <- marker_temp2 %>% dplyr::select(-gene)
  } else {
    marker_temp <- marker_temp %>% mutate(n = 1:cut_num)
    marker_temp2 <- as.data.frame(tidyr::spread(marker_temp, key = "cluster", value = "gene") %>% dplyr::select(-n))
  }
  marker_temp2
}

#' calculate adjusted p-values for hypergeometric test of gene lists
#' or jaccard index
#'
#' @param bin_mat binarized single-cell expression matrix
#' @param marker_m matrix or dataframe of candidate genes for each cluster
#' @param n number of genes in the genome
#' @param metric adjusted p-value for hypergeometric test, or jaccard index
#'
#' @export
compare_lists <- function(bin_mat, marker_m, n = 30000, metric = "hyper"){
  # "expressed" genes per single cell data cluster
  if (metric == "hyper"){
    out <- lapply(colnames(bin_mat),
                function(x){
                  per_col <- lapply(colnames(marker_m),
                                    function(y){
                                      marker_list <- unlist(marker_m[,y],use.names = FALSE)
                                      bin_temp <- bin_mat[,x][bin_mat[,x] == 1]
                                      list_top <- names(bin_temp)

                                      t <- length(intersect(list_top, marker_list))
                                      a <- max(length(list_top),length(marker_list))
                                      b <- min(length(list_top),length(marker_list))
                                      sum(dhyper(t:b, a, n - a, b))})
                  do.call(cbind, as.list(p.adjust(per_col)))
                })
  }

  if (metric == "jaccard"){
    out <- lapply(colnames(bin_mat),
                  function(x){
                    per_col <- lapply(colnames(marker_m),
                                      function(y){
                                        marker_list <- unlist(marker_m[,y],use.names = FALSE)
                                        bin_temp <- bin_mat[,x][bin_mat[,x] == 1]
                                        list_top <- names(bin_temp)

                                        I <- length(intersect(list_top, marker_list))
                                        I/(length(list_top) + length(marker_list) - I)
                                        })
                    do.call(cbind, per_col)
                  })
  }
  res <- do.call(rbind, out)
  rownames(res) <- colnames(bin_mat)
  colnames(res) <- colnames(marker_m)
  res
}
