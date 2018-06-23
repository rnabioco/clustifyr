#' Average expression values per cluster
#' @param mat expression matrix
#' @param cluster_info data.frame with cells
#' @param log_scale input data is natural log, averaging will be done
#' @param cell_col column in cluster_info with cell ids
#' @param cluster_col column in cluster_info with cluster number
#' @export
average_clusters <- function(mat, cluster_info,
                             log_scale = T,
                             cell_col = "rn",
                             cluster_col = "cluster"){

  cluster_ids <- split(cluster_info[[cell_col]],
                       cluster_info[[cluster_col]])

  out <- lapply(cluster_ids,
                function(cell_ids){
                  if(!all(cell_ids %in% colnames(mat))){
                    stop("cell ids not found in input matrix")
                  }
                  if (log_scale) {
                    mat_data <- expm1(mat[, cell_ids])
                  } else {
                    mat_data <- mat[, cell_ids]
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
