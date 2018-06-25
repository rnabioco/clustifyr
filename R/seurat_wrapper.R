#' Wrapper to integate into Seurat workflow
#'
#' @param seurat_object seruat_object after tsne projections and clustering
#' @param bulk_mat bulk expression matrix, or ranked list
#' @param query_gene_list A vector of genes of interest.
#' @param per_cell if true run per cell, otherwise per cluster.
#' @param num_perm number of permutations, set to 0 by default
#' @param return_full Return full results includeing scores instead
#'  of correlation coefficient only.
#' @param compute_method method(s) for computing similarity scores
#' @param output either output a new seurat object, or just a dataframe of calls
#' @param ... additional arguments to pass to compute_method function
#'
#' @export
#'
clustify_seurat <- function(seurat_object,
                            bulk_mat,
                            query_gene_list,
                            per_cell = F,
                            num_perm = 0,
                            cluster_col = "cluster",
                            compute_method = corr_coef,
                            output = "object"){
  expr_mat <- seurat_object@data

  metadata <- data.table::copy(seurat_object@meta.data)
  tsne <- data.table::copy(seurat_object@dr$tsne@cell.embeddings)
  metadata_tibble <- data.table::setDT(data.frame(tsne), keep.rownames = TRUE)
  metadata_tibble2 <- data.table::setDT(metadata, keep.rownames = TRUE)
  as_tibble(metadata_tibble2)
  as_tibble(metadata_tibble) -> metadata_tibble
  metadata <- inner_join(metadata_tibble, metadata_tibble2 %>% select(rn, cluster_col), by = "rn")

  res <- run_cor(expr_mat,
                 metadata,
                 bulk_mat,
                 query_gene_list,
                 per_cell = F,
                 num_perm = 0,
                 cluster_col = "cluster",
                 compute_method = corr_coef,
                 return_full = FALSE)

  df_temp <- as.data.frame(t(apply(res, 1, function(x) x - max(x))))
  df_call <- data.table::copy(df_temp)
  df_call[df_call==0]="1"
  df_call[df_call!="1"]="0"

  calls <- c()
  for(name in rownames(df_call)){
    calls <- c(calls, get_best_str(name, df_call, df_temp))
  }
  calls_df <- data.frame("cluster" = rownames(df_call), "call" = calls, stringsAsFactors=FALSE)
  calls_df
  metadata <- data.table::copy(seurat_object@meta.data)
  calls <- left_join(metadata, calls_df, by = "cluster") %>% select(call)
  rownames(calls) <- rownames(metadata)

  if (output == "object"){
  Seurat::AddMetaData(seurat_object, calls, "calls")
  } else if (output == "df"){
    calls
  }
}

#' Function to make call and attach score
#'
#' @param name name of row to query
#' @param best_mat binarized call matrix
#' @param cor_mat correlation matrix
#'
#' @export
#'
get_best_str <- function(name, best_mat, cor_mat) {
  if (sum(as.numeric(best_mat[name,])) > 0) {
    best.names <- colnames(best_mat)[which(best_mat[name,]==1)]
    best.cor <- round(cor_mat[name, which(best_mat[name,]==1)],2)
    for (i in 1:length(best.cor)) {
      if (i == 1) {
        str <- paste0(best.names[i], " (", best.cor[i], ") ")
      } else {
        str <- paste0(str, "; ", best.names[i], " (", best.cor[i], ") ")
      }
    }
  } else {
    str <- ""
  }
  return(str)
}

