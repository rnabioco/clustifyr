#' Function to convert labelled seurat object to avg expression matrix
#'
#' @param seurat_object seurat_object after tsne projections and clustering
#' @param cluster_col column name where classified cluster names are stored in seurat meta data, cannot be "rn"
#' @param var.genes_only whether to keep only var.genes in the final matrix output
#' @param assay any additional assay data, such as ADT, to include
#'
#' @export
use_seurat_comp <- function(seurat_object,
                            cluster_col = "classified",
                            var.genes_only = FALSE,
                            assay = NULL) {
  temp_mat <- seurat_object@data

  if (!is.null(assay)){
    temp_mat2 <- seurat_object@assay[[assay]]@raw.data
    temp_mat <- rbind(temp_mat, as.matrix(temp_mat2))
  }

  if (var.genes_only == TRUE) {
    temp_mat <- temp_mat[seurat_object@var.genes, ]
  }


  temp_res <- average_clusters(temp_mat,
                               seurat_object@meta.data,
                               log_scale = TRUE,
                               cluster_col = cluster_col)

  temp_res
}
