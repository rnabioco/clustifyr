#' Function to convert labelled seurat object to avg expression matrix
#'
#' @param seurat_object seurat_object after tsne projections and clustering
#' @param cluster_col column name where classified cluster names are stored in seurat meta data, cannot be "rn"
#' @param var.genes_only whether to keep only var.genes in the final matrix output, could also look up genes used for PCA
#' @param assay_name any additional assay data, such as ADT, to include. If more than 1, pass a vector of names
#'
#' @export
use_seurat_comp <- function(seurat_object,
                            cluster_col = "classified",
                            var.genes_only = FALSE,
                            assay_name = NULL) {
  temp_mat <- seurat_object@data

  if (!is.null(assay_name)){
    if (!is.vector(assay_name)) {
      temp_mat2 <- seurat_object@assay[[assay_name]]@raw.data
      temp_mat <- rbind(temp_mat, as.matrix(temp_mat2))
    } else {
      for (element in assay_name) {
        temp_mat2 <- seurat_object@assay[[element]]@raw.data
        temp_mat <- rbind(temp_mat, as.matrix(temp_mat2))
      }
    }
  }

  if (var.genes_only == TRUE) {
    temp_mat <- temp_mat[seurat_object@var.genes, ]
  } else if (var.genes_only == "PCA") {
    temp_mat <- temp_mat[rownames(seurat_object@dr$pca@gene.loadings), ]
  }

  temp_res <- average_clusters(temp_mat,
                               seurat_object@meta.data,
                               log_scale = TRUE,
                               cluster_col = cluster_col)

  temp_res
}

#' Function to convert labelled seurat object to fully prepared metadata
#'
#' @param seurat_object seurat_object after tsne projections and clustering
#' @param dr dimension reduction method
#' @export
use_seurat_meta <- function(seurat_object,
                            dr = "tsne") {
  temp_dr <- as.data.frame(seurat_object@dr[[dr]]@cell.embeddings)
  temp_dr <- rownames_to_column(temp_dr, "rn")
  temp_meta <- rownames_to_column(seurat_object@meta.data, "rn")
  temp <- left_join(temp_meta, temp_dr, by = "rn")
  column_to_rownames(temp, "rn")
}
