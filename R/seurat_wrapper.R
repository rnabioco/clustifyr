#' Function to convert labelled seurat object to avg expression matrix
#'
#' @param seurat_object seurat_object after tsne projections and clustering
#' @param cluster_col column name where classified cluster names are stored in seurat meta data, cannot be "rn"
#' @param var.genes_only whether to keep only var.genes in the final matrix output, could also look up genes used for PCA
#' @param assay_name any additional assay data, such as ADT, to include. If more than 1, pass a vector of names
#' @param method whether to take mean (default) or median
#' @param seurat3 if using newest version

#'
#' @export
use_seurat_comp <- function(seurat_object,
                            cluster_col = "classified",
                            var.genes_only = FALSE,
                            assay_name = NULL,
                            method = "mean") {
  if (class(seurat_object) == "Seurat") {
    return(use_seurat3_comp(seurat_object,
                            cluster_col,
                            var.genes_only,
                            assay_name,
                            method
                            ))
  }

  temp_mat <- seurat_object@data
  if (var.genes_only == TRUE) {
    temp_mat <- temp_mat[seurat_object@var.genes, ]
  } else if (var.genes_only == "PCA") {
    temp_mat <- temp_mat[rownames(seurat_object@dr$pca@gene.loadings), ]
  }

  if (!is.null(assay_name)) {
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

  temp_res <- average_clusters(temp_mat,
    seurat_object@meta.data,
    log_scale = TRUE,
    cluster_col = cluster_col,
    method = method
  )

  temp_res
}

#' Function to convert labelled seurat3 object to avg expression matrix
#'
#' @param seurat_object seurat_object after tsne projections and clustering
#' @param cluster_col column name where classified cluster names are stored in seurat meta data, cannot be "rn"
#' @param var.genes_only whether to keep only var.genes in the final matrix output, could also look up genes used for PCA
#' @param assay_name any additional assay data, such as ADT, to include. If more than 1, pass a vector of names
#' @param method whether to take mean (default) or median
#' @param seurat3 if using newest version

#'
#' @export
use_seurat3_comp <- function(seurat_object,
                            cluster_col = "classified",
                            var.genes_only = FALSE,
                            assay_name = NULL,
                            method = "mean") {
  temp_mat <- seurat_object@assays$RNA@data

  if (var.genes_only == TRUE) {
    temp_mat <- temp_mat[seurat_object@assays$RNA@var.features, ]
  } else if (var.genes_only == "PCA") {
    temp_mat <- temp_mat[rownames(seurat_object@reductions$pca@feature.loadings), ]
  }

  if (!is.null(assay_name)) {
    if (!is.vector(assay_name)) {
      temp_mat2 <- seurat_object@assays[[assay_name]]@counts
      temp_mat <- rbind(temp_mat, as.matrix(temp_mat2))
    } else {
      for (element in assay_name) {
        temp_mat2 <- seurat_object@assays[[element]]@counts
        temp_mat <- rbind(temp_mat, as.matrix(temp_mat2))
      }
    }
  }

  temp_res <- average_clusters(temp_mat,
                               seurat_object@meta.data,
                               log_scale = TRUE,
                               cluster_col = cluster_col,
                               method = method
  )

  temp_res
}

#' Function to convert labelled seurat object to fully prepared metadata
#'
#' @param seurat_object seurat_object after tsne projections and clustering
#' @param dr dimension reduction method
#' @param seurat3 if using newest version
#' @export
use_seurat_meta <- function(seurat_object,
                            dr = "tsne",
                            seurat3 = F) {
  if (class(seurat_object) == "Seurat") {
    seurat3 = T
  }
  if (seurat3 == F) {
    temp_dr <- as.data.frame(seurat_object@dr[[dr]]@cell.embeddings)
  } else {
    temp_dr <- as.data.frame(seurat_object@reductions[[dr]]@cell.embeddings)
  }
  temp_dr <- tibble::rownames_to_column(temp_dr, "rn")
  temp_meta <- tibble::rownames_to_column(seurat_object@meta.data, "rn")
  temp <- dplyr::left_join(temp_meta, temp_dr, by = "rn")
  tibble::column_to_rownames(temp, "rn")
}
