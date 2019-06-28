#' Function to convert labelled seurat object to avg expression matrix
#' @return reference expression matrix, with genes as row names, and cell types as column names
#' @export
seurat_ref <- function(seurat_object, ...) {
  UseMethod("seurat_ref", seurat_object)
}

#' @rdname seurat_ref
#' @param seurat_object seurat_object after tsne or umap projections and clustering
#' @param cluster_col column name where classified cluster names are stored in seurat meta data, cannot be "rn"
#' @param var_genes_only whether to keep only var_genes in the final matrix output, could also look up genes used for PCA
#' @param assay_name any additional assay data, such as ADT, to include. If more than 1, pass a vector of names
#' @param method whether to take mean (default) or median
#'
#' @export
seurat_ref.seurat <- function(seurat_object,
                              cluster_col = "classified",
                              var_genes_only = FALSE,
                              assay_name = NULL,
                              method = "mean") {
  temp_mat <- seurat_object@data
  if (var_genes_only == TRUE) {
    temp_mat <- temp_mat[seurat_object@var.genes, ]
  } else if (var_genes_only == "PCA") {
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
    if_log = TRUE,
    cluster_col = cluster_col,
    method = method
  )

  temp_res
}

#' @rdname seurat_ref
#' @param seurat_object seurat_object after tsne or umap projections and clustering
#' @param cluster_col column name where classified cluster names are stored in seurat meta data, cannot be "rn"
#' @param var_genes_only whether to keep only var_genes in the final matrix output, could also look up genes used for PCA
#' @param assay_name any additional assay data, such as ADT, to include. If more than 1, pass a vector of names
#' @param method whether to take mean (default) or median
#'
#' @export
seurat_ref.Seurat <- function(seurat_object,
                              cluster_col = "classified",
                              var_genes_only = FALSE,
                              assay_name = NULL,
                              method = "mean") {
  if (class(seurat_object) == "Seurat") {
    temp_mat <- seurat_object@assays$RNA@data

    if (var_genes_only == TRUE) {
      temp_mat <- temp_mat[seurat_object@assays$RNA@var.features, ]
    } else if (var_genes_only == "PCA") {
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
  } else {
    stop("warning, not seurat3 object")
  }

  temp_res <- average_clusters(temp_mat,
    seurat_object@meta.data,
    if_log = TRUE,
    cluster_col = cluster_col,
    method = method
  )

  temp_res
}

#' Function to convert labelled seurat object to fully prepared metadata
#' @return dataframe of metadata, including dimension reduction plotting info
#' @export
seurat_meta <- function(seurat_object, ...) {
  UseMethod("seurat_meta", seurat_object)
}

#' @rdname seurat_meta
#' @param seurat_object seurat_object after tsne or umap projections and clustering
#' @param dr dimension reduction method
#' @export
seurat_meta.seurat <- function(seurat_object,
                               dr = "umap") {
  dr2 <- dr
  temp_dr <- as.data.frame(seurat_object@dr[[dr2]]@cell.embeddings)

  temp_dr <- tibble::rownames_to_column(temp_dr, "rn")
  temp_meta <- tibble::rownames_to_column(seurat_object@meta.data, "rn")
  temp <- dplyr::left_join(temp_meta, temp_dr, by = "rn")
  tibble::column_to_rownames(temp, "rn")
}

#' @rdname seurat_meta
#' @param seurat_object seurat_object after tsne or umap projections and clustering
#' @param dr dimension reduction method
#' @export
seurat_meta.Seurat <- function(seurat_object,
                               dr = "umap") {
  dr2 <- dr
  temp_dr <- as.data.frame(seurat_object@reductions[[dr2]]@cell.embeddings)

  temp_dr <- tibble::rownames_to_column(temp_dr, "rn")
  temp_meta <- tibble::rownames_to_column(seurat_object@meta.data, "rn")
  temp <- dplyr::left_join(temp_meta, temp_dr, by = "rn")
  tibble::column_to_rownames(temp, "rn")
}

#' Function to convert labelled object to avg expression matrix
#'
#' @param input object after tsne or umap projections and clustering
#' @param cluster_col column name where classified cluster names are stored in seurat meta data, cannot be "rn"
#' @param var_genes_only whether to keep only var.genes in the final matrix output, could also look up genes used for PCA
#' @param assay_name any additional assay data, such as ADT, to include. If more than 1, pass a vector of names
#' @param method whether to take mean (default) or median
#' @param lookuptable if not supplied, will look in built-in table for object parsing
#' @return reference expression matrix, with genes as row names, and cell types as column names

#' @export
object_ref <- function(input,
                       cluster_col = NULL,
                       var_genes_only = FALSE,
                       assay_name = NULL,
                       method = "mean",
                       lookuptable = NULL) {
  if (!(stringr::str_detect(class(input), "seurat"))) {
    input_original <- input
    temp <- parse_loc_object(input,
      type = class(input),
      expr_loc = NULL,
      meta_loc = NULL,
      var_loc = NULL,
      cluster_col = cluster_col,
      lookuptable = lookuptable
    )
    if (!(is.null(temp[["expr"]]))) {
      message(paste0("recognized object type - ", class(input)))
    }
    input <- temp[["expr"]]
    metadata <- temp[["meta"]]
    query_genes <- temp[["var"]]
    if (is.null(cluster_col)) {
      cluster_col <- temp[["col"]]
    }
  }

  temp_mat <- input
  if (var_genes_only == TRUE) {
    temp_mat <- temp_mat[query_genes, ]
  }

  temp_res <- average_clusters(temp_mat,
    metadata,
    if_log = TRUE,
    cluster_col = cluster_col,
    method = method
  )

  temp_res
}
