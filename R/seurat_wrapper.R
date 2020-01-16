#' Function to convert labelled seurat object to avg expression matrix
#' @return reference expression matrix, with genes as row names,
#' and cell types as column names
#' @export
seurat_ref <- function(seurat_object, ...) {
    UseMethod("seurat_ref", seurat_object)
}

#' @rdname seurat_ref
#' @param seurat_object seurat_object after tsne or umap projections
#'  and clustering
#' @param cluster_col column name where classified cluster names
#'  are stored in  seurat meta data, cannot be "rn"
#' @param var_genes_only whether to keep only var_genes in the final
#' matrix output, could also look up genes used for PCA
#' @param assay_name any additional assay data, such as ADT, to include.
#'  If more than 1, pass a vector of names
#' @param method whether to take mean (default) or median
#' @param subclusterpower whether to get multiple averages per
#'  original cluster
#' @param if_log input data is natural log,
#' averaging will be done on unlogged data
#' @param ... additional arguments
#' @examples
#' ref <- seurat_ref(
#'     seurat_object = s_small,
#'     cluster_col = "res.1",
#'     var_genes_only = TRUE
#' )
#' ref[1:3, 1:3]
#' @export
seurat_ref.seurat <- function(seurat_object,
                              cluster_col = "classified",
                              var_genes_only = FALSE,
                              assay_name = NULL,
                              method = "mean",
                              subclusterpower = 0,
                              if_log = TRUE,
                              ...) {
    temp_mat <- seurat_object@data
    if (is.logical(var_genes_only) && var_genes_only) {
        temp_mat <- temp_mat[seurat_object@var.genes, ]
    } else if (var_genes_only == "PCA") {
        temp_mat <- temp_mat[rownames(seurat_object@dr$pca@gene.loadings), ]
    }

    if (!is.null(assay_name)) {
        temp_mat <- temp_mat[0, ]
        for (element in assay_name) {
            temp_mat2 <- seurat_object@assay[[element]]@raw.data
            temp_mat <- rbind(temp_mat, as.matrix(temp_mat2))
        }
    }

    temp_res <- average_clusters(
        temp_mat,
        seurat_object@meta.data,
        cluster_col = cluster_col,
        method = method,
        subclusterpower = subclusterpower,
        if_log = if_log
    )

    temp_res
}

#' @rdname seurat_ref
#' @export
seurat_ref.Seurat <- function(seurat_object,
                              cluster_col = "classified",
                              var_genes_only = FALSE,
                              assay_name = NULL,
                              method = "mean",
                              subclusterpower = 0,
                              if_log = TRUE,
                              ...) {
    if (is(seurat_object, "Seurat")) {
        temp_mat <- seurat_object@assays$RNA@data

        if (is.logical(var_genes_only) && var_genes_only) {
            temp_mat <- temp_mat[seurat_object@assays$RNA@var.features, ]
        } else if (var_genes_only == "PCA") {
            temp_mat <-
              temp_mat[rownames(seurat_object@reductions$pca@feature.loadings),
                       ]
        }

        if (!is.null(assay_name)) {
            temp_mat <- temp_mat[0, ]
            for (element in assay_name) {
                temp_mat2 <- seurat_object@assays[[element]]@counts
                temp_mat <- rbind(temp_mat, as.matrix(temp_mat2))
            }
        }
    } else {
        stop("warning, not seurat3 object")
    }

    temp_res <- average_clusters(
        temp_mat,
        seurat_object@meta.data,
        cluster_col = cluster_col,
        method = method,
        subclusterpower = subclusterpower,
        if_log = if_log
    )

    temp_res
}

#' Function to convert labelled seurat object to fully prepared metadata
#' @return dataframe of metadata, including dimension reduction plotting info
#' @examples
#' \dontrun{
#' seurat_meta(s_small)
#' }
#' @export
seurat_meta <- function(seurat_object, ...) {
    UseMethod("seurat_meta", seurat_object)
}

#' @rdname seurat_meta
#' @param seurat_object seurat_object after tsne or
#'  umap projections and clustering
#' @param dr dimension reduction method
#' @param ... additional arguments
#' @export
seurat_meta.seurat <- function(seurat_object,
                               dr = "umap",
                               ...) {
    dr2 <- dr
    temp_dr <-
        tryCatch(
            as.data.frame(seurat_object@dr[[dr2]]@cell.embeddings),
            error = function(e) {
                message("cannot find dr info")
                return(NA)
            }
        )
    if (!is.data.frame(temp_dr)) {
        return(seurat_object@meta.data)
    } else {
        temp_dr <- tibble::rownames_to_column(temp_dr, "rn")
        temp_meta <-
            tibble::rownames_to_column(seurat_object@meta.data, "rn")
        temp <- dplyr::left_join(temp_meta, temp_dr, by = "rn")
        return(tibble::column_to_rownames(temp, "rn"))
    }
}

#' @rdname seurat_meta
#' @export
seurat_meta.Seurat <- function(seurat_object,
                               dr = "umap",
                               ...) {
    dr2 <- dr

    mdata <- seurat_object@meta.data
    temp_col_id <- get_unique_column(mdata, "rn")

    temp_dr <-
        tryCatch(
            as.data.frame(seurat_object@reductions[[dr2]]@cell.embeddings),
            error = function(e) {
                message("cannot find dr info")
                return(NA)
            }
        )
    if (!is.data.frame(temp_dr)) {
        return(mdata)
    } else {
        temp_dr <- tibble::rownames_to_column(temp_dr, temp_col_id)
        temp_meta <- tibble::rownames_to_column(mdata, temp_col_id)
        temp <- dplyr::left_join(temp_meta, temp_dr, by = temp_col_id)
        return(tibble::column_to_rownames(temp, temp_col_id))
    }
}

#' Function to convert labelled object to avg expression matrix
#'
#' @param input object after tsne or umap projections and clustering
#' @param cluster_col column name where classified cluster names
#' are stored in seurat meta data, cannot be "rn"
#' @param var_genes_only whether to keep only var.genes in the
#' final matrix output, could also look up genes used for PCA
#' @param assay_name any additional assay data, such as ADT, to include.
#'  If more than 1, pass a vector of names
#' @param method whether to take mean (default) or median
#' @param lookuptable if not supplied, will look
#' in built-in table for object parsing
#' @param if_log input data is natural log,
#' averaging will be done on unlogged data
#' @return reference expression matrix, with genes as row names,
#'  and cell types as column names
#' @examples
#' object_ref(
#'     s_small3,
#'     cluster_col = "RNA_snn_res.1"
#' )
#' @export
object_ref <- function(input,
                       cluster_col = NULL,
                       var_genes_only = FALSE,
                       assay_name = NULL,
                       method = "mean",
                       lookuptable = NULL,
                       if_log = TRUE) {
    if (!is(input, "seurat")) {
        input_original <- input
        temp <- parse_loc_object(
            input,
            type = class(input),
            expr_loc = NULL,
            meta_loc = NULL,
            var_loc = NULL,
            cluster_col = cluster_col,
            lookuptable = lookuptable
        )
        if (!(is.null(temp[["expr"]]))) {
            message("recognized object type - ", class(input))
        }
        input <- temp[["expr"]]
        metadata <- temp[["meta"]]
        query_genes <- temp[["var"]]
        if (is.null(cluster_col)) {
            cluster_col <- temp[["col"]]
        }
    }

    temp_mat <- input
    if (is.logical(var_genes_only) && var_genes_only) {
        temp_mat <- temp_mat[query_genes, ]
    }

    temp_res <- average_clusters(
        temp_mat,
        metadata,
        cluster_col = cluster_col,
        method = method,
        if_log = if_log
    )

    temp_res
}
