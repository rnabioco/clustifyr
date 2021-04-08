#' Function to access object data
#' @return expression matrix, with genes as row names,
#' and cell types as column names
#' @export
object_data <- function(object, ...) {
    UseMethod("object_data", object)
}

#' @rdname object_data
#' @param object object after tsne or umap projections
#'  and clustering
#' @param slot data to access
#' @param n_genes number of genes limit for Seurat variable genes, by default 1000,
#'   set to 0 to use all variable genes (generally not recommended)
#' @param ... additional arguments
#' @examples
#' mat <- object_data(
#'     object = s_small3,
#'     slot = "data"
#' )
#' mat[1:3, 1:3]
#' @export
object_data.Seurat <- function(object,
    slot,
    n_genes = 1000,
    ...) {
    if (slot == "data") {
        temp <- tryCatch(object@assays$RNA@data,
                         error = function(e) {
                             message("detected spatial data, using raw counts")
                             object@assays$Spatial@counts
                             })
        return(temp)
    } else if (slot == "meta.data") {
        return(object@meta.data)
    } else if (slot == "var.genes") {
        vars <- tryCatch(object@assays$RNA@var.features,
                         error = function(e) {
                             object@assays$SCT@var.features 
                         })
        if (length(vars) == 0) {
            message("found variable genes in SCT slot")
            vars <- object@assays$SCT@var.features 
        } 
        if ((length(vars) > n_genes) & (n_genes > 0)) {
            vars <- vars[seq_len(n_genes)]
        }
        return(vars)
    } else if (slot == "pca") {
        return(object@reductions$pca@feature.loadings)
    }
}

#' @rdname object_data
#' @param object object after tsne or umap projections
#'  and clustering
#' @param slot data to access
#' @param ... additional arguments
#' @importFrom SingleCellExperiment logcounts colData
#' @examples
#' mat <- object_data(
#'     object = sce_small,
#'     slot = "data"
#' )
#' mat[1:3, 1:3]
#' @export
object_data.SingleCellExperiment <- function(object,
    slot,
    ...) {
    if (slot == "data") {
        return(SingleCellExperiment::logcounts(object))
    } else if (slot == "meta.data") {
        return(as.data.frame(SingleCellExperiment::colData(object)))
    }
}

#' Function to write metadata to object
#' @return object with newly inserted metadata columns
#' @export
write_meta <- function(object, ...) {
    UseMethod("write_meta", object)
}

#' @rdname write_meta
#' @param object object after tsne or umap projections
#'  and clustering
#' @param meta new metadata dataframe
#' @param ... additional arguments
#' @examples
#' obj <- write_meta(
#'     object = s_small3,
#'     meta = seurat_meta(s_small3)
#' )
#' @export
write_meta.Seurat <- function(object,
    meta,
    ...) {
    if ("Seurat" %in% loadedNamespaces()) {
        object_new <- object
        object_new@meta.data <- meta
        return(object_new)
    } else {
        message("Seurat not loaded")
    }
}

#' @rdname write_meta
#' @param object object after tsne or umap projections
#'  and clustering
#' @param meta new metadata dataframe
#' @param ... additional arguments
#' @importFrom SingleCellExperiment colData
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment `colData<-`
#' @examples
#' obj <- write_meta(
#'     object = sce_small,
#'     meta = object_data(sce_small, "meta.data")
#' )
#' @export
write_meta.SingleCellExperiment <- function(object,
    meta,
    ...) {
    if ("SingleCellExperiment" %in% loadedNamespaces()) {
        colData(object) <- S4Vectors::DataFrame(meta)
        return(object)
    } else {
        message("SingleCellExperiment not loaded")
    }
}

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
        temp_mat <- object_data(seurat_object, "data")

        if (is.logical(var_genes_only) && var_genes_only) {
            temp_mat <- temp_mat[object_data(seurat_object, "var.genes"), ]
        } else if (var_genes_only == "PCA") {
            temp_mat <-
                temp_mat[rownames(object_data(seurat_object, "pca")), ]
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
        object_data(seurat_object, "meta.data"),
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
#' seurat_meta(s_small3)
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
seurat_meta.Seurat <- function(seurat_object,
    dr = "umap",
    ...) {
    dr2 <- dr

    mdata <- object_data(seurat_object, "meta.data")
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
        if (tibble::has_rownames(temp)) {
            temp <- tibble::remove_rownames(temp)
        }
        return(tibble::column_to_rownames(temp, temp_col_id))
    }
}

#' Function to convert labelled object to avg expression matrix
#' @return reference expression matrix, with genes as row names,
#'  and cell types as column names
#' @export
object_ref <- function(input, ...) {
    UseMethod("object_ref", input)
}

#' @rdname object_ref
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
#' @param ... additional arguments
#' @examples
#' object_ref(
#'     s_small3,
#'     cluster_col = "RNA_snn_res.1"
#' )
#' @export
object_ref.default <- function(input,
    cluster_col = NULL,
    var_genes_only = FALSE,
    assay_name = NULL,
    method = "mean",
    lookuptable = NULL,
    if_log = TRUE,
    ...) {
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

#' @rdname object_ref
#' @export
object_ref.Seurat <- function(input,
    cluster_col = NULL,
    var_genes_only = FALSE,
    assay_name = NULL,
    method = "mean",
    lookuptable = NULL,
    if_log = TRUE,
    ...) {
    temp_mat <- object_data(input, "data")
    metadata <- object_data(input, "meta.data")
    query_genes <- object_data(input, "var.genes")
    if (is.null(cluster_col)) {
        message("please indicate metadata column containing cell identities")
    }

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

#' @rdname object_ref
#' @export
object_ref.SingleCellExperiment <- function(input,
    cluster_col = NULL,
    var_genes_only = FALSE,
    assay_name = NULL,
    method = "mean",
    lookuptable = NULL,
    if_log = TRUE,
    ...) {
    temp_mat <- object_data(input, "data")
    metadata <- object_data(input, "meta.data")
    if (is.null(cluster_col)) {
        message("please indicate metadata column containing cell identities")
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
