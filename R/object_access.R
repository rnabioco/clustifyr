#' An example Seurat object
#' 
#' @return a SingleCellExperiment object populated with data
#' from the [pbmc_matrix_small] scRNA-seq dataset, additionally
#' annotated with cluster assignments.
#' 
#' @importFrom SeuratObject CreateSeuratObject CreateDimReducObject VariableFeatures
#' @export
so_pbmc <- function() {
  x <- pbmc_example_data()
  so <- SeuratObject::CreateSeuratObject(x$mat,
                                         meta.data = x$metadata) 
  umap_dr <- SeuratObject::CreateDimReducObject(embeddings = x$umap,
                                                key = "umap_",
                                                assay = "RNA")
  if(is_seurat_v5()) {
    so <- SeuratObject::SetAssayData(so,
      "data", 
      SeuratObject::LayerData(so, layer = "counts")
    )
  } else {
    so <- SeuratObject::SetAssayData(so,
      "data", 
      SeuratObject::GetAssayData(so, slot = "counts")
    )
  }
  so[["umap"]] <- umap_dr
  SeuratObject::VariableFeatures(so) <- x$vargenes
  so
}
  
#' An example SingleCellExperiment object
#' 
#' @return a SingleCellExperiment object populated with data
#' from the [pbmc_matrix_small] scRNA-seq dataset, additionally
#' annotated with cluster assignments.
#' 
#' @export
sce_pbmc <- function() {
  x <- pbmc_example_data()
  md <- x$metadata[, c(1:5, 7)]
  # rename to more sce-like names
  colnames(md) <- c("cell_source", 
                    "sum",
                    "detected",
                    "subsets_Mito_percent",
                    "clusters",
                    "cell_type")
  SingleCellExperiment::SingleCellExperiment(list(counts = x$mat,
                                                  logcounts = x$mat), 
                                             colData = md,
                                             reducedDims = list(
                                               UMAP = x$umap
                                             ))
}

pbmc_example_data <- function() {
  mat <- clustifyr::pbmc_matrix_small
  md <- clustifyr::pbmc_meta
  umap_cols <- c("UMAP_1", "UMAP_2")
  umap <- as.matrix(md[, umap_cols])
  md <- md[, setdiff(colnames(md), umap_cols)]
  vargenes <- clustifyr::pbmc_vargenes
  
  list(mat = mat,
       metadata = md,
       umap = umap,
       vargenes = vargenes)
}

#' Function to access object data
#' @return expression matrix, with genes as row names,
#' and cell types as column names
#' @export
object_data <- function(object, ...) {
    UseMethod("object_data", object)
}

#' @rdname object_data
#' @param object SingleCellExperiment or Seurat object 
#' @param slot data to access
#' @param n_genes number of genes limit for Seurat variable genes, by default 1000,
#'   set to 0 to use all variable genes (generally not recommended)
#' @param ... additional arguments
#' @examples
#' so <- so_pbmc()
#' mat <- object_data(
#'     object = so,
#'     slot = "data"
#' )
#' mat[1:3, 1:3]
#' @export
object_data.Seurat <- function(object,
    slot,
    n_genes = 1000,
    ...) {
    if (slot == "data") {
        temp <- get_seurat_matrix(object, ...)
        return(temp)
    } else if (slot == "meta.data") {
        return(object@meta.data)
    } else if (slot == "var.genes") {
        vars <- SeuratObject::VariableFeatures(object)
       
        if (is.null(vars) || length(vars) <= 1) {
            message("variable genes not found, please manually specify with query_genes argument")
        }
        if ((length(vars) > n_genes) & (n_genes > 0)) {
            vars <- vars[seq_len(n_genes)]
        }
        
        return(vars)
    } else {
      stop(slot, " access method not implemented")
    }
}

#' @importFrom utils packageVersion
is_seurat_v5 <- function() {
  utils::packageVersion("SeuratObject") >= '5.0.0'
}

extract_v5_matrix <- function(x, ...) {
  ob_layers <- SeuratObject::Layers(x)
  if("data" %in% ob_layers) {
    res <- SeuratObject::LayerData(x, layer = "data", ...)
  } else if ("counts" %in% ob_layers) {
    message("Unable to find 'data' layer, using 'count' layer instead")
    res <- SeuratObject::LayerData(x, layer = "counts", ...)
  } else {
    da <- DefaultAssay(x)
    stop("\nUnable to find data or count layer in ", da, " Assay of SeuratObject\n",
         "Extracting data from V5 objects with multiple count\n",
         "or data layers is not supported")
  }
  res
}

extract_v4_matrix <- function(x) {
  res <- SeuratObject::GetAssayData(x, layer = "data")
  
  if(length(res) == 0) {
    message("Unable to find 'data' slot, using 'count' slot instead")
    res <- SeuratObject::GetAssayData(x, layer = "count")
  }
  
  res
}

get_seurat_matrix <- function(x, warn = TRUE) {
  
  ob_assay <- SeuratObject::DefaultAssay(x)
  if(warn && ob_assay != "RNA") {
    warning("Default assay of input Seurat object is ", ob_assay, "\n",
            "Data will be used from this assay rather than RNA")
  }
  
  if(is_seurat_v5()) {
    res <- extract_v5_matrix(x)
  } else {
    res <- extract_v4_matrix(x)
  }
  res
}

#' @rdname object_data
#' @param object object after tsne or umap projections
#'  and clustering
#' @param slot data to access
#' @param ... additional arguments
#' @importFrom SingleCellExperiment logcounts colData
#' @examples
#' sce <- sce_pbmc()
#' mat <- object_data(
#'     object = sce,
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
    } else {
        stop(slot, " access method not implemented")
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
#' so <- so_pbmc()
#' obj <- write_meta(
#'     object = so,
#'     meta = seurat_meta(so)
#' )
#' @export
write_meta.Seurat <- function(object,
    meta,
    ...) {
    object_new <- object
    object_new@meta.data <- meta
    object_new
}

#' @rdname write_meta
#' @param object object after tsne or umap projections
#'  and clustering
#' @param meta new metadata dataframe
#' @param ... additional arguments
#' @importFrom SingleCellExperiment colData
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData<-
#' @examples
#' sce <- sce_pbmc()
#' obj <- write_meta(
#'     object = sce,
#'     meta = object_data(sce, "meta.data")
#' )
#' @export
write_meta.SingleCellExperiment <- function(object,
    meta,
    ...) {
    colData(object) <- S4Vectors::DataFrame(meta)
    object
}

#' Function to convert labelled seurat object to avg expression matrix
#' @return reference expression matrix, with genes as row names,
#' and cell types as column names
#' @examples 
#' so <- so_pbmc()
#' ref <- seurat_ref(so, cluster_col = "seurat_clusters")
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
            og_assay <- SeuratObject::DefaultAssay(seurat_object)
            assay_name <- setdiff(assay_name, og_assay)
            temp_mat <- temp_mat[0, ]
            for (element in assay_name) {
                SeuratObject::DefaultAssay(seurat_object) <- element
                temp_mat2 <- object_data(seurat_object, "data", warn = FALSE)
                temp_mat <- rbind(temp_mat, as.matrix(temp_mat2))
            }
            SeuratObject::DefaultAssay(seurat_object) <- og_assay
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
#' so <- so_pbmc()
#' m <- seurat_meta(so)
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
#' so <- so_pbmc()
#' object_ref(
#'     so,
#'     cluster_col = "seurat_clusters"
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
