#' Wrapper to integate into Seurat workflow
#'
#' @param seurat_object seruat_object after tsne projections and clustering
#' @param bulk_mat bulk expression matrix, or matrixized ranked list
#' @param query_genes A vector of genes of interest, or by default loads var.genes from seurat object, a number can be given to subset var.genes
#' @param per_cell if true run per cell, otherwise per cluster.
#' @param num_perm number of permutations, set to 0 by default
#' @param return_full Return full results includeing scores instead
#'  of correlation coefficient only.
#' @param compute_method method(s) for computing similarity scores
#' @param output either output a new seurat object, or just a dataframe of calls
#' @param ... additional arguments to pass to compute_method function
#'
#' @export
clustify_seurat <- function(seurat_object,
                            bulk_mat,
                            query_genes = "var.genes",
                            per_cell = F,
                            num_perm = 0,
                            cluster_col = "cluster",
                            compute_method = "pearson",
                            output = "object",
                            carry_cor = FALSE) {
  expr_mat <- seurat_object@data

  metadata <- seurat_object@meta.data %>%
    tibble::rownames_to_column("rn") %>%
    tibble::as_tibble()
  tsne <- as.data.frame(seurat_object@dr$tsne@cell.embeddings) %>%
    tibble::rownames_to_column("rn") %>%
    tibble::as_tibble()

  metadata <- inner_join(metadata,
                         tsne %>% select(cell, cluster_col),
                         by = "rn")

  # if query_genes == NULL, use var.genes saved in seurat object
  if (query_genes == "var.genes") {
    query_genes <- seurat_object@var.genes
  } else if (typeof(query_genes) == "double") {
    query_genes <- (seurat_object@hvg.info %>% rownames_to_column() %>% arrange(desc(gene.dispersion.scaled)) %>% pull(rowname))[1:query_genes]
  }

  # if per_cell, cluster_col should default to "rn"
  if ((per_cell == TRUE) & (dim(unique(metadata[, cluster_col]))[1] != dim(metadata)[1])) {
    cluster_orig <- cluster_col
    cluster_col <- "rn"
  }

  res <- clustify(expr_mat,
    metadata,
    bulk_mat,
    query_genes,
    per_cell = per_cell,
    num_perm = num_perm,
    cluster_col = cluster_col,
    compute_method = compute_method,
    return_full = FALSE
  )

  df_call <- as.data.frame(t(apply(res, 1, function(x) x - max(x))))
  df_call[df_call == 0] <- "1"
  df_call[df_call != "1"] <- "0"

  calls <- c()
  for (name in rownames(df_call)) {
    calls <- c(calls, get_best_str(name, df_call, res, carry_cor))
  }
  calls_df <- data.frame("tempid" = rownames(df_call),
                         "call" = calls,
                         stringsAsFactors = FALSE)
  calls_df

  metadata <- seurat_object@meta.data
  calls <- left_join(metadata, calls_df,
                     by = structure(names = cluster_col, "
                                    tempid")) %>%
    select(call)
  rownames(calls) <- rownames(seurat_object@meta.data)

  if (output == "object") {
    seurat_object@meta.data$call <- calls
    return(seurat_object)
  } else if (output == "df") {
    return(calls)
  }
}

#' Function to convert labelled seurat object to avg expression matrix
#'
#' @param seurat_object seruat_object after tsne projections and clustering
#' @param cluster_col column name where classified cluster names are stored in seurat meta data, cannot be "rn"
#' @param var.genes_only whether to keep only var.genes in the final matrix output
#'
#' @export
use_seurat_comp <- function(seurat_object,
                            cluster_col = "classified",
                            var.genes_only = FALSE) {
  temp_mat <- average_clusters(seurat_object@data,
    seurat_object@meta.data,
    log_scale = TRUE,
    cluster_col = cluster_col
  )
  if (var.genes_only == TRUE) {
    temp_mat <- temp_mat[seurat_object@var.genes, ]
  }
  temp_mat
}
