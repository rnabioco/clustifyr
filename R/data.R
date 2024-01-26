#' Matrix of single-cell RNA-seq PBMCs.
#'
#' Count matrix of 3k pbmcs from Seurat3 tutorial, with only var.features
#'
#' @format A sparseMatrix with genes as rows and cells as columns.
#'
#' @family data
#'
#' @source \url{https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html}
"pbmc_matrix_small"

#' Meta-data for single-cell RNA-seq PBMCs.
#'
#' Metadata, including umap, of 3k pbmcs from Seurat3 tutorial
#'
#' @family data
#' @source `[pbmc_matrix]` processed by Seurat
"pbmc_meta"

#' Marker genes identified by Seurat from single-cell RNA-seq PBMCs.
#'
#' Dataframe of markers from Seurat FindAllMarkers function
#'
#' @family data
#' @source `[pbmc_matrix]` processed by Seurat
"pbmc_markers"

#' Marker genes identified by M3Drop from single-cell RNA-seq PBMCs.
#'
#' Selected features of 3k pbmcs from Seurat3 tutorial
#'
#' @format A data frame with 3 variables:
#'
#' @family data
#' @source `[pbmc_matrix]` processed by `[M3Drop]`
"pbmc_markers_M3Drop"

#' Variable genes identified by Seurat from single-cell RNA-seq PBMCs.
#'
#' Top 2000 variable genes from 3k pbmcs from Seurat3 tutorial
#'
#' @family data
#' @source `[pbmc_matrix]` processed by Seurat
"pbmc_vargenes"

#' reference matrix from seurat citeseq CBMC tutorial
#'
#' @family data
#' @source \url{https://satijalab.org/seurat/v3.0/multimodal_vignette.html#identify-differentially-expressed-proteins-between-clusters}
"cbmc_ref"

#' reference marker matrix from seurat citeseq CBMC tutorial
#'
#' @family data
#' @source \url{https://satijalab.org/seurat/v3.0/multimodal_vignette.html#identify-differentially-expressed-proteins-between-clusters}
"cbmc_m"

#' table of references stored in clustifyrdata
#'
#' @family data
#' @source  various packages
"downrefs"

#' Vector of human genes for 10x cellranger pipeline
#'
#' @family data
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest}
"human_genes_10x"

#' Vector of mouse genes for 10x cellranger pipeline
#'
#' @family data
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest}
"mouse_genes_10x"
